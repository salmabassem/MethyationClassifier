
###################################################################################

Prepare_Training = function ( MethData, Epitypes, Query)  {
  
  rownames(MethData) <- MethData[, 1]
  drop_names <- intersect(colnames(MethData), c("sample.ID", "Epitype.Color2"))
  if (length(drop_names)) {
    MethData <- MethData[, setdiff(colnames(MethData), drop_names), drop = FALSE]
  }
  
  cpg <- intersect(colnames(MethData), colnames(Query))
  if (length(cpg) == 0L) stop("No overlapping CpGs between training and query.")
  MethData <- MethData[, cpg, drop = FALSE]
  
  rownames(Epitypes) <- Epitypes$sample.ID
  table(Epitypes$Subtype, Epitypes$Epitype)
  
  rownames(Epitypes) <- Epitypes$sample.ID
  Epitypes <- Epitypes["Subtype"]
  iPlexData_Epi <- transform(merge(MethData, Epitypes, by = 0), row.names = 1)
  iPlexData_Epi$Subtype <- as.factor(iPlexData_Epi$Subtype)
  
  df <- data.frame(iPlexData_Epi, check.names = FALSE)
  
  return(df)

} 


###################################################################################

nested_cv_calibrated_rf <- function(
    data, label_col = "Subtype",
    ntrees = 500, mtry = 6,
    outer_folds = 3, inner_folds = 3,
    seed = 123
){
  stopifnot(label_col %in% names(data))
  y <- data[[label_col]]
  feats <- setdiff(names(data), label_col)
  
  set.seed(seed)
  outer_idx <- caret::createFolds(y, k = outer_folds, list = TRUE, returnTrain = TRUE)
  
  misclassification_errors <- c()
  auc_values <- c()
  
  # helper for balanced downsampling
  make_sampsize <- function(yv) {
    tab <- table(yv)
    s <- rep(min(tab), length(tab)); names(s) <- names(tab); s
  }
  
  last_calibration_probs  <- NULL
  last_calibration_labels <- NULL
  
  for (i in seq_along(outer_idx)) {
    train_indices <- outer_idx[[i]]
    train_data <- data[train_indices, , drop = FALSE]
    test_data  <- data[-train_indices, , drop = FALSE]
    
    class_distribution <- table(train_data[[label_col]])
    sampsize <- rep(min(class_distribution), length(class_distribution)); names(sampsize) <- names(class_distribution)
    
    set.seed(seed + i)
    rf_model <- randomForest::randomForest(
      x = train_data[, feats, drop = FALSE],
      y = train_data[[label_col]],
      ntree = ntrees, mtry = mtry,
      sampsize = sampsize,
      importance = TRUE, proximity = TRUE, oob.prox = TRUE, keep.forest = TRUE
    )
    
    inner_folds_idx <- caret::createFolds(train_data[[label_col]], k = inner_folds)
    
    calibration_probs  <- NULL
    calibration_labels <- NULL
    
    for (j in seq_along(inner_folds_idx)) {
      inner_train <- train_data[ inner_folds_idx[[j]], , drop = FALSE]
      inner_test  <- train_data[-inner_folds_idx[[j]], , drop = FALSE]
      
      class_distribution <- table(inner_train[[label_col]])
      sampsize <- rep(min(class_distribution), length(class_distribution)); names(sampsize) <- names(class_distribution)
      
      inner_rf <- randomForest::randomForest(
        x = inner_train[, feats, drop = FALSE],
        y = inner_train[[label_col]],
        ntree = ntrees, mtry = mtry,
        sampsize = sampsize,
        importance = TRUE, proximity = TRUE, oob.prox = TRUE, keep.forest = TRUE
      )
      
      inner_probs <- predict(inner_rf, inner_test[, feats, drop = FALSE], type = "prob")
      calibration_probs  <- rbind(calibration_probs, inner_probs)
      calibration_labels <- c(calibration_labels, as.character(inner_test[[label_col]]))
    }
    
    calibration_model <- glmnet::cv.glmnet(
      x = as.matrix(calibration_probs),
      y = as.factor(calibration_labels),
      family = "multinomial", alpha = 0
    )
    
    test_probs <- predict(rf_model, test_data[, feats, drop = FALSE], type = "prob")
    calibrated_scores <- predict(calibration_model, newx = as.matrix(test_probs), type = "response")
    if (length(dim(calibrated_scores)) == 3) calibrated_scores <- calibrated_scores[, , 1, drop = TRUE]
    colnames(calibrated_scores) <- colnames(test_probs)
    calibrated_scores <- as.data.frame(calibrated_scores, check.names = FALSE)
    
    predicted_classes <- colnames(calibrated_scores)[max.col(calibrated_scores)]
    misclassification_errors <- c(misclassification_errors, mean(as.character(predicted_classes) != as.character(test_data[[label_col]])))
    
    auc_values <- c(auc_values, as.numeric(pROC::multiclass.roc(test_data[[label_col]], calibrated_scores)$auc))
    
    # keep the *last* fold's inner stacks, like your script ends up using
    last_calibration_probs  <- calibration_probs
    last_calibration_labels <- calibration_labels
  }
  
  list(
    mean_error = mean(misclassification_errors),
    mean_auc   = mean(auc_values),
    per_fold_error = misclassification_errors,
    per_fold_auc   = auc_values,
    last_calibration_probs  = last_calibration_probs,
    last_calibration_labels = last_calibration_labels
  )
}


##################################################################################################

Fit_RF_model <-  function(data, cv, label_col = "Subtype",
                          ntrees = 500, mtry = 6, seed = 123) {
  
feats <- setdiff(names(data), label_col)
class_distribution <- table(data$Subtype)
sampsize <- rep(min(class_distribution), length(class_distribution))

set.seed(seed)
rf_model <- randomForest::randomForest(
  x = data[, feats, drop = FALSE],
  y = data[[label_col]],
  ntree = ntrees, mtry = mtry,
  sampsize = sampsize,
  importance = TRUE, proximity = TRUE, oob.prox = TRUE, keep.forest = TRUE
)

rf_model

conf_matrix(
  df.true = data[[label_col]],
  df.pred = rf_model$predicted,
  title   = "Training Data Confusion Matrix - No calibration"
)

calibration_model <- glmnet::cv.glmnet(
  x = as.matrix(cv$last_calibration_probs),
  y = as.factor(cv$last_calibration_labels),
  family = "multinomial",
  alpha  = 0
)


test_probabilities <- predict(rf_model, data, type = "prob")
calibrated_probs <- predict(calibration_model, newx = as.matrix(test_probabilities), type = "response")
names <- colnames(calibrated_probs)
calibrated_probs_df <- as.data.frame(calibrated_probs)
colnames(calibrated_probs_df) <- names
datatable(calibrated_probs_df)

calibrated_probs_df$calls <-  colnames(calibrated_probs_df)[max.col(calibrated_probs_df)]
conf_matrix(
  df.true = data[[label_col]],
  df.pred = calibrated_probs_df$calls,
  title   = "Conf. Matrix on Training Data - Calibrated Calls"
)



list(
  rf_model          = rf_model,
  calibration_model = calibration_model,
  calibrated_table  = calibrated_probs_df,
  features          = feats,
  classes           = levels(data[[label_col]])
)
}






####################################################################################



Run_Pure_Model <- function( pure_rf_model = models$rf_model, pure_calibration_model = models$calibration_model, test_data = test, blanks_threshold = 21){
  
  test_data <- test_data[test_data$Blanks < blanks_threshold , ]
  m <- which(colnames(Query) == "Blanks")
  test_data <- test_data[, -m]


  test_probabilities <- predict(pure_rf_model, test_data, type = "prob")
  test_calls <- as.data.frame(predict(pure_rf_model, test_data))
  colnames(test_calls) <- "calls"
  test_df <- transform(merge(test_probabilities, test_calls, by = 0), row.names = 1)
  
  calibrated_probs <- predict(pure_calibration_model, newx = as.matrix(test_probabilities), type = "response")
  names <- colnames(calibrated_probs)
  calibrated_probs_df <- as.data.frame(calibrated_probs)
  colnames(calibrated_probs_df) <- names
  datatable(calibrated_probs_df)
  
  calibrated_probs_df$calls <-  colnames(calibrated_probs_df)[max.col(calibrated_probs_df)]
  
  Query <- Query[rownames(test_df), ]
  
  Pure_Probs <- list(pure_model = as.data.frame(test_df), calibrated_model = calibrated_probs_df, Query = Query)
  
  return(Pure_Probs)
} 



###################################################################################

