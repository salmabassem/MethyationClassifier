
library(dplyr)
library(randomForest)
library(Metrics)
library(caret)
library(DT)
library(ComplexHeatmap)
library(randomForest)
library(glmnet)
library(caret)
library(pROC)

source("R/confusion_matrix.R")
source("R/All_functions.R")


set.seed(123)


### Add Data 
MethData <- read.csv("Data/TrainingData.csv")
Epitypes <- as.data.frame(read.csv("Data/Epitypes.csv"))
Query <- read.csv("Data/TestData.csv", row.names = 1)


## Prepare Training Data 
iPlexData_Epi <-  Prepare_Training(MethData, Epitypes, Query)


## Train the models
cv <- nested_cv_calibrated_rf(iPlexData_Epi, ntrees = 500, mtry = 6, outer_folds = 3, inner_folds = 3, seed = 123)
print(cv$mean_error)
print(cv$mean_auc)


## Fit the Models 
models <- Fit_RF_model(iPlexData_Epi, cv = cv)


## Look at Confusion Matrix 
conf_matrix(df.true = iPlexData_Epi$Subtype, df.pred = models$rf_model$predicted,
           title = "Training Data Confusion Matrix - No calibration")

conf_matrix(df.true = iPlexData_Epi$Subtype, df.pred = models$calibrated_table$calls,
            title = "Training Data Confusion Matrix - No calibration")

calibrated_probs_df <- models$calibrated_table
DT::datatable(calibrated_probs_df)


# Run_Pure_Model on Query
preds <- Run_Pure_Model( test_data = Query)

conf_matrix(df.true = preds$Query$Subtype, df.pred = preds$pure_model$calls,
            title = "Conf. Matrix on Query - Calibrated Cells")

conf_matrix(df.true = preds$Query$Subtype, df.pred = preds$calibrated_model$calls,
            title = "Conf. Matrix on Query - Calibrated Cells")

