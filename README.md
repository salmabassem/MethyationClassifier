# AML iPLEX Subtyping Classifier (Random Forest + Calibration)

This repository contains a clean, reproducible R pipeline for AML iPLEX CpG subtyping using pure and mixed (impure) models. It:

- aligns the *training* CpGs to the *current query*
- trains a **Random Forest** (“Pure model”),
- fits a **multinomial ridge (glmnet) calibrator** to improve probability estimates,
- evaluates with **nested cross-validation**, and
- predicts on new data with **both** uncalibrated and calibrated outputs.

---

## Repository structure

```
.
├─ R/
│  ├─ All_functions.R          # Prepare_Training, nested_cv_calibrated_rf, Fit_RF_model, Run_Pure_Model
│  └─ confusion_matrix.R       # helper used by examples (conf_matrix)
├─ Data/
│  ├─ TrainingData.csv
│  ├─ Epitypes.csv
│  └─ TestData.csv
├─ outputs/                    # created by example script
└─ example_demo.R              # end-to-end example script
```

> Place your CSVs in `data/`. The example script reads from there.

---

## Requirements

- **R** ≥ 4.1
- Packages: `randomForest`, `glmnet`, `caret`, `pROC`, `DT` (optional, for viewing tables)

Install packages:
```r
install.packages(c("randomForest", "glmnet", "caret", "pROC", "DT"))
```

---

## Data expectations

- **Training matrix** — `data/TrainingData.csv`  
  Rows = samples; **first column = sample ID**. Remaining columns = CpGs (numeric).  
  Non-CpG columns like `sample.ID`, `Epitype.Color2` are automatically dropped.

- **Labels** — `data/Epitypes.csv`  
  Must contain `sample.ID` and `Subtype` (factor target).

- **Query** — `data/TestData.csv`  
  Rows = samples; columns = CpGs present in your assay. Optional columns:
  - `Subtype` (if you have ground truth for evaluation)
  - `Blanks` (QC column; rows with high `Blanks` can be filtered before prediction)

---

## Functions (in `R/All_functions.R`)

### `Prepare_Training(MethData, Epitypes, Query)`
Builds the training frame using **only CpGs shared** between training and query and attaches `Subtype`.  
**Returns:** `data.frame` with aligned CpGs + `Subtype` (factor).

### `nested_cv_calibrated_rf(data, ntrees=500, mtry=6, outer_folds=3, inner_folds=3, seed=123, label_col="Subtype")`
Outer K-fold CV; within each outer fold, fits an inner-fold **glmnet** calibrator and evaluates the **calibrated** outer test.  
**Returns:** `mean_error`, `mean_auc`, `per_fold_error`, `per_fold_auc`, and the **last** outer-fold’s `calibration_probs/labels` (to mirror the original Rmd workflow).

### `Fit_RF_model(data, cv, ntrees=500, mtry=6, seed=123, label_col="Subtype", show_table=FALSE)`
Trains the **final RF** on all rows, shows **OOB (uncalibrated)** confusion, fits the **final** glmnet calibrator using `cv$last_calibration_*`, and computes **calibrated** training predictions & confusion.  
**Returns:** `rf_model`, `calibration_model`, `calibrated_table`, `features`, `classes`.

### `Run_Pure_Model(test_data, pure_rf_model, pure_calibration_model, blanks_threshold=21, label_col="Subtype")`
Filters by `Blanks` (if present), aligns columns to the RF feature set, and returns **both** uncalibrated and calibrated predictions.  
**Returns:** list with `pure_model` (uncalibrated probs + `calls`) and `calibrated_model` (calibrated probs + `calls`).

---

## Quickstart

Run the full demo (train → evaluate → predict):

```bash
Rscript example_demo.R
```

Outputs:
- `outputs/predictions_uncalibrated.csv`
- `outputs/predictions_calibrated.csv`

---

## Example (inline)

```r
library(randomForest); library(glmnet); library(caret); library(pROC); library(DT)
source("R/confusion_matrix.R")
source("R/All_functions.R")

set.seed(123)

# 1) Load data (preserve CpG names)
MethData <- read.csv("data/TrainingData.csv", check.names = FALSE)
Epitypes <- read.csv("data/Epitypes.csv", check.names = FALSE)
Query    <- read.csv("data/TestData.csv",
                     row.names = 1, check.names = FALSE)

# 2) Prepare training aligned to Query CpGs
iPlexData_Epi <- Prepare_Training(MethData, Epitypes, Query)

# 3) Nested CV metrics (per-fold calibrator)
cv <- nested_cv_calibrated_rf(iPlexData_Epi,
                              ntrees = 500, mtry = 6,
                              outer_folds = 3, inner_folds = 3,
                              seed = 123)
cat("Mean misclassification error:", cv$mean_error, "\n")
cat("Mean multiclass AUC:",        cv$mean_auc,   "\n")

# 4) Fit final RF + final calibrator
models <- Fit_RF_model(iPlexData_Epi, cv = cv, ntrees = 500, mtry = 6, seed = 123)

# 5) Confusion matrices on training
conf_matrix(df.true = iPlexData_Epi$Subtype,
            df.pred = models$rf_model$predicted,
            title   = "Training — OOB (Uncalibrated)")

conf_matrix(df.true = iPlexData_Epi$Subtype,
            df.pred = models$calibrated_table$calls,
            title   = "Training — Calibrated")

# 6) Predict on Query (uncalibrated & calibrated)
preds <- Run_Pure_Model(
  test_data              = Query,
  pure_rf_model          = models$rf_model,
  pure_calibration_model = models$calibration_model,
  blanks_threshold       = 21
)

dir.create("outputs", showWarnings = FALSE)
write.csv(preds$pure_model,       "outputs/predictions_uncalibrated.csv", row.names = TRUE)
write.csv(preds$calibrated_model, "outputs/predictions_calibrated.csv",   row.names = TRUE)
```

---

## Tips

- Always read CSVs with `check.names = FALSE` to preserve CpG names.
- If you see **“No overlapping CpGs between training and query”**, ensure headers match and that you called `Prepare_Training()` with the same Query you will predict on.

---

