# =============================================================================
# #18 -- Drug-Resistance Prediction from Genotype (Staphylococcus aureus) -- R twin
# Reads the REAL BV-BRC-derived bridge tables 01_build.ipynb exports into data_py/:
#   meth_frame.csv  -- methicillin modeling frame (genome_id, <AMR gene cols>, label, mlst)
#   ery_frame.csv   -- erythromycin modeling frame, same shape
# Run the Python notebook through Step 2b FIRST so these exist.
# Study PROJECT_NARRATIVE.md and WORKFLOW.md first -- this is the SECOND real
# dataset pivot in this project (K. pneumoniae -> S. aureus, see Section 7b) and
# Step 2's feature matrix is a broad, unfiltered set of curated AMR genes, not a
# hand-picked mecA/erm list -- see the Python notebook's Step 2 markdown for why.
# ~80-90% of this file is blanks (____) for you to fill in -- each blank has the
# exact real code as a trailing comment hint, same as 01_build.ipynb.
# =============================================================================

## ---- Step 0 -- load packages + real data (GIVEN) --------------------------
install.packages(c("randomForest", "pROC", "dplyr"))
library(dplyr)
library(randomForest)
library(pROC)

setwd("D:/GIT HUB PROJECTS/amr-genotype-prediction")
dir.create("results_R", showWarnings = FALSE)

meth_frame <- read.csv("data_py/meth_frame.csv", stringsAsFactors = FALSE)
ery_frame  <- read.csv("data_py/ery_frame.csv", stringsAsFactors = FALSE)

cat(nrow(meth_frame), "methicillin isolates,", nrow(ery_frame), "erythromycin isolates\n")
cat("methicillin columns:", length(colnames(meth_frame)), "total (mostly AMR gene features)\n")
cat("erythromycin columns:", length(colnames(ery_frame)), "total (mostly AMR gene features)\n")

## ---- Step 1 -- identify feature columns  <- WRITE IT -----------------------
## WHAT: the feature columns are the AMR-gene presence/absence columns -- every
## curated AMR gene seen in enough genomes (Step 2 of the Python notebook), NOT
## a hand-picked mecA/erm list -- everything in the frame except genome_id,
## label, and mlst.
## WHY: evaluate_model() in Step 3 needs to know exactly which columns are real
## predictors vs. bookkeeping columns (id/label/split-key).
## Docs: base R setdiff() to exclude known non-feature column names from colnames().
get_feature_cols <- function(frame) {
  setdiff(colnames(frame),c("genome_id","label","mlst"))                                        # setdiff(colnames(frame), c("genome_id", "label", "mlst"))
}

meth_feature_cols <- get_feature_cols(meth_frame)                     # get_feature_cols(meth_frame)
ery_feature_cols  <- get_feature_cols(ery_frame)                     # get_feature_cols(ery_frame)
cat("methicillin features:", length(meth_feature_cols), "AMR genes\n")
cat("erythromycin features:", length(ery_feature_cols), "AMR genes\n")

get_feature_cols <- function(frame) {
  setdiff(colnames(frame), c("genome_id", "label", "mlst"))
}

meth_feature_cols <- get_feature_cols(meth_frame)
ery_feature_cols  <- get_feature_cols(ery_frame)
cat("methicillin features:", length(meth_feature_cols), "AMR genes\n")
cat("erythromycin features:", length(ery_feature_cols), "AMR genes\n")
## ---- Step 2 -- MLST-aware split vs random split  <- WRITE IT ---------------
## WHAT: build TWO splits per drug -- (a) MLST-aware: hold out WHOLE sequence
## types for test; (b) plain random: for direct comparison in Step 3.
## WHY: isolates within one MLST sequence type share genetic background for
## reasons unrelated to any single resistance gene -- the same structure-
## leakage risk as TB lineage (project #1) or QSAR's scaffold split (#15).
## globally disseminated MRSA clones (e.g. ST5, ST8/USA300, ST22/EMRSA-15) mean
## a random split can let the model partly learn "which clone is this" as a
## shortcut instead of the actual resistance mechanism.
## Docs: base R `sample()` on unique(frame$mlst) to pick held-out sequence
## types; boolean indexing with %in% to split rows.
mlst_aware_split <- function(frame, test_frac = 0.2, seed = 0) {
  set.seed(seed)
  all_mlsts <- unique(frame$mlst)                           # unique(frame$mlst)
  test_mlsts <- sample(all_mlsts, size = round(test_frac * length(all_mlsts)))                          # sample(all_mlsts, size = round(test_frac * length(all_mlsts)))
  is_test <- frame$mlst %in% test_mlsts                             # frame$mlst %in% test_mlsts
  list(train = frame[!is_test, ], test = frame[is_test, ])
}

random_split <- function(frame, test_frac = 0.2, seed = 0) {
  set.seed(seed)
  n <- nrow(frame)
  test_idx <- sample(seq_len(n), size = round(test_frac * n))                            # sample(seq_len(n), size = round(test_frac * n))
  list(train = frame[-test_idx, ], test = frame[test_idx, ])
}

meth_mlst_split   <- mlst_aware_split(meth_frame)                     # mlst_aware_split(meth_frame)
meth_random_split <- random_split(meth_frame)                    # random_split(meth_frame)
ery_mlst_split    <- mlst_aware_split(ery_frame)                     # mlst_aware_split(ery_frame)
ery_random_split  <- random_split(ery_frame)                     # random_split(ery_frame)

## ---- Step 3 -- train + evaluate per drug, per split, per model  <- WRITE IT
## WHAT: fit BOTH a logistic regression (glm) and a random forest, per drug,
## per split. Report sensitivity, specificity, AUROC -- same metrics as Python.
## WHY: a false negative (resistant isolate called susceptible) risks
## treatment failure -- sensitivity and specificity matter separately, not
## just overall accuracy.
## Docs: stats::glm(label ~ ., data = train[, c(feature_cols, "label")], family = "binomial");
## randomForest::randomForest(as.factor(label) ~ ., data = train[, c(feature_cols, "label")]);
## pROC::roc(response, predictor) for AUROC; table(predicted, actual) for the
## confusion matrix underlying sensitivity/specificity.
evaluate_model <- function(train, test, feature_cols, model_type = c("logistic", "rf")) {
  model_type <- match.arg(model_type)
  train_sub <- train[, c(feature_cols, "label")]
  test_sub  <- test[, c(feature_cols, "label")]

  if (model_type == "logistic") {
    model <- glm(label ~ ., data = train_sub, family = "binomial")                             # glm(label ~ ., data = train_sub, family = "binomial")
    pred_prob <- predict(model, newdata = test_sub, type = "response")                         # predict(model, newdata = test_sub, type = "response")
  } else {
    train_sub$label <- as.factor(train_sub$label)
    model <- randomForest(label ~ ., data = train_sub)                             # randomForest(label ~ ., data = train_sub)
    pred_prob <- predict(model, newdata = test_sub, type = "prob")[, "1"]                         # predict(model, newdata = test_sub, type = "prob")[, "1"]
  }

  pred_class <- ifelse(pred_prob >= 0.5, 1, 0)                          # ifelse(pred_prob >= 0.5, 1, 0)
  cm <- cm <- table(predicted = factor(pred_class, levels = c(0, 1)),
                    actual = factor(test_sub$label, levels = c(0, 1)))                                  # table(predicted = pred_class, actual = test_sub$label)
  tp <- cm["1", "1"]                                  # cm["1", "1"]  -- guard: 0 if that cell is missing
  tn <- cm["0", "0"]                                  # cm["0", "0"]
  fp <- cm["1", "0"]                                  # cm["1", "0"]
  fn <- cm["0", "1"]                                  # cm["0", "1"]
  sensitivity <- tp / (tp + fn)                         # tp / (tp + fn)
  specificity <- tn / (tn + fp)                         # tn / (tn + fp)
  auroc <- as.numeric(pROC::roc(test_sub$label, pred_prob, quiet = TRUE)$auc)                               # as.numeric(pROC::roc(test_sub$label, pred_prob, quiet = TRUE)$auc)

  data.frame(sensitivity = sensitivity, specificity = specificity, auroc = auroc)
}

drug_splits <- list(
  methicillin  = list(`mlst-aware` = meth_mlst_split, random = meth_random_split),
  erythromycin = list(`mlst-aware` = ery_mlst_split, random = ery_random_split)
)
feature_cols_by_drug <- list(methicillin = meth_feature_cols, erythromycin = ery_feature_cols)
model_types <- c("logistic", "rf")

results_list <- list()
 for (drug in names(drug_splits)) {
   for (split_name in names(drug_splits[[drug]])) {
    split <- drug_splits[[drug]][[split_name]]
     for (model_type in model_types) {
       metrics <- evaluate_model(split$train, split$test, feature_cols_by_drug[[drug]], model_type)
       metrics$drug <- drug; metrics$split <- split_name; metrics$model <- model_type
       results_list[[length(results_list) + 1]] <- metrics
     }
   }
 }

results_df <- do.call(rbind, results_list)                            # do.call(rbind, results_list)
print(results_df)
write.csv(results_df, "results_R/model_performance_R.csv", row.names = FALSE)

## ---- Step 4 -- feature importance: does the model rediscover mecA/erm?  <- WRITE IT (HERO)
## WHAT: for each drug, extract + print which AMR-gene features the
## mlst-aware-split models weighted most heavily.
## WHY: this is the real payoff of the broad, unfiltered feature matrix -- the
## model was never told which gene family matters. If methicillin's top
## features are mecA/mecA_1/mecA_2/mecI/mecR1, and erythromycin's top features
## are from the erm family (Erm(A)/Erm(B)/Erm(C)/ErmA/ErmC/ErmT/...), that's a
## real, independent confirmation of published AMR biology -- not something
## assumed into the features ahead of time the way the original K. pneumoniae
## design did.
## Docs: coef(model) for glm (drop the intercept); importance(model) for randomForest.
get_feature_importance <- function(fitted_model, model_type) {
  if (model_type == "logistic") {
    coefs <- coef(fitted_model)[-1]                              # coef(fitted_model)[-1]  -- drop (Intercept)
    sort(coefs, decreasing = TRUE)                                      # sort(coefs, decreasing = TRUE)
  } else {
    imp <- randomForest::importance(fitted_model)[, "MeanDecreaseGini"]                               # randomForest::importance(fitted_model)[, "MeanDecreaseGini"]
    sort(imp, decreasing = TRUE)                                      # sort(imp, decreasing = TRUE)
  }
}

feature_importance_records <- list()
 for (drug in names(drug_splits)) {
   split <- drug_splits[[drug]][["mlst-aware"]]
   feature_cols <- feature_cols_by_drug[[drug]]
   train_sub <- split$train[, c(feature_cols, "label")]
   for (model_type in model_types) {
     if (model_type == "logistic") {
       fitted <- glm(label ~ ., data = train_sub, family = "binomial")
     } else {
       train_sub_rf <- train_sub; train_sub_rf$label <- as.factor(train_sub_rf$label)
       fitted <- randomForest(label ~ ., data = train_sub_rf)
     }
     ranked <- get_feature_importance(fitted, model_type)
     cat("\n", drug, "/", model_type, "top 10 features:\n"); print(head(ranked, 10))
     feature_importance_records[[length(feature_importance_records) + 1]] <- data.frame(
       drug = drug, model = model_type, feature = names(ranked), importance = as.numeric(ranked)
     )
   }
 }
# compare the top features above against the EXPECTED_GENE_MARKERS dict in
# 01_build.ipynb's Step 2 -- did mecA/erm-family genes actually rise to the top?

feature_importance_df <- do.call(rbind, feature_importance_records)                 # do.call(rbind, feature_importance_records)
write.csv(feature_importance_df, "results_R/feature_importance_R.csv", row.names = FALSE)

# =============================================================================
# INTERPRET (write as comments, same 5 questions as WORKFLOW.md Section 7):
# 1. Computed?  2. Expected?  3. Match (your real sensitivity/specificity/AUROC
#    per drug per split)?  4. Confirms (did feature importance -- run on a
#    BROAD, unfiltered gene set -- actually rediscover mecA-family genes for
#    methicillin and erm-family genes for erythromycin, on its own)?
#    5. Caveat -- compare your R results to Python's; note any differences and
#    whether they trace to glm vs LogisticRegression regularization defaults,
#    randomForest vs scikit-learn RF hyperparameter defaults, or genuine
#    RNG/split differences.
# =============================================================================
