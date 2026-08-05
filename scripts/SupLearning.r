# Supervised learning
# NBDA project - Leonardo Miculan, mat. 256108

library(tidyverse)
library(edgeR)
library(caret)
library(randomForest)
library(glmnet)
library(pROC)
library(DESeq2)
library(ggplot2)
library(ggvenn)
library(patchwork)

rm(list = ls())
gc()

load("filtered_vst_matrix.rData")

# 1. Prepare data (samples as rows, genes as columns)
model_data <- as.data.frame(t(filtered_vst_matrix))
model_data$Class <- as.factor(make.names(colDataDF_down$tcga.cgc_sample_sample_type))

pos_class <- levels(model_data$Class)[1]
neg_class <- levels(model_data$Class)[2]

# 2. Outer CV folds (stratified, preserves class ratio in each fold)
set.seed(42)
k_outer <- 10
outer_folds <- createFolds(model_data$Class, k = k_outer, list = TRUE, returnTrain = FALSE)

# Inner tuning control — identical to before, just now nested inside each outer fold
inner_control <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

# Containers for pooled out-of-fold predictions
pred_rf     <- list()
pred_svm    <- list()
pred_glmnet <- list()

for (i in seq_along(outer_folds)) {
  cat(sprintf("\n===== Outer fold %d / %d =====\n", i, k_outer))
  
  test_idx  <- outer_folds[[i]]
  train_idx <- setdiff(seq_len(nrow(model_data)), test_idx)
  
  outer_train <- model_data[train_idx, ]
  outer_test  <- model_data[test_idx, ]
  
  # ---- Random Forest ----
  set.seed(123)
  cat("  Training Random Forest...\n")
  p_features <- length(rownames(filtered_vst_matrix))
  tune_grid_rf <- expand.grid(
    mtry = c(
      floor(0.5 * sqrt(p_features)),
      floor(sqrt(p_features)),
      floor(2 * sqrt(p_features)),
      floor(0.05 * p_features),
      floor(0.10 * p_features)))
  
  fold_rf <- caret::train(Class ~ .,
                        data = outer_train,
                        tuneGrid = tune_grid_rf,
                        method = "rf",
                        trControl = inner_control,
                        metric = "ROC",
                        ntree = 1000)
  
  # ---- SVM (Radial) — sigma re-estimated on THIS fold's training data only ----
  set.seed(123)
  cat("  Training SVM...\n")
  sigma_est <- kernlab::sigest(as.matrix(outer_train %>% dplyr::select(-Class)))
  sigma_med <- as.numeric(sigma_est[2])
  C_vals    <- 2^(-3:3)
  tuneGrid  <- expand.grid(sigma = sigma_med, C = C_vals)
  
  fold_svm <- caret::train(Class ~ ., data = outer_train, method = "svmRadial",
                            trControl = inner_control, metric = "ROC",
                            tuneGrid = tuneGrid,
                            preProc = c("center", "scale"))
  
  # ---- GLMNET ----
  set.seed(123)
  cat("  Training GLMNET...\n")
  fold_glmnet <- caret::train(Class ~ ., data = outer_train, method = "glmnet",
                               trControl = inner_control, metric = "ROC",
                               preProc = c("center", "scale"))
  
  # ---- Predict on the untouched outer test fold ----
  prob_rf     <- predict(fold_rf,     newdata = outer_test, type = "prob")[[pos_class]]
  prob_svm    <- predict(fold_svm,    newdata = outer_test, type = "prob")[[pos_class]]
  prob_glmnet <- predict(fold_glmnet, newdata = outer_test, type = "prob")[[pos_class]]
  
  pred_rf[[i]]     <- data.frame(obs = outer_test$Class, prob = prob_rf,     fold = i)
  pred_svm[[i]]    <- data.frame(obs = outer_test$Class, prob = prob_svm,    fold = i)
  pred_glmnet[[i]] <- data.frame(obs = outer_test$Class, prob = prob_glmnet, fold = i)
}

# 3. Pool predictions across all outer folds (this is now effectively your full ~100 samples)
pooled_rf     <- bind_rows(pred_rf)
pooled_svm    <- bind_rows(pred_svm)
pooled_glmnet <- bind_rows(pred_glmnet)

# 4. Build ROC objects from pooled, never-trained-on predictions
roc_rf     <- roc(pooled_rf$obs,     pooled_rf$prob,     levels = c(neg_class, pos_class))
roc_svm    <- roc(pooled_svm$obs,    pooled_svm$prob,    levels = c(neg_class, pos_class))
roc_glmnet <- roc(pooled_glmnet$obs, pooled_glmnet$prob, levels = c(neg_class, pos_class))

auc_rf     <- round(auc(roc_rf), 3)
auc_svm    <- round(auc(roc_svm), 3)
auc_glmnet <- round(auc(roc_glmnet), 3)

# Confidence intervals — worth reporting given the earlier discussion of estimate stability
ci_rf     <- ci.auc(roc_rf)
ci_svm    <- ci.auc(roc_svm)
ci_glmnet <- ci.auc(roc_glmnet)

cat("\n--- Nested CV AUCs (pooled across", k_outer, "outer folds) ---\n")
cat(sprintf("RF:     %.3f  [%.3f, %.3f]\n", auc_rf, ci_rf[1], ci_rf[3]))
cat(sprintf("SVM:    %.3f  [%.3f, %.3f]\n", auc_svm, ci_svm[1], ci_svm[3]))
cat(sprintf("GLMNET: %.3f  [%.3f, %.3f]\n", auc_glmnet, ci_glmnet[1], ci_glmnet[3]))

# 5. Publication-ready combined plot — same styling as before, now from pooled nested-CV predictions
roc_plot <- ggroc(
  list(
    "Random Forest" = roc_rf,
    "SVM (Radial)"  = roc_svm,
    "GLMNET"        = roc_glmnet
  ),
  legacy.axes = TRUE,
  size = 1.1
) +
  geom_segment(
    aes(x = 0, xend = 1, y = 0, yend = 1),
    color = "grey50",
    linetype = "dashed"
  ) +
  scale_color_manual(
    values = c(
      "Random Forest" = "#E64B35",
      "SVM (Radial)"  = "#4DBBD5",
      "GLMNET"        = "#00A087"
    ),
    labels = c(
      paste0("Random Forest (AUC = ", auc_rf, ")"),
      paste0("SVM (Radial) (AUC = ", auc_svm, ")"),
      paste0("GLMNET (AUC = ", auc_glmnet, ")")
    )
  ) +
  labs(
    title = "ROC Curves — Nested Cross-Validation",
    subtitle = paste0("Target Phenotype: ", pos_class, " | Pooled predictions across ", k_outer, " outer folds"),
    x = "1 - Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)",
    color = "Model Architecture"
  ) +
  theme_minimal() +
  theme(
    legend.position = c(0.70, 0.25),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
    plot.title = element_text(face = "bold", size = 13),
    axis.title = element_text(face = "bold")
  )

print(roc_plot)
ggsave("results/ROC_Curves_NestedCV.pdf", plot = roc_plot, width = 7, height = 6, dpi = 300)

# 6. Final prediction models - fit on the FULL dataset for downstream variable importance / DEG extraction
# (nested CV above answers "how good is this approach"; these give you the actual models to extract genes from)
set.seed(123)
model_rf_final <- caret::train(Class ~ ., data = model_data, method = "rf",
                                trControl = inner_control, tuneGrid = tune_grid_rf, metric = "ROC")

set.seed(123)
sigma_est_final <- kernlab::sigest(as.matrix(model_data %>% dplyr::select(-Class)))
sigma_med_final <- as.numeric(sigma_est_final[2])
model_svm_final <- caret::train(Class ~ ., data = model_data, method = "svmRadial",
                                 trControl = inner_control, metric = "ROC",
                                 tuneGrid = expand.grid(sigma = sigma_med_final, C = 2^(-3:3)),
                                 preProc = c("center", "scale"))

set.seed(123)
model_glmnet_final <- caret::train(Class ~ ., data = model_data, method = "glmnet",
                                    trControl = inner_control, metric = "ROC",
                                    preProc = c("center", "scale"))

# Save final model values
capture.output(print(model_rf_final), file = "results/model_rf_final.txt")
capture.output(print(model_svm_final), file = "results/model_svm_final.txt")
capture.output(print(model_glmnet_final), file = "results/model_glmnet_final.txt")

save(list = c("model_rf_final","model_svm_final","model_glmnet_final"), file = "final_models.rData")

varImp_rf     <- varImp(model_rf_final)$importance
varImp_svm    <- varImp(model_svm_final)$importance
varImp_glmnet <- varImp(model_glmnet_final)$importance

# 1. Robust function to extract top N genes
# Handles models with an "Overall" column (RF, GLMNET) 
# and models with class-specific importance columns (SVM)
get_top_genes <- function(varImp_df, top_n = 500) {
  if ("Overall" %in% colnames(varImp_df)) {
    importance_vec <- varImp_df$Overall
  } else {
    # Compute rowMeans across class-specific columns to derive a unified importance score
    importance_vec <- varImp_df %>% dplyr::select(where(is.numeric)) %>% rowMeans(na.rm = TRUE)
  }

  df <- data.frame(gene = rownames(varImp_df), importance = importance_vec) %>%
    arrange(desc(importance)) %>%
    slice_head(n = top_n)
  
  return(df)
}

# Extract the top 500 features for each model
top_rf     <- get_top_genes(varImp_rf, top_n = 750)
top_svm    <- get_top_genes(varImp_svm, top_n = 750)
top_glmnet <- get_top_genes(varImp_glmnet, top_n = 750)

# 2. Assign directionality (Up/Down) using log2FC without hard filtering
classify_direction <- function(top_genes_df, log2fc_df) {
  merged <- top_genes_df %>%
    left_join(log2fc_df, by = "gene") %>%
    mutate(Regulation = ifelse(log2FC > 0, "Up", "Down"))
  
  return(merged)
}

# Import log2 fold changes for all genes in the filtered VST matrix
load("log2FCDF.rData")

# Apply classification to all three models
DEGs_rf     <- classify_direction(top_rf, log2FCdf)
DEGs_svm    <- classify_direction(top_svm, log2FCdf)
DEGs_glmnet <- classify_direction(top_glmnet, log2FCdf)

# 3. Extract gene vectors for downstream consensus intersection
up_rf     <- DEGs_rf$gene[DEGs_rf$Regulation == "Up"]
up_svm    <- DEGs_svm$gene[DEGs_svm$Regulation == "Up"]
up_glmnet <- DEGs_glmnet$gene[DEGs_glmnet$Regulation == "Up"]

down_rf     <- DEGs_rf$gene[DEGs_rf$Regulation == "Down"]
down_svm    <- DEGs_svm$gene[DEGs_svm$Regulation == "Down"]
down_glmnet <- DEGs_glmnet$gene[DEGs_glmnet$Regulation == "Down"]

# 7. Venn diagrams for Up- and Down-regulated intersections
venn_up <- ggvenn(
  list(RF = up_rf, SVM = up_svm, GLMNET = up_glmnet),
  fill_color = c("#E64B35", "#4DBBD5", "#00A087"),
  stroke_size = 0.5, set_name_size = 4
) + labs(title = "Up-regulated DEG overlap across models")

venn_down <- ggvenn(
  list(RF = down_rf, SVM = down_svm, GLMNET = down_glmnet),
  fill_color = c("#E64B35", "#4DBBD5", "#00A087"),
  stroke_size = 0.5, set_name_size = 4
) + labs(title = "Down-regulated DEG overlap across models")

venn_up + venn_down

ggsave("results/DEG_Venn_Diagrams.pdf", width = 12, height = 6, dpi = 300)

# 8. Final consensus DEG lists: genes flagged consistently by all three models
consensus_up_majority <- Reduce(union, list(
  intersect(up_rf, up_svm),
  intersect(up_rf, up_glmnet),
  intersect(up_svm, up_glmnet)
))
consensus_down_majority <- Reduce(union, list(
  intersect(down_rf, down_svm),
  intersect(down_rf, down_glmnet),
  intersect(down_svm, down_glmnet)
))

consensus_DEGs <- bind_rows(
  data.frame(gene = consensus_up_majority, Regulation = "Up"),
  data.frame(gene = consensus_down_majority, Regulation = "Down")
)

consensus_DEGs
save(list = c("consensus_DEGs", "log2fc_df"), file = "consensus_DEGs.rData")
