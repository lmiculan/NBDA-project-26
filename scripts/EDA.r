# Exploratoryy data analysis
# NBDA project - Leonardo Miculan, mat. 256108

library(recount3)
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(pheatmap)

# Create output directory for plots
if (!dir.exists("results/plots")) {
  dir.create("results/plots", recursive = TRUE)
}

# Download and prepare the data
projects <- available_projects()
proj_info <- subset(projects, project == "PRAD" & project_type == "data_sources" & file_source == "tcga")
rse <- create_rse(proj_info, type = "gene")

# Extract the count matrix (Rows = Genes, Columns = Samples)
countsDF <- as.data.frame(rse@assays@data$raw_counts)

# Extract metadata (Sample information)
colDataDF <- as.data.frame(rse@colData)
colDataDF <- colDataDF[colDataDF$tcga.cgc_sample_sample_type != "Metastatic", ]

# ------------------------------------
# 1. Full dataset health checks and plots
# ------------------------------------
raw_dims <- dim(countsDF)
raw_na <- sum(is.na(countsDF))
raw_zero_frac <- mean(countsDF == 0)
raw_gene_zero95 <- sum(rowMeans(countsDF == 0) > 0.95)
raw_sample_zero95 <- sum(colMeans(countsDF == 0) > 0.95)
raw_sample_type_counts <- colDataDF %>% dplyr::count(tcga.cgc_sample_sample_type)

cat("Full dataset dimensions:", raw_dims[1], "genes x", raw_dims[2], "samples
")
cat("Metadata dimensions:", dim(colDataDF)[1], "samples x", dim(colDataDF)[2], "fields
")
cat("NA values in raw counts:", raw_na, "
")
cat("Fraction of zeros in raw counts:", raw_zero_frac, "
")
cat("Genes with >95% zeros:", raw_gene_zero95, "
")
cat("Samples with >95% zeros:", raw_sample_zero95, "
")
cat("Sample type counts:
")
print(raw_sample_type_counts)

raw_row_sums <- rowSums(countsDF)
raw_col_sums <- colSums(countsDF)
raw_summary <- list(
  gene_sum_summary = summary(raw_row_sums),
  sample_sum_summary = summary(raw_col_sums)
)
print(raw_summary)

pdf("results/plots/00_FullDataset_RawCounts_Boxplot.pdf", width = 14, height = 6)
boxplot(log2(countsDF + 1), outline = FALSE, las = 2,
        main = "Full dataset raw count distribution (log2(count + 1))",
        ylab = "log2(count + 1)")
dev.off()

# ------------------------------------
# 2. Subset and normalize the data
# ------------------------------------
set.seed(123)
target_n <- 150
colDataDF_down <- colDataDF %>%
  group_by(tcga.cgc_sample_sample_type) %>%
  group_modify(~ {
    class_prop <- nrow(.x) / nrow(colDataDF)
    class_n <- round(target_n * class_prop)
    slice_sample(.x, n = class_n)
  }) %>%
  ungroup()

stopifnot(all(colDataDF_down$external_id %in% colnames(countsDF)))
countsDF_down <- countsDF[, colDataDF_down$external_id]
filtered_countsDF <- countsDF_down[rowSums(countsDF_down >= 10) >= 10, ]

dds <- DESeqDataSetFromMatrix(countData = filtered_countsDF,
                              colData = colDataDF_down,
                              design = ~ tcga.cgc_sample_sample_type)
vst_counts <- vst(dds, blind = TRUE)

# Calculate log2FC
dds <- DESeq(dds)
res <- results(dds, contrast = c("tcga.cgc_sample_sample_type", "Primary Tumor", "Solid Tissue Normal"))

# Extract normalized matrix only
vst_matrix <- assay(vst_counts)

# Keep the top 10,000 most variable genes in the subset (safe to handle <10000 genes)
var_by_gene <- apply(vst_matrix, 1, var)
n_top <- min(10000, length(var_by_gene))
top_idx <- order(var_by_gene, decreasing = TRUE)[1:n_top]
filtered_vst_matrix <- vst_matrix[top_idx, , drop = FALSE]

# Subset the DESeq2 results to the same genes
top_genes <- rownames(filtered_vst_matrix)
res_filtered <- res[top_genes, , drop = FALSE]
log2FCdf <- as.data.frame(res_filtered) |>
  rownames_to_column("gene") |>
  rename(log2FC = log2FoldChange) |>
  select(gene, log2FC)

save(list = c("filtered_vst_matrix", "colDataDF_down", "pca_data"), file = "filtered_vst_matrix.rData")
save(list = c("log2FCdf"), file = "log2FCDF.rData")

pca_data <- prcomp(t(filtered_vst_matrix))

# ------------------------------------
# 3. Plots for the subsetted normalized dataset
# ------------------------------------
# Normalization verification boxplots
vst_long <- as.data.frame(filtered_vst_matrix) |>
  rownames_to_column("gene") |>
  pivot_longer(cols = -gene, names_to = "sample", values_to = "expression") |>
  left_join(
    colDataDF_down |>
      select(external_id, tcga.cgc_sample_sample_type),
    by = c("sample" = "external_id")
  )

p1 <- ggplot(vst_long, aes(x = tcga.cgc_sample_sample_type, y = expression, fill = tcga.cgc_sample_sample_type)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  labs(title = "Subset VST expression distribution by sample type",
       x = "Sample type",
       y = "Expression (VST)",
       fill = "Sample type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("results/plots/01_Subset_VST_Boxplot_BySampleType.pdf", plot = p1, width = 8, height = 6, dpi = 300)

p2 <- ggplot(vst_long, aes(x = reorder(sample, expression, FUN = median), y = expression)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7, outlier.shape = NA) +
  labs(title = "Subset VST expression distribution across samples",
       x = "Sample",
       y = "Expression (VST)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))

ggsave("results/plots/02_Subset_VST_Boxplot_AllSamples.pdf", plot = p2, width = 14, height = 6, dpi = 300)

# Scree plot from PCA on subset data
variance_explained <- summary(pca_data)$importance[2, ]
cumulative_variance <- summary(pca_data)$importance[3, ]

scree_df <- tibble(
  PC = factor(1:length(variance_explained), levels = seq_along(variance_explained)),
  variance = variance_explained,
  cumulative = cumulative_variance
)

p3 <- ggplot(scree_df[1:30, ], aes(x = PC)) +
  geom_point(aes(y = variance), size = 3, color = "steelblue") +
  geom_line(aes(y = variance, group = 1), color = "steelblue", linewidth = 1) +
  geom_line(aes(y = cumulative, group = 1), color = "coral", linewidth = 1, linetype = "dashed") +
  labs(title = "Subset scree plot: variance explained by PCs",
       x = "Principal component",
       y = "Proportion of variance explained",
       subtitle = "Solid = individual variance, dashed = cumulative variance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("results/plots/03_Subset_Scree_Plot.pdf", plot = p3, width = 10, height = 6, dpi = 300)

# PCA scatterplot on subset samples
pca_df <- as.data.frame(pca_data$x) %>%
  rownames_to_column("sample") %>%
  left_join(
    colDataDF_down |>
      select(external_id, tcga.cgc_sample_sample_type),
    by = c("sample" = "external_id")
  )

p4 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = tcga.cgc_sample_sample_type)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "PCA of VST-normalized subset samples (PC1 vs PC2)",
       x = paste0("PC1 (", round(variance_explained[1] * 100, 2), "%)"),
       y = paste0("PC2 (", round(variance_explained[2] * 100, 2), "%)"),
       color = "Sample type") +
  theme_minimal()

ggsave("results/plots/04_Subset_PCA_PC1_PC2.pdf", plot = p4, width = 8, height = 6, dpi = 300)

# Histograms on the subsetted VST distribution
hist_df <- vst_long %>% sample_n(min(nrow(vst_long), 50000))

p5 <- ggplot(hist_df, aes(x = expression)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, color = "black") +
  labs(title = "Histogram of subset VST expression values",
       x = "Expression (VST)",
       y = "Frequency") +
  theme_minimal()

ggsave("results/plots/05_Subset_VST_Histogram.pdf", plot = p5, width = 8, height = 6, dpi = 300)

p6 <- ggplot(hist_df, aes(x = expression, fill = tcga.cgc_sample_sample_type)) +
  geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
  facet_wrap(~ tcga.cgc_sample_sample_type) +
  labs(title = "Histogram of subset VST expression by sample type",
       x = "Expression (VST)",
       y = "Frequency",
       fill = "Sample type") +
  theme_minimal()

ggsave("results/plots/06_Subset_VST_Histogram_ByType.pdf", plot = p6, width = 10, height = 6, dpi = 300)

# Scatterplot of two top variable genes in the subset
var_genes <- sort(var_by_gene, decreasing = TRUE)
top_genes <- names(var_genes)[1:2]
scatter_df <- tibble(
  gene1 = filtered_vst_matrix[top_genes[1], ],
  gene2 = filtered_vst_matrix[top_genes[2], ],
  sample_type = colDataDF_down$tcga.cgc_sample_sample_type
)

p7 <- ggplot(scatter_df, aes(x = gene1, y = gene2, color = sample_type)) +
  geom_point(size = 3, alpha = 0.75) +
  geom_smooth(method = "lm", se = TRUE, color = "black", alpha = 0.2) +
  labs(title = paste("Scatterplot of top two variable genes:", top_genes[1], "vs", top_genes[2]),
       x = paste0("Expression: ", top_genes[1]),
       y = paste0("Expression: ", top_genes[2]),
       color = "Sample type") +
  theme_minimal()

ggsave("results/plots/07_Subset_GeneScatter_Top2.pdf", plot = p7, width = 8, height = 6, dpi = 300)

# Heatmaps for the subsetted data
heatmap_n <- 50
heatmap_genes <- order(apply(filtered_vst_matrix, 1, var), decreasing = TRUE)[1:heatmap_n]
heatmap_data <- filtered_vst_matrix[heatmap_genes, ]
annotation_df <- data.frame(
  SampleType = colDataDF_down$tcga.cgc_sample_sample_type,
  row.names = colnames(heatmap_data)
)
annotation_colors <- list(
  SampleType = c(
    "Primary Tumor" = "#F8766D",
    "Solid Tissue Normal" = "#00BFC4"
  )
)

pheatmap(heatmap_data,
         annotation_col = annotation_df,
         annotation_colors = annotation_colors,
         scale = "row",
         clustering_distance_cols = "euclidean",
         clustering_distance_rows = "euclidean",
         clustering_method = "complete",
         main = "Heatmap of top 50 variable genes in subset",
         fontsize = 8,
         show_rownames = FALSE,
         show_colnames = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         filename = "results/plots/08_Subset_Heatmap_Top50VariableGenes.pdf",
         width = 12,
         height = 10)

p_values <- apply(filtered_vst_matrix, 1, function(gene_values) {
  groups <- split(gene_values, colDataDF_down$tcga.cgc_sample_sample_type)
  t.test(groups[[1]], groups[[2]])$p.value
})

ttest_genes <- order(p_values)[1:heatmap_n]
heatmap_data_ttest <- filtered_vst_matrix[ttest_genes, ]


pheatmap(heatmap_data_ttest,
         annotation_col = annotation_df,
         annotation_colors = annotation_colors,
         scale = "row",
         clustering_distance_cols = "euclidean",
         clustering_distance_rows = "euclidean",
         clustering_method = "complete",
         main = "Heatmap of top 50 significant genes in subset",
         fontsize = 8,
         show_rownames = FALSE,
         show_colnames = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         filename = "results/plots/09_Subset_Heatmap_Top50SignificantGenes.pdf",
         width = 12,
         height = 10)


# Additional QC summary for the subset
summary_stats <- vst_long %>%
  group_by(tcga.cgc_sample_sample_type) %>%
  summarise(
    mean = mean(expression, na.rm = TRUE),
    median = median(expression, na.rm = TRUE),
    sd = sd(expression, na.rm = TRUE),
    q25 = quantile(expression, 0.25, na.rm = TRUE),
    q75 = quantile(expression, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

p8 <- ggplot(vst_long, aes(x = tcga.cgc_sample_sample_type, y = expression, fill = tcga.cgc_sample_sample_type)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.2, alpha = 0.8) +
  labs(title = "Subset VST expression by sample type (violin + boxplot)",
       x = "Sample type",
       y = "Expression (VST)",
       fill = "Sample type") +
  theme_minimal()

ggsave("results/plots/10_Subset_Violin_Boxplot.pdf", plot = p8, width = 8, height = 6, dpi = 300)

summary_stats
