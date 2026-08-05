# WCGNA analysis for DEGs list
# NBDA project - Leonardo Miculan, mat. 256108

invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE))
library(WGCNA)
library(tidyverse)
library(ggplot2)
library(patchwork)

rm(list = ls())

# Enable multi-threading for WGCNA functions
load("filtered_vst_matrix.rData")
load("consensus_DEGs.rData")

allowWGCNAThreads()

# 1. Prepare expression data (transpose to rows = samples, columns = genes)
datExpr <- as.matrix(t(filtered_vst_matrix))
storage.mode(datExpr) <- "double"

# then run goodSamplesGenes and the rest as before
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0) print("Removing bad genes...")
  if (sum(!gsg$goodSamples) > 0) print("Removing bad samples...")
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# align sample labels to filtered samples
sample_labels <- colDataDF_down$tcga.cgc_sample_sample_type
sample_labels <- sample_labels[gsg$goodSamples]

# 2. Choose a soft-thresholding power (beta)
powers <- c(c(1:10), seq(from = 8, to = 20, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the scale-free topology fit index
df_fit <- data.frame(
  power = sft$fitIndices[, 1],
  fit = -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  meanConnectivity = sft$fitIndices[, 5]
)

softthresholding_plot <- ggplot(df_fit, aes(x = power, y = fit)) +
  geom_point(color = "red") +
  geom_point(data = df_fit[df_fit$power == sft$powerEstimate, ], aes(x = power, y = fit), pch=21, fill=NA, size=4, colour="black", stroke=1) +
  geom_text(data = df_fit[df_fit$power == sft$powerEstimate, ], aes(x = power, y = fit, label = paste(sft$powerEstimate)), vjust = -1.5, color = "black") +
  geom_hline(yintercept = 0.90, color = "red") +
  labs(x = "Soft Threshold (power)", y = "Scale Free Topology Model Fit, signed R^2", title = "Scale independence") +
  theme_minimal()


meanconn <- ggplot(df_fit, aes(x = power, y = meanConnectivity)) +
  geom_point(color = "blue") +
  labs(x = "Soft Threshold (power)", y = "Mean Connectivity", title = "Mean connectivity") +
  theme_minimal()

softthresholding_plot + meanconn

ggsave("results/plots/softThresh_curves.pdf", width = 10, height = 6, dpi = 300)

softPower <- sft$powerEstimate
if (is.na(softPower)) softPower <- 6 

# 3. Construct the gene network and identify modules

net <- blockwiseModules(
  datExpr,
  power = softPower,
  networkType = "signed",
  TOMType = "signed",
  corType = "pearson",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 3
)

# Convert labels to colors for plotting
moduleColors <- labels2colors(net$colors)
colors_for_block <- as.matrix(unname(moduleColors[net$blockGenes[[1]]]))

plotDendroAndColors(
  net$dendrograms[[1]],
  colors_for_block,
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  save = TRUE
)


# Grab module eigengenes from the blockwiseModules output
MEs <- net$MEs
colnames(MEs) <- paste0("ME", labels2colors(as.numeric(gsub("ME", "", colnames(MEs)))))

# Binarize clinical traits correctly
sample_labels <- as.factor(sample_labels)
traitData <- as.data.frame(model.matrix(~ 0 + sample_labels))
colnames(traitData) <- gsub("^sample_labels", "", colnames(traitData))
rownames(traitData) <- rownames(datExpr)

moduleTraitCor <- cor(MEs, traitData, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# Prepare data for ggplot2 multi-trait heatmap
cor_df <- as.data.frame(moduleTraitCor) %>% 
  rownames_to_column(var = "Module") %>%
  pivot_longer(cols = -Module, names_to = "Trait", values_to = "Correlation")

pval_df <- as.data.frame(moduleTraitPvalue) %>% 
  rownames_to_column(var = "Module") %>%
  pivot_longer(cols = -Module, names_to = "Trait", values_to = "Pvalue")

plot_df <- left_join(cor_df, pval_df, by = c("Module", "Trait")) %>%
  mutate(
    Module = gsub("ME", "", Module),
    Significance = paste0(round(Correlation, 2), "\n(", signif(Pvalue, 2), ")")
  )

ggplot(plot_df, aes(x = Trait, y = Module, fill = Correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Significance), size = 3, color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1), name="Pearson r") +
  labs(title = "Module-Trait Relationships across Groups", x = "Clinical Phenotype / Class", y = "WGCNA Module") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold")
  )

ggsave("results/plots/module_heatmap.pdf", width = 12, height = 10, dpi = 300)

# 4. Extract Gene Significance and Module Membership
geneTraitSignificance <- as.data.frame(cor(datExpr, traitData, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(datExpr)))

names(geneTraitSignificance) <- paste0("GS.", names(traitData))
names(GSPvalue) <- paste0("p.GS.", names(traitData))

modNames <- substring(names(MEs), 3)
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(datExpr)))

names(geneModuleMembership) <- paste0("MM.", modNames)
names(MMPvalue) <- paste0("p.MM.", modNames)

# 5. Hypergeometric Overrepresentation Analysis (ML Consensus vs WGCNA Modules)
universe_genes <- colnames(datExpr)
ml_query_genes <- consensus_DEGs$gene

gene_module_df <- data.frame(
  gene = universe_genes,
  Module = moduleColors
)

unique_modules <- unique(moduleColors)
unique_modules <- setdiff(unique_modules, "grey")

enrichment_results <- data.frame(
  Module = character(), 
  Pvalue = numeric(), 
  OddsRatio = numeric(), 
  OverlapCount = numeric(),
  ModuleSize = numeric(),
  stringsAsFactors = FALSE
)

total_universe_size <- length(universe_genes)
total_query_size <- length(intersect(ml_query_genes, universe_genes))

for (mod in unique_modules) {
  mod_genes <- gene_module_df$gene[gene_module_df$Module == mod]
  module_size <- length(mod_genes)
  
  overlap_genes <- intersect(mod_genes, ml_query_genes)
  a <- length(overlap_genes)
  b <- module_size - a
  c <- total_query_size - a
  d <- total_universe_size - module_size - c
  
  mat <- matrix(c(a, c, b, d), nrow = 2)
  ft <- fisher.test(mat, alternative = "greater")
  
  enrichment_results <- rbind(enrichment_results, data.frame(
    Module = mod,
    Pvalue = ft$p.value,
    OddsRatio = as.numeric(ft$estimate),
    OverlapCount = a,
    ModuleSize = module_size
  ))
}

enrichment_results$FDR <- p.adjust(enrichment_results$Pvalue, method = "BH")
enrichment_results <- enrichment_results %>% arrange(FDR)
print(enrichment_results)

# 6. Filter Significant Modules and Extract High-Confidence Hubs
significant_modules <- enrichment_results$Module[enrichment_results$FDR < 0.05]

hub_threshold_list <- c()

for (mod in significant_modules) {
  mm_col <- paste0("MM.", mod)
  gs_col <- "GS.Primary Tumor"
  
  if (mm_col %in% names(geneModuleMembership) && gs_col %in% names(geneTraitSignificance)) {
    module_hubs <- rownames(geneModuleMembership)[
      moduleColors == mod &
      abs(geneModuleMembership[[mm_col]]) > 0.80 &   
      abs(geneTraitSignificance[[gs_col]]) > 0.30     
    ]
    hub_threshold_list <- c(hub_threshold_list, module_hubs)
  }
}

selected_hub_genes <- unique(hub_threshold_list)

selected_hub_genes

# Save the finalized, filtered hub gene list for downstream enrichment analysis
saveRDS(selected_hub_genes, file = "selected_hub_genes.rds")