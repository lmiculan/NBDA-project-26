# Canonical EA with GO and KEGG
# NBDA project - Leonardo Miculan, mat. 256108
library(clusterProfiler)
library(ggplot2)
library(org.Hs.eg.db)
library(enrichplot)

rm(list = ls())
selected_hub_genes <- readRDS("selected_hub_genes.rds")
load("filtered_vst_matrix.rData")

selected_hub_genes_clean <- sub("\\..*$", "", selected_hub_genes)
universe_genes_clean <- sub("\\..*$", "", rownames(filtered_vst_matrix))

selected_hub_genes_symbol <- bitr(
  selected_hub_genes_clean,
  fromType = "ENSEMBL",
  toType = "SYMBOL",
  OrgDb = org.Hs.eg.db
)

universe_genes_symbol <- bitr(
  universe_genes_clean,
  fromType = "ENSEMBL",
  toType = "SYMBOL",
  OrgDb = org.Hs.eg.db
)

save(list = c("selected_hub_genes_symbol"), file = "EA_input_data.rData")

# GO Enrichment Analysis
go_enrichment_bp <- enrichGO(gene = selected_hub_genes_clean,
                          keyType = "ENSEMBL",
                          OrgDb = org.Hs.eg.db,
                          universe = universe_genes_clean,
                          ont = "BP",  # Biological Process
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2,
                          readable = TRUE)

go_enrichment_mf <- enrichGO(gene = selected_hub_genes_clean,
                          keyType = "ENSEMBL",
                          universe = universe_genes_clean,
                          OrgDb = org.Hs.eg.db,
                          ont = "MF",  # Molecular Function
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2,
                          readable = TRUE)

go_enrichment_cc <- enrichGO(gene = selected_hub_genes_clean,
                          universe = universe_genes_clean,
                          keyType = "ENSEMBL",
                          OrgDb = org.Hs.eg.db,
                          ont = "CC",  # Cellular Component
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2,
                          readable = TRUE)

# KEGG Pathway Enrichment Analysis

selected_hub_genes_entrez <- bitr(selected_hub_genes_clean, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
universe_genes_entrez <- bitr(universe_genes_clean, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

kegg_enrichment <- enrichKEGG(gene = selected_hub_genes_entrez$ENTREZID,
                              universe = universe_genes_entrez$ENTREZID,
                              organism = 'hsa',
                              pAdjustMethod = "BH",
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.2)

# Visualize the results
# GO Enrichment plots
go_bp_plot <- enrichplot::dotplot(go_enrichment_bp) + ggtitle("GO Biological Process Enrichment")
ggsave("results/plots/GO_bp.pdf", plot = go_bp_plot, width = 7, height = 6, dpi = 300)
go_mf_plot <- enrichplot::dotplot(go_enrichment_mf) + ggtitle("GO Molecular Function Enrichment")
ggsave("results/plots/GO_mf.pdf", plot = go_mf_plot, width = 7, height = 6, dpi = 300)
go_cc_plot <- enrichplot::dotplot(go_enrichment_cc) + ggtitle("GO Cellular Component Enrichment")
ggsave("results/plots/GO_cc.pdf", plot = go_cc_plot, width = 7, height = 6, dpi = 300)

# KEGG Enrichment plot
kegg_plot <- enrichplot::dotplot(kegg_enrichment) + ggtitle("KEGG Pathway Enrichment")
ggsave("results/plots/GO_kegg.pdf", plot = kegg_plot, width = 7, height = 6, dpi = 300)
