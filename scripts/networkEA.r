# Network enrichment analysis data preparation for STRING and EnrichNet analysis via Cytoscape
# NBDA project - Leonardo Miculan, mat. 256108

library(STRINGdb)
library(enrichR)

# Load the list of hub genes identified from WGCNA
rm(list = ls())
load("EA_input_data.rData")  # Load the previously saved mappings for hub and universe genes

# Write tsv gene list file for Cytoscape STRINGapp
write.table(selected_hub_genes_symbol$SYMBOL, file = "results/hub_genes.tsv", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

# STRING PPI Analysis done in Cytoscape using the STRINGapp.

# EnrichNet Analysis
enrichr_results <- enrichr(selected_hub_genes_symbol$SYMBOL, databases = c("GO_Biological_Process_2026","GO_Molecular_Function_2026", "GO_Cellular_Component_2026", "KEGG_2026"))

# Function to transform enrichR dataframe to EnrichmentMap Generic Format
format_enrichr_for_em <- function(enrichr_df, phenotype = 1) {
  em_df <- enrichr_df %>%
    # Filter for significant terms
    filter(Adjusted.P.value < 0.05) %>%
    # Select and rename columns matching Cytoscape EnrichmentMap specifications
    transmute(
      GeneSet     = Term,
      Description = Term,
      p.value     = P.value,
      FDR         = Adjusted.P.value,
      Phenotype   = phenotype,
      # Convert enrichR ';' gene separator to commas or tabs for EnrichmentMap parsing
      Genes       = gsub(";", ",", Genes)
    )
  
  return(em_df)
}

# 1. Format individual enrichR database results
em_go_bp <- format_enrichr_for_em(enrichr_results$GO_Biological_Process_2026, phenotype = 1)
em_go_mf <- format_enrichr_for_em(enrichr_results$GO_Molecular_Function_2026, phenotype = 1)
em_go_cc <- format_enrichr_for_em(enrichr_results$GO_Cellular_Component_2026, phenotype = 1)
em_kegg  <- format_enrichr_for_em(enrichr_results$KEGG_2026, phenotype = 1)

# 2. Export reformatted data to TSV files
write.table(em_go_bp, file = "results/GO_BP_EnrichmentMap_input.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(em_go_mf, file = "results/GO_MF_EnrichmentMap_input.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(em_go_cc, file = "results/GO_CC_EnrichmentMap_input.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(em_kegg,  file = "results/KEGG_EnrichmentMap_input.txt",  sep = "\t", quote = FALSE, row.names = FALSE)
