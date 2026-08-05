# Unsupervised learning
# NBDA project - Leonardo Miculan, mat. 256108

rm(list = ls())
load("filtered_vst_matrix.rData")
library(ClusterR)

# Hierarchical clustering
# Compute the distance matrix
dist_matrix <- dist(t(filtered_vst_matrix), method = "euclidean")
# Perform hierarchical clustering
hc <- hclust(dist_matrix, method = "complete")

# Plot the dendrogram
pdf("results/plots/Dendrogram_TCGA-PRAD.pdf", width = 14, height = 8)
plot(hc, main = "Hierarchical Clustering Dendrogram of TCGA-PRAD Samples", 
     xlab = "", sub = "", cex = 0.6, labels = colDataDF_down$external_id)
dev.off()

cluster_assignments <- cutree(hc, k = 2)  # Cut the dendrogram into 2 clusters

# PCA visulalization of the clusters
ggplot(data = as.data.frame(pca_data$x), aes(x = PC1, y = PC2)) +
  geom_point(aes(color = as.factor(cluster_assignments)), size = 3) +
  labs(title = "PCA of TCGA-PRAD Samples with Hierarchical Clustering", x = paste("PC1 (", round(summary(pca_data)$importance[2, 1] * 100, 2), "%)", sep = ""), y = paste("PC2 (", round(summary(pca_data)$importance[2, 2] * 100, 2), "%)", sep = "")) +
  theme_minimal() +
  scale_color_discrete(name = "Cluster", labels = c("Cluster 1", "Cluster 2"))

ggsave("results/plots/PCA_Hierarchical_Clustering_TCGA-PRAD_Samples.pdf", width = 8, height = 6, dpi = 300)
