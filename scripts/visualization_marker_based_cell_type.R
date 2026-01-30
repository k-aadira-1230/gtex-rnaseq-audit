# Marker-based cell-type gradient visualization
# Maps Ensembl IDs to gene symbols
ens_to_symbol <- symbol_to_ens %>%
  filter(ENSEMBL %in% present_ens) %>%
  distinct(ENSEMBL, SYMBOL)

marker_mat_renamed <- marker_mat[, present_ens]
colnames(marker_mat_renamed) <- ens_to_symbol$SYMBOL[
  match(colnames(marker_mat_renamed), ens_to_symbol$ENSEMBL)
]

# Computes astrocyte and neuronal scores
marker_scaled <- scale(marker_mat_renamed)

umap_df$astro_score  <- rowMeans(marker_scaled[, c("GFAP", "AQP4")])
umap_df$neuron_score <- rowMeans(marker_scaled[, c("RBFOX3", "MAP2", "TUBB3")])

# Visualizes gradients on UMAP and PCA
ggplot(umap_df, aes(UMAP1, UMAP2, color = astro_score)) +
  geom_point(size = 3) +
  scale_color_viridis_c() +
  theme_classic() +
  labs(title = "Astrocyte Gradient on Cell-Type UMAP")
ggsave("results/pca/Astrocytegradient.png", width=6, height=5)

ggplot(umap_df, aes(UMAP1, UMAP2, color = neuron_score)) +
  geom_point(size = 3) +
  scale_color_viridis_c() +
  theme_classic() +
  labs(title = "Neuronal Gradient on Cell-Type UMAP")
ggsave("results/pca/Neuronalgradient.png", width=6, height=5)

# Quantifies association with major PCs
cor(merged_df$astro_score, merged_df$PC1)
cor(merged_df$astro_score, merged_df$PC2)

#PCA-Astro score

library(ggplot2)

p_astro_pca <- ggplot(merged_df, aes(x = PC1, y = PC2, color = astro_score)) +
  geom_point(size = 3, alpha = 0.85) +
  scale_color_viridis_c(name = "Astrocyte score") +
  theme_classic(base_size = 14) +
  labs(
    title = "PCA of Brain Cortex Samples",
    subtitle = "Color indicates astrocyte marker expression",
    x = "PC1",
    y = "PC2"
  )

p_astro_pca

ggsave(
  filename = "results/pca/PCA_colored_by_astro_score.png",
  plot = p_astro_pca,
  width = 6,
  height = 5,
  dpi = 300
)
