#PCA and UMAP
library(tidyverse)
library(ggplot2)
library(uwot)

log_cpm <- readRDS("data/processed/logCPM_cortex_60.rds")
meta <- readRDS("data/processed/metadata_cortex_60.rds")

pca <- prcomp(t(log_cpm), center = TRUE, scale. = FALSE)

pca_df <- as.data.frame(pca$x[, 1:5])
pca_df$SAMPID <- rownames(pca_df)
pca_df <- left_join(pca_df, meta, by = "SAMPID")

#PCA for sex
ggplot(pca_df, aes(PC1, PC2, color = factor(SEX))) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "PCA of GTEx Brain Cortex RNA-seq", color = "Sex") +
  theme_classic()
ggsave("results/pca/PCA_sex.png", width = 6, height = 5)

#PCA for donor effect
ggplot(pca_df, aes(PC1, PC2, color = SUBJID)) +
  geom_point(size = 3) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(title = "PCA colored by donor")
ggsave("results/pca/PCA_donor.png", width = 6, height = 5)

#Running UMAP on PCA space
set.seed(123)

umap_coords <- umap(pca$x[, 1:10])
umap_df <- as.data.frame(umap_coords)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$SAMPID <- rownames(pca$x)
umap_df <- left_join(umap_df, meta, by = "SAMPID")

#UMAP plot for sex
ggplot(umap_df, aes(UMAP1, UMAP2, color = factor(SEX))) +
  geom_point(size = 3, alpha = 0.8) +
  theme_classic() +
  labs(title = "UMAP of Brain Cortex Samples", color = "Sex")
ggsave("results/pca/UMAP_sex.png", width = 6, height = 5)

#UMAP plot for donor eefect (exploratory)
ggplot(umap_df, aes(UMAP1, UMAP2, color = SUBJID)) +
  geom_point(size = 3) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(title = "UMAP colored by donor (exploratory)")
ggsave("results/pca/UMAP_donor.png", width = 6, height = 5)
