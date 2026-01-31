#Anomaly Detection (Isolation Forest)
library(tidyverse)
library(isotree)

#loading PCA results
pca <- prcomp(t(log_cpm), scale. = TRUE)

pca_df <- as.data.frame(pca$x[, 1:10])
pca_df$SAMPID <- colnames(log_cpm)

#Running isolation forest
set.seed(123)

iso_model <- isolation.forest(
  pca_df %>% select(starts_with("PC")),
  ntrees = 500,
  sample_size = 256
)

pca_df$if_score <- predict(iso_model, newdata = pca_df, type = "score")

#setting threshold
threshold <- quantile(pca_df$if_score, 0.95)

pca_df <- pca_df %>%
  mutate(
    anomaly = if_else(if_score >= threshold, "Outlier", "Typical")
  )

table(pca_df$anomaly)

#visualizing anomalies on PCA
pca_plot_df <- pca_df %>%
  left_join(cortex_expr_meta, by = "SAMPID")

ggplot(pca_plot_df, aes(PC1, PC2, color = anomaly)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = c("Typical" = "grey70", "Outlier" = "red")) +
  theme_classic(base_size = 14) +
  labs(
    title = "Isolation Forest–Detected Anomalies",
    subtitle = "Applied to PCA space (PC1–PC10)",
    color = "Sample type"
  )

saveRDS(pca_df, "results/ml/isolation_forest_scores.rds")

ggsave("results/pca/isolation_forest_pca.png", width = 6, height = 5, dpi = 300)
