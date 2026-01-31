#what do outliers correlate with?
outliers <- pca_df %>%
  filter(anomaly == "Outlier") %>%
  select(SAMPID, starts_with("PC"), if_score)

outliers

#correlation between IF score and cell type
cor.test(pca_df$if_score, pca_df$astro_score, method = "spearman")
cor.test(pca_df$if_score, pca_df$neuron_score, method = "spearman")

#plotting IF score and cell type
ggplot(pca_df, aes(astro_score, if_score, color = anomaly)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(x = "Astro score", y = "Isolation Forest score")
ggsave("results/pca/astro_vs_ifscore")

ggplot(pca_df, aes(neuron_score, if_score, color = anomaly)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(x = "Neuron score", y = "Isolation Forest score")
ggsave("results/pca/neuron_vs_ifscore")

#Wilcox test to see if correlation is statistically significant
wilcox.test(pca_df$astro_score[pca_df$anomaly == "Typical"],
            pca_df$astro_score[pca_df$anomaly == "Outlier"], exact = FALSE)
wilcox.test(pca_df$neuron_score[pca_df$anomaly == "Typical"],
            pca_df$neuron_score[pca_df$anomaly == "Outlier"], exact = FALSE)
