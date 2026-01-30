#Overlaying marker expression on existing PCA/UMAP
markers <- c(
  "RBFOX3","MAP2","TUBB3",
  "GFAP","AQP4",
  "GAD1","GAD2","PVALB","SST","VIP",
  "CUX2","SATB2",
  "BCL11B","FEZF2","FOXP2"
)

marker_expr <- log_cpm[rownames(log_cpm) %in% markers, ]

#Confirmation
# Are marker genes present?
head(rownames(log_cpm))
rownames(log_cpm) <- gsub("\\..*$", "", rownames(log_cpm))

library(org.Hs.eg.db)
library(AnnotationDbi)

marker_genes <- c(
  "RBFOX3","MAP2","TUBB3","GFAP","AQP4",
  "GAD1","GAD2","PVALB","SST","VIP",
  "CUX2","SATB2","BCL11B","FEZF2","FOXP2"
)

marker_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = marker_genes,
  keytype = "SYMBOL",
  columns = c("SYMBOL", "ENSEMBL")
)

marker_map <- marker_map %>%
  filter(!is.na(ENSEMBL)) %>%
  distinct(ENSEMBL, .keep_all = TRUE)

present_markers <- intersect(marker_map$ENSEMBL, rownames(log_cpm))
present_markers
length(present_markers)

marker_expr <- log_cpm[present_markers, ]
dim(marker_expr)

pca_markers <- prcomp(t(marker_expr), scale. = TRUE)

pca_df <- data.frame(
  PC1 = pca_markers$x[,1],
  PC2 = pca_markers$x[,2],
  SEX = factor(cortex_expr_meta$SEX),
  DONOR = cortex_expr_meta$SUBJID)

library(ggplot2)

ggplot(pca_df, aes(PC1, PC2, color = SEX)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCA of Cell-Type Marker Expression", subtitle = "GTEx Brain Cortex Samples", color = "Sex")
ggsave("results/pca/PCA_cell_comp_sex.png", width=6, height=5)

ggplot(pca_df, aes(PC1, PC2, color = DONOR)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "PCA of Cell-Type Marker Expression", subtitle = "GTEx Brain Cortex Samples", color = "Donor")
ggsave("results/pca/PCA_cell_comp_donor.png", width=6, height=5)


#UMAP for cell composition (exploratory)
marker_mat_scaled <- scale(marker_mat)
library(uwot)
set.seed(123)

umap_res <- umap(
  marker_mat_scaled,
  n_neighbors = 15,   
  min_dist = 0.3,
  metric = "euclidean")

umap_df <- data.frame(
  UMAP1 = umap_res[,1],
  UMAP2 = umap_res[,2]
)

umap_df <- cbind(umap_df, cortex_meta)

library(ggplot2)

ggplot(umap_df, aes(UMAP1, UMAP2, color = factor(SEX))) +
  geom_point(size = 3, alpha = 0.85) +
  labs(
    title = "UMAP of Cell-Type Marker Expression",
    subtitle = "GTEx Brain Cortex Samples",
    color = "Sex"
  ) +
  theme_classic()
ggsave("results/pca/UMAP_cell_comp_sex.png", width=6, height=5)

ggplot(umap_df, aes(UMAP1, UMAP2, color = SUBJID)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(
    title = "UMAP of Cell-Type Marker Expression",
    subtitle = "Colored by Donor",
    color = "Donor",
  ) +
  theme_classic() +
  theme(legend.position="None")
ggsave("results/pca/UMAP_cell_comp_donor.png", width=6, height=5)
