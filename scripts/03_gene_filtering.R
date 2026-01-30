library(edgeR)

# Computing CPM
cpm_mat <- cpm(counts_mat)

# Filtering lowly expressed genes
cpm_threshold <- 1
min_samples <- 6

keep_genes <- rowSums(cpm_mat > cpm_threshold) >= min_samples
filtered_counts <- counts_mat[keep_genes, ]

# Saving filtered counts
saveRDS(filtered_counts, "data/processed/filtered_counts_cortex_60.rds")
