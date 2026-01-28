library(edgeR)

# Creating DGEList
dge <- DGEList(counts = filtered_counts)

# TMM normalization
dge <- calcNormFactors(dge, method = "TMM")

# Log-CPM transformation
log_cpm <- cpm(dge, log = TRUE, prior.count = 1)

# Saving normalized expression
saveRDS(log_cpm, "data/processed/logCPM_cortex_60.rds")
