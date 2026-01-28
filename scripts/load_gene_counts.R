library(data.table)
library(tidyverse)

# Loading gene counts
counts <- fread("data/raw/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz",
  skip = 2)

# Extracting cortex sample IDs
cortex_ids <- cortex_70$SAMPID
valid_ids <- intersect(cortex_ids, colnames(counts))

# Subset the count table
counts_cortex <- counts %>%
  select(Name, Description, all_of(valid_ids)) %>%
  filter(!grepl("^_", Name))

# Syncing metadata
cortex_expr_meta <- cortex_70 %>%
  filter(SAMPID %in% valid_ids) %>%
  arrange(match(SAMPID, valid_ids))

# Creating count matrix
gene_ids <- counts_cortex$Name
counts_mat <- counts_cortex %>%
  select(-Name, -Description) %>%
  as.matrix()

rownames(counts_mat) <- gene_ids
storage.mode(counts_mat) <- "numeric"

# QC check
summary(colSums(counts_mat))

# Saving processed data
saveRDS(counts_mat, "data/processed/counts_mat_cortex_70.rds")
saveRDS(cortex_expr_meta, "data/processed/metadata_cortex_70.rds")
