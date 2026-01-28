#loading library
library(tidyverse)

#loading metadata
sample_meta<-read_tsv("C:/Users/aadirakr2020/Desktop/Audit/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", show_col_types=FALSE)
subject_meta<-read_tsv("C:/Users/aadirakr2020/Desktop/Audit/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt", show_col_types=FALSE)

#Filtering brain-cortex samples
cortex_samples <- sample_meta %>%
  filter(SMTSD == "Brain - Cortex")

stopifnot(nrow(cortex_samples) > 0)

#Deriving SAMPID from SUBJID
cortex_samples <- cortex_samples %>%
  mutate( SUBJID = sub("(^GTEX-[^-]+).*", "\\1", SAMPID))

#Joining donor metadata (sex,age)
cortex_meta<-cortex_samples %>%
  left_join(subject_meta, by="SUBJID")

#Confirmation
colnames(cortex_meta)
table(cortex_meta$SEX)
length(unique(cortex_meta$SUBJID))

#Removing samples with missing sex
cortex_meta<-cortex_meta %>%
  filter(!is.na(SEX))

#Checking number of cortex samples
nrow(cortex_meta)

# # Randomly selected 70 Brain Cortex samples for focused analysis
set.seed(123)
cortex_70<-cortex_meta %>%
  sample_n(70)

#Saving filtered metadata
write_tsv(cortex_70, "C:/Users/aadirakr2020/Desktop/Audit/GTEx_Brain_Cortex_metadata_70samples.tsv")
