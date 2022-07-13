library(tidyverse)
# Some files on the GDC portal are not available anymore
# Combine metadata of both datasets

metadata_old = read.table("downloads/TCGA/HLA/metadata_202001.tsv",
                          sep = "\t",
                          header = T)

metadata_new = read.table("downloads/TCGA/HLA/metadata.tsv",
                      sep = "\t",
                      header = T)

conversion_table <- list(old=metadata_old, new=metadata_new) %>%
  map_dfr(function(df) {
    select(df,
           file_id,
           data_type = experimental_strategy,
           sample_id = cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id)
  }, .id = "version")

conversion_version <- conversion_table %>%
  group_by(file_id, data_type, sample_id) %>%
  summarise(v = list(version))

# Check whether file ID -> (data_type, aliquot ID) mapping is unique
conversion_version %>%
  group_by(file_id) %>%
  add_count() %>%
  filter(n > 1) %>%
  nrow
# => OK

# As mapping does not conflict for the two versions:
# Take union of two versions, remove version information
conversion_version %>%
  select(-v) %>%
  saveRDS("downloads/TCGA/HLA/file_ids_to_aliquots.rds")
