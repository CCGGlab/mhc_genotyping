library(tidyverse)
# Read gold standard data to match the cell line IDs
gs <- readRDS("downloads/pub/adams_2005/hla_calls.rds")

# ----
# Abaan (DNA)
run_info <- readr::read_csv("downloads/pub/abaan_2013/run_info.csv")

# Determine overlap between cell lines of prediction and gold standard
setdiff(run_info$SampleName, gs$`Cell Line`)
# => 15 cell lines of the prediction were not present in the gold standard (with an exact match)

# Some cell lines are present with a slightly different name
matched_sample_ids <- fuzzyjoin::stringdist_right_join(
  distinct(run_info, sample_id.old = SampleName),
  distinct(gs, sample_id.gs = `Cell Line`),
  by = 1,
  distance_col = "dist"
) %>%
  group_by(sample_id.old) %>%
  slice_min(order_by = dist, n = 1)

# Matching of sample IDs was checked manually and considered to be correct
matched_sample_ids %>%
  filter(dist > 0)
# For all 58 cell lines present in the gold standard data a match was found
# (Remove the dist column as this is annoying when matching on ID)
matched_sample_ids <- select(matched_sample_ids, -dist)

# Create the final mapping table from SRR accession number
# to the cell line names used in the gold standard
mapping_table_abaan <- select(run_info, sample_id = Run, sample_id.old = SampleName) %>%
  left_join(matched_sample_ids, by = "sample_id.old") %>%
  transmute(sample_id.srr = sample_id, sample_id.new = if_else(is.na(sample_id.gs), sample_id.old, sample_id.gs))

# ----
# Idem for Reinhold (RNA)
reinhold_meta <- read_delim(
  "downloads/pub/reinhold_2019/reinhold_meta.txt",
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

setdiff(reinhold_meta$library_name, gs$`Cell Line`) %>% length
# => 16 cell lines of the prediction were not present in the gold standard (with an exact match)

# Some cell lines are present with a slightly different name
matched_sample_ids <- fuzzyjoin::stringdist_right_join(
  distinct(reinhold_meta, sample_id.old = library_name),
  distinct(gs, sample_id.gs = `Cell Line`),
  by = 1,
  distance_col = "dist"
) %>%
  group_by(sample_id.old) %>%
  slice_min(order_by = dist, n = 1)

# Matching of sample IDs was checked manually and considered to be correct
matched_sample_ids %>%
  filter(dist > 0)

matched_sample_ids <- select(matched_sample_ids, -dist)

# Create the final mapping table from SRR accession number
# to the cell line names used in the gold standard
mapping_table_reinhold <- select(reinhold_meta, sample_id = run_accession, sample_id.old = library_name) %>%
  left_join(matched_sample_ids, by = "sample_id.old") %>%
  transmute(sample_id.srr = sample_id, sample_id.new = if_else(is.na(sample_id.gs), sample_id.old, sample_id.gs))

mapping_table <- bind_rows(mapping_table_abaan, mapping_table_reinhold)
# ----

saveRDS(mapping_table, "data/sample_names_nci60.rds")