library(tidyverse)
source("scripts/functions/ggroup_mapper.R")

readRDS("downloads/pub/adams_2005/hla_calls.rds") %>%
  rename(sample_id = `Cell Line`) %>%
  select(-ID, -Tissue) %>%
  allele_map %>%
  saveRDS(str_glue("data/gold_standard_nci60.rds"))
