library(tidyverse)

mapping = readRDS("data/ggroup_mapping.rds") %>%
  select(allele, group) %>%
  deframe

ga_map <- function(gene, allele) {
  gene_allele = paste0(gene, '*', allele)
  
  # Return allele itself if no G-group was defined
  if_else(gene_allele %in% names(mapping), str_remove(mapping[gene_allele], '^[^*]*\\*'), allele)
}

allele_col_map <- function(x) {
  # This function receives one column at a time
  # Retrieve gene name from column name:
  gene_idx = simplify2array(str_split(cur_column(), "\\."))
  gene = gene_idx[[1,1]]
  # Perform the G-group mapping on this column
  ga_map(gene, x)
}

allele_map <- function(prediction) {
  mutate(prediction, across(-sample_id, allele_col_map))
}