library(tidyverse)

source("scripts/functions/extract_dataset.R")

# Input: all alleles predicted for a particular (data_type, sample_id, gene)
meta_prediction <- function(df, tie_solver, tie_group, selected_tools=NULL) {
  tscore_df <- df %>%
    # Only consider a subset of tools
    filter(is.null(selected_tools) | tool %in% selected_tools) %>%
    group_by(data_type, tool, gene, sample_id) %>%
    # Sort alleles: fix order of allele 1 and 2
    arrange(prediction) %>%
    # Alleles to columns
    mutate(idx=paste0("a", row_number())) %>%
    spread(idx, prediction) %>%
    # Do not consider predictions containing an NA
    drop_na
  
  # TODO: Check what happens when no prediction was made for a sample
  # by any tool
  typing <- tscore_df %>%
    # Count how many tools voted for a certain genotype (per sample and gene)
    # (Also include mhc_class as this column needs to be kept in the data)
    group_by(data_type, sample_id, gene, mhc_class, a1, a2) %>%
    summarise(value=n(), tools=list(tool), .groups="keep") %>%
    # Select allele with highest typing score
    # (We need to drop the alleles from the grouping first)
    group_by(data_type, sample_id, gene) %>%
    slice_max(value, n=1, with_ties=T) %>%
    # In case of ties, take the prediction of the best tool in this list
    # Merge with the tool list to obtain the score for the current data_type, gene
    unnest(tools) %>%
    inner_join(tie_solver, by=c("data_type", tie_group, "tools"="tool")) %>%
    slice_max(S, n=1, with_ties=T) %>%
    # If there is still a tie, select the prediction made by the tool
    # with the highest average score for the entire MHC class
    slice_max(S_class, n=1, with_ties=T)
  
  typing <- typing %>%
    # In case there are still ties, take a random prediction
    slice_sample(n = 1) %>%
    ungroup %>%
    # Rename the allele columns to get a uniform format
    select(data_type, sample_id, gene, `1` = a1, `2` = a2)
  
  typing
}
