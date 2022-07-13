library(tidyverse)

exact_genotype_check <- function(pred, gs) {
  a <- list(pred, gs) %>%
    map(function(.) {
      # Convert both to long table format
      # Split column name to gene name and allele index
      pivot_longer(., !sample_id) %>%
        separate(name, c("gene", "idx")) %>%
        group_by(sample_id, gene) %>%
        # Sort alleles
        arrange(value, .by_group = T) %>%
        mutate(idx = seq_along(value)) %>%
        # Drop sample_id, gene pairs where the prediction is NA
        # for any of the 2 alleles
        group_by(sample_id, gene) %>%
        filter(!any(is.na(c_across(value))))
    })
  
  r <- inner_join(a[[1]], a[[2]], by = c("sample_id", "gene", "idx")) %>%
    group_by(sample_id, gene) %>%
    # Sample is correct if both alleles are correct
    mutate(is_correct = all(value.x == value.y)) %>%
    group_by(gene) %>%
    # Calculate accuracy (hard accuracy)
    summarise(H = mean(is_correct))
}
