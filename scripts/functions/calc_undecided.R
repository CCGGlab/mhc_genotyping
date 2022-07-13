# Calculate the number of undecided predictions
calc_undecided <- function(pred) {
  pred %>%
    # Average number of undecided per gene (of the total number of predictions made)
    # (Note: This also equals the average of (sum(is.na(allele)) per sample, gene) / 2)
    group_by(gene, .add = TRUE) %>%
    summarise(undecided_ratio=mean(is.na(allele)), .groups="drop")
}