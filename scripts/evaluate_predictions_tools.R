library(tidyverse)

source("scripts/functions/evaluate_predictions.R")

gold_standard_files = c(
  "1kg" = "1kg",
  "NCI-60" = "nci60"
)

# Read nested lists with all predictions
predictions_all <- readRDS("data/hla_predictions_ggroup.rds")

# Evaluate predictions for which we have a gold standard
res <- imap(gold_standard_files, ~process_dataset(.y, .x, predictions_all))

# Extract the is_correct / metrics data frames
evaluated_lst <- res %>%
  map_depth(2, 1)

metrics_lst <- res %>%
  map_depth(2, 2)

# Save to RDS files
saveRDS(evaluated_lst, "data/iscorrect_tools.rds")
saveRDS(metrics_lst, "data/performance_tools.rds")