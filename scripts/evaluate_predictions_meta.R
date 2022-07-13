library(tidyverse)

source("scripts/functions/evaluate_predictions.R")

# Read metapredictions
predictions <- readRDS("data/hla_predictions_ggroup_meta.rds")

# Mapping of data set name to the suffix of the gold standard RDS files
gold_standard_files = c(
  "1kg" = "1kg",
  "NCI-60" = "nci60"
)

# Process both the best compact version
# and the metaprediction using all supported tools
for (version in c("compact", "alltools")) {
  pred_met_fmt <- predictions[[version]]
  resmeta <- imap(gold_standard_files["1kg"], ~process_dataset(.y, .x, pred_met_fmt))
  
  # Extract the is_correct / metrics data frames
  m_evaluated_lst <- resmeta %>%
    map_depth(2, 1)
  
  m_metrics_lst <- resmeta %>%
    map_depth(2, 2)
  
  # Save to RDS files
  saveRDS(m_evaluated_lst, str_glue("data/iscorrect_meta_{version}.rds"))
  saveRDS(m_metrics_lst, str_glue("data/performance_meta_{version}.rds"))
}
