library(tidyverse)
source("scripts/functions/extract_dataset.R")

# Read predictions
pred_met <- extract_dataset("1kg", readRDS("data/performance_tools.rds"))

# Calculate additional score.
# Using terminology of paper: accuracy = (1 - undecided_ratio) * correctness
# (Note that the column "accuracy" below corresponds with "correctness")
tools_performance_pred <- pred_met %>%
  mutate(S=(1-undecided_ratio) * accuracy) %>%
  arrange(data_type, gene, desc(S)) %>%
  select(data_type, gene, tool, S)

# Average performance of tool per gene
avg_tools_performance <- tools_performance_pred %>%
  group_by(data_type, tool) %>%
  summarise(S=mean(S, na.rm = T), .groups="drop")

# For DPA1 and DPB1 we have no accuracy,
# assume the tools would perform as for an average gene
tools_performance <- list("DPA1", "DPB1") %>%
  map( ~ mutate(avg_tools_performance, gene = .) %>%
         # Only add scores for tools that actually support these genes
         filter((gene != "DPA1") |
                  !(
                    tool %in% c("Optitype", "Polysolver", "PHLAT", "xHLA")
                  )) %>%
         filter((gene != "DPB1") |
                  !(
                    tool %in% c("Optitype", "Polysolver", "PHLAT")
                  ))) %>%
  bind_rows %>%
  bind_rows(tools_performance_pred) %>%
  drop_na %>%
  mutate(mhc_class = if_else(gene %in% c("A", "B", "C"), "I", "II"))

# Calculate mean score of tool for the current MHC class and data type
# Ignore unsupported genes
# -> Will still produce NaN when every gene of this MHC class is not supported
# Not in tie_solver -> score 0
# (means that every gene in this MHC class is not supported)
tie_breaker <- tools_performance %>%
  group_by(data_type, mhc_class, tool) %>%
  mutate(S_class=coalesce(mean(S, na.rm=T), 0))

saveRDS(tie_breaker, "data/metaclassifier_dependencies/tie_solver.rds")