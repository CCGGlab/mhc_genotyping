library(tidyverse)
# Plot of metaclasifier using all tools
# compared with results of best tool for that gene and data type

# Functions
source("scripts/functions/extract_dataset.R")

# Read accuracies of tools and metaclassifier
tool_acc <- readRDS("data/performance_tools.rds")
meta_acc <- readRDS("data/performance_meta_alltools.rds")


# Get the accuracy of the best performing tool for that gene and data type
best_acc <- extract_dataset("1kg", tool_acc) %>%
  mutate(S = (1 - undecided_ratio) * accuracy) %>%
  group_by(data_type, gene) %>%
  slice_max(S, with_ties = F)

# Construct data frame for plot
# - Calculate the accuracy
# - Add a column with the score of the best tool
# - Format the accuracy as percentages
plt_df <- extract_dataset("1kg", meta_acc) %>%
  select(-tool) %>%
  mutate(S = (1 - undecided_ratio) * accuracy) %>%
  left_join(
    select(best_acc, data_type, gene, S),
    by = c("data_type", "gene"),
    suffix = c(".meta", ".best")
  ) %>%
  pivot_longer(
    starts_with("S."),
    names_to = "tool",
    names_prefix = "S.",
    values_to = "S"
  ) %>%
  # Label that is displayed on top of the bars
  # Round to 0.1%
  mutate(lbl = scales::percent_format(accuracy = 0.1)(S)) %>%
  # For the order of "meta" and "best" in the legend
  mutate(tool = fct_relevel(tool, "meta", "best")) %>%
  mutate(gene = paste0("HLA-", gene))

plt <-
  ggplot(plt_df, aes(
    x = gene,
    y = S,
    fill = tool,
    label = lbl
    # group = tool
  )) +
  ggpubr::theme_pubclean() +
  theme(
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(linetype = "dotted")
  ) +
  # Make barplot:
  # - Two separate bars for meta / best classifier
  # - Take different transparency for meta / best
  # geom_col(aes(alpha = tool), position = position_dodge2()) +
  geom_col(position = position_dodge2()) +
  geom_text(
    aes(y = 0.9, colour = tool),
    position = position_dodge2(width = 0.9),
    vjust = 0,
    hjust = 0.5,
    size = 3
  ) +
  # Two different panels for DNA / RNA
  facet_wrap( ~ data_type) +
  coord_cartesian(ylim = c(0.9, 1)) +
  # Breaks each 5% for the accuracy
  scale_y_continuous(
    "accuracy",
    breaks = seq(0, 1, 0.01),
    # y-axis labels formatted as percentages, round to 1%
    labels = scales::label_percent(accuracy = 1),
    # Reduce space below the plot
    expand = expansion(mult = c(0.01, 0.05))
  ) +
  # Change colour of the labels for meta / best
  scale_colour_manual(values = c("white", "black"))

ggsave(
  "results/figs/supp_metaclassifier_alltools.pdf",
  plt,
  scale = 2,
  width = 21,
  units = "cm"
)
