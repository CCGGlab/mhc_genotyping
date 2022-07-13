library(tidyverse)
source("scripts/functions/extract_dataset.R")

# => Two first selected tools for optimization of metaclassifier

# Read data frames
performance_tools <- extract_dataset("1kg", readRDS("data/performance_tools.rds")) %>%
  mutate(mhc_class = if_else(gene %in% c("A", "B", "C"), "I", "II")) %>%
  mutate(S = (1 - undecided_ratio) * accuracy) %>%
  group_by(data_type, tool, mhc_class) %>%
  mutate(S_class = coalesce(mean(S, na.rm=T), 0)) %>%
  ungroup

is_correct <- extract_dataset("1kg", readRDS("data/iscorrect_tools.rds"))

tools <- distinct(performance_tools, data_type, tool)

determine_initial_tools <- function(group_var, excluded_tools) {
  # Calculate average performance per gene / MHC class and select the best tool
  best_tools <- performance_tools %>%
    anti_join(excluded_tools, by=c("data_type", "tool")) %>%
    group_by(data_type, across(all_of(group_var)), tool) %>%
    summarise(across(c(undecided_ratio, accuracy, S, S_class), mean, na.rm=T),
              .groups="drop_last") %>%
    slice_max(S) %>%
    # If there are ties, select the tool with the highest average performance
    slice_max(S_class) %>%
    # Map of (data_type, mhc_class) -> best tool
    select(data_type, tool, !!group_var)
  
  # Number of correct allele predictions by tool for sample and gene
  # -> 2 is a perfect score for that sample and gene
  num_correct <- is_correct %>%
    anti_join(excluded_tools, by=c("data_type", "tool")) %>%
    group_by(data_type, tool, sample_id, gene) %>%
    summarise(n = sum(is_correct, na.rm = T), .groups="drop") %>%
    mutate(mhc_class = if_else(gene %in% c("A", "B", "C"), "I", "II"))
  
  # Take the best tool and the tool that is the most complementary with it
  result <- num_correct %>%
    # Cross the data frame of the best tools with all other tools
    semi_join(best_tools, by = c("data_type", "tool", group_var)) %>%
    left_join(num_correct, by = c("data_type", "sample_id", group_var)) %>%
    # We do not want to compare tools with themselves
    filter(tool.x != tool.y) %>%
    # Find the most complementary tool of the best tool
    # - Which tool made 2 correct predictions on the samples where the best tool was not perfect
    mutate(is_better = ((n.x < 2) & (n.y == 2))) %>%
    # - Calculate per gene / MHC class the number of better predictions by the second tool
    group_by(data_type, across(all_of(group_var)), tool.x, tool.y) %>%
    summarise(num_better = sum(is_better), .groups="drop_last") %>%
    # - Select the tool with the highest number of "better predictions"
    slice_max(num_better) %>%
    # - Select the best tool if there is still a tie
    left_join(performance_tools, by=c("data_type", "tool.y" = "tool", group_var)) %>%
    slice_max(S) %>%
    slice_max(S_class) %>%
    # Convert the result to a convenient format
    select(data_type, !!group_var, tool.x, tool.y) %>%
    group_by(data_type, across(all_of(group_var))) %>%
    summarise(tools=list(c_across()), .groups="drop")
  
  result
}

exclusion_dfs = list(others = tibble(data_type = NA_character_, tool = NA_character_),
                     tcga = tibble(data_type = "RNA-Seq", tool = "HLA-HD"))

initial_tools <- list()
for (dataset in names(exclusion_dfs)) {
  for (group_var in c("gene", "mhc_class")) {
    initial_tools[[dataset]][[group_var]] <-
      determine_initial_tools(group_var, exclusion_dfs[[dataset]])
  }
}

saveRDS(initial_tools, "data/metaclassifier_dependencies/initial_tools.rds")
