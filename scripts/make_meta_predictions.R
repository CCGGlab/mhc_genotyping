library(tidyverse)
future::plan(future::multicore)

# Functions
source("scripts/functions/extract_dataset.R")
source("scripts/functions/build_group_arg_list.R")
source("scripts/functions/meta_prediction.R")

tie_solver <- readRDS("data/metaclassifier_dependencies/tie_solver.rds")

# Which tools are selected for which data type and MHC class
tool_selector_compact <- readRDS("data/meta_optimized.rds") %>%
  # use the per MHC class prediction
  chuck("class") %>%
  mutate(num_tools = map_int(tools, length)) %>%
  # 4 tool model preferred based on plot made by plot_metaclassifier.R
  filter(
    #(data_type == "RNA-Seq" & num_tools == 2) |
           (data_type == "DNA-Seq" & mhc_class == "I" & num_tools == 4) |
           (data_type == "DNA-Seq" & mhc_class == "II" & num_tools == 4)) %>%
  distinct(data_type, mhc_class, tools) %>%
  # Convert to nested list: data_type -> MHC class -> tools
  select(data_type, mhc_class, tools) %>%
  nest(data = c(mhc_class, tools)) %>%
  deframe %>%
  map(deframe)

# tool_selector_compact_tcga <- readRDS("data/meta_optimized.rds") %>%
#   # use the per MHC class prediction
#   chuck("class_tcga") %>%
#   mutate(num_tools = map_int(tools, length)) %>%
#   # Manually determined that this was optimal for all data_types and MHC classes
#   filter((data_type == "DNA-Seq" & mhc_class == "I" & num_tools == 7) |
#            (data_type == "DNA-Seq" & mhc_class == "II" & num_tools == 4) |
#            (data_type == "RNA-Seq" & mhc_class == "I" & num_tools == 2) |
#            (data_type == "RNA-Seq" & mhc_class == "II" & num_tools == 4)) %>%
#   distinct(data_type, mhc_class, tools) %>%
#   # Convert to nested list: data_type -> MHC class -> tools
#   select(data_type, mhc_class, tools) %>%
#   nest(data = c(mhc_class, tools)) %>%
#   deframe %>%
#   map(deframe)

# As long as we do not make a compact meta-prediction for RNA we do not need a separate model for TCGA
tool_selector_compact_tcga <- tool_selector_compact

# Wrapper for meta_prediction function
# Select the correct tools and run the metaclassifier
make_meta_prediction <- function(df, key, tie_solver, tool_selector) {
  if (is.null(tool_selector)) {
    selected_tools <- NULL
  } else {
    selected_tools <- tool_selector[[key$data_type]][[key$mhc_class]]
  }
  
  current_meta_pred <- df %>%
    meta_prediction(tie_solver, "gene", selected_tools) %>%
    ungroup
}

# Convert HLA predictions to wide format: column sample ID and column per allele
# Nested per data type with a single item "meta" (the "tool")
# data_type -> tool -> predictions in wide format
format_results <- function(results) {
  results %>%
    gather("idx", "allele", `1`, `2`, -data_type) %>%
    mutate(locus=str_glue("{gene}.{idx}")) %>%
    # Remove index and gene to prevent NA while converting to wide format
    # Also remove data_type:
    # <- would cause issues with the evaluation function
    select(-idx, -gene) %>%
    spread(locus, allele) %>%
    group_by(data_type) %>%
    summarise(data=list(list(meta=cur_data())), .groups="drop_last") %>%
    deframe
}

predictions_all <- readRDS("data/hla_predictions_ggroup.rds")
predictions_meta <- list()

for (ds in c("1kg", "NCI-60", "TCGA", "Liu", "Riaz")) {
  pred_ds <- extract_dataset(ds, predictions_all)
  
  # Convert predictions to long format
  pred_long <- pred_ds %>%
    gather("locus", "prediction", -data_type, -tool, -sample_id) %>%
    separate(locus, c("gene", "idx"))
  
  # Perform metaprediction for both "version":
  # 1) the optimal compact metaclassifier and 2) with all tools
  for (version in c("compact", "alltools")) {
    # Use no tool selector when all tools must be used
    tool_selector <-  switch(version,
                             "compact" = switch(ds,
                                                "TCGA" = tool_selector_compact_tcga,
                                                tool_selector_compact),
                             "alltools" = NULL)
    
    # Build list of arguments to call "make_meta_predictions" in parallel
    # Make compact consensus classification for DNA only
    if (version == "compact") {
      meta_p_args <- pred_long %>%
        mutate(mhc_class = if_else(gene %in% c("A", "B", "C"), "I", "II")) %>%
        # Only consider DNA-Seq
        filter(data_type == "DNA-Seq") %>%
        group_by(data_type, mhc_class) %>%
        build_group_arg_list("df", "key")
    } else {
      meta_p_args <- pred_long %>%
        mutate(mhc_class = if_else(gene %in% c("A", "B", "C"), "I", "II")) %>%
        group_by(data_type, mhc_class) %>%
        build_group_arg_list("df", "key")
    }
    
    # Make meta predictions
    results <- furrr::future_pmap_dfr(meta_p_args, make_meta_prediction, tie_solver, tool_selector)
    # Store in the list, per version (compact / all tools) and data type
    predictions_meta[[version]][[ds]] <- format_results(results)
  }
}

saveRDS(predictions_meta, "data/hla_predictions_ggroup_meta.rds")