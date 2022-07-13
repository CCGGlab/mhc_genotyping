library(tidyverse)
future::plan(future::multicore(workers=144))

# Helper functions
source("scripts/functions/ggroup_mapper.R")

# Constants:
# Path to conversion scripts
# -> based on the work of the design project students of 2020-2021
conversion_folder <- "scripts/functions/Conversion"

# Define which datasets are processed in this script
# - Keys correspond to the suffix of the tool output folders
# - Values correspond to clean names
dataset_names <- c(
  "1kg" = "1kg",
  "tcga" = "TCGA",
  "nci60_1kg" = "NCI-60"
)
dataset_suffix = names(dataset_names)

# Functions:
tool_to_script_prefix <- function(tool_output, mapped_tools) {
  tool <- basename(tool_output)
  # Map tools
  if (tool %in% names(mapped_tools)) {
    tool <- mapped_tools[[tool]]
  }
  tool
}

# Load in the needed conversion scripts
# Store as named list
load_conversion_functions <- function(conversion_folder) {
  list.files(path = conversion_folder) %>%
    as_tibble_col("script") %>%
    mutate(name = str_remove(script, "_.*"),
           fun = map(script, ~source(paste0(conversion_folder, "/", .), local = T)$value)) %>%
    select(name, fun) %>%
    deframe
}

# Function to add an index to the locus name
gene_name <- function(x) {
  x %>%
    as_tibble_col() %>%
    group_by(value) %>%
    mutate(locus = paste0(value, ".", row_number()), value) %>%
    pull(locus)
}

parse_tool_output <- function(result_folder, conversion_functions) {
  tool <- tool_to_script_prefix(result_folder, mapped_tools)
  print(tool)
  if (!tool %in% names(conversion_functions)) {
    print(paste0("Following tool was not recognized: ", tool))
    result <- tibble()
  } else {
    result <- conversion_functions[[tool]](result_folder)
  }
  
  result <- result[order(row.names(result)),]
  result <- as_tibble(result, rownames = "sample_id", .name_repair = gene_name)
  
  # Return as dataframe with new column for the tool
  res <- mutate(result, tool = tool)
  stopifnot(nrow(res) > 0)
  res
}

# Sample names in the 1kg output folders were not always cleaned
# Sample ID always corresponds to the first 7 characters of the extracted file base name
map_sample_name_1kg <- function(res_df) {
  res_df %>%
    mutate(sample_id = str_sub(sample_id, 1, 7))
}

# Map SRR to cell line names
sample_names_nci60 <- readRDS("data/sample_names_nci60.rds")
map_sample_name_nci60 <- function(res_df) {
  res_df <- res_df %>%
    left_join(sample_names_nci60, by = c("sample_id" = "sample_id.srr")) %>%
    mutate(sample_id = coalesce(sample_id.new, sample_id)) %>%
    select(-sample_id.new)
}

map_sample_name_tcga <- function(res_df) {
  # Post-process res_df
  # For xHLA the sample IDs were derived from the sample_id field in the JSON file
  # As a result the sample_id is still erroneous
  res_df = res_df %>%
    group_by(tool) %>%
    group_modify(function(df, tool) {
      if (tool == "xHLA") {
        df = mutate(df, sample_id=str_remove(sample_id, "([^_]*_){3}"))
      }
      df
    })
  
  # Check whether the sample IDs are valid GUIDs
  num_notcorrect <- res_df$sample_id %>%
    keep(~str_length(.) != 36) %>%
    length
  stopifnot(num_notcorrect == 0)
  
  # Now we have correct GUIDs, adapt data type and sample name later
  res_df
}

# Use case submitter IDs for TCGA data
# Also rename "gathered" to the correct data_type
sample_names_tcga = readRDS("data/sample_names_tcga.rds")
map_sample_name_tcga2 <- function(res_df, key) {
  if (key$dataset_suffix == "tcga") {
    res_df %>%
      select(file_id = sample_id, !all_of("data_type")) %>%
      left_join(sample_names_tcga, by="file_id") %>%
      select(-file_id)
  } else {
    res_df
  }
}

# Create map from dataset suffixes to the correct sample renaming function
sample_funmap <- c("1kg" = map_sample_name_1kg,
                   "nci60_1kg" = map_sample_name_nci60,
                   "tcga" = map_sample_name_tcga)
map_sample_name <- function(res_df, dataset) {
  if(dataset %in% names(sample_funmap)) {
    sample_funmap[[dataset]](res_df)
  } else {
    res_df
  }
}

process_dataset <- function(dataset_suffix, data_type) {
  # Load in tool output as data frames
  tool_output_folder = str_glue("temp/output_tools_{dataset_suffix}/")
  resultsLoc <- paste0(tool_output_folder, data_type)
  
  # Some datasets contain data for only one data type
  if (!dir.exists(resultsLoc)) {
    return(tibble())
  }
  
  resultFolders <- dir(path = resultsLoc, include.dirs = T, full.names = T)
  
  cat("Reading in and converting result files to data frame:\n")
  
  map_dfr(resultFolders, parse_tool_output, conversion_functions) %>%
    # Crop to 4 digit resolution
    mutate(across(-tool, str_replace, "^([0-9]+:[0-9]+).*", "\\1")) %>%
    map_sample_name(dataset_suffix)
}

# Maps:
# Map directory name to script filename prefix if different.
mapped_tools <- c(
  "HLAminer" = "HLAminerHPRA",
  "HLAScan" = "HLAscan",
  "Kourami" = "kourami",
  "HLA-VBSeq" = "HLAVBseq"
)

conversion_functions <- load_conversion_functions(conversion_folder)

predictions_cdf <- list("dataset_suffix" = dataset_suffix,
     "data_type" = c("DNA-Seq", "RNA-Seq")) %>%
  # Create a dataframe with all (dataset, data type) combinations
  cross_df()
  
# predictions_cdf <- tibble(dataset_suffix = "nci60_1kg", data_type = c("DNA-Seq", "RNA-Seq"))

predictions_df = predictions_cdf %>%
  # The TCGA dataset is an exception
  # Here a slightly different folder structure was used
  mutate(data_type = if_else(dataset_suffix == "tcga", "gathered", data_type)) %>%
  distinct() %>%
  # Process all predictions per dataset and data type
  mutate(data = furrr::future_map2(dataset_suffix, data_type, process_dataset)) %>%
  unnest(data) %>%
  # Special case for TCGA: rename sample ID and set data_type correctly
  group_by(dataset_suffix) %>%
  group_modify(map_sample_name_tcga2) %>%
  ungroup %>%
  # Cleanup the names of the dataset
  mutate(dataset = dataset_names[dataset_suffix], dataset_suffix = NULL) %>%
  # Move all "metadata" columns to the left (unnecessary, but looks better)
  select(dataset, data_type, tool, sample_id, everything())

# Add NAs to predictions for samples that were not predicted by a tool,
# but were predicted by others
existing_tools <- predictions_df %>%
  distinct(dataset, data_type, tool)
existing_samples <- predictions_df %>%
  distinct(dataset, data_type, sample_id)
tools_all_samples <-
  left_join(existing_tools, existing_samples, by = c("dataset", "data_type"))
predictions_df <- predictions_df %>%
  full_join(tools_all_samples,
            by = c("dataset", "data_type", "tool", "sample_id"))

# Save as a nested list
predictions_list <- predictions_df %>%
  # Nest value columns
  nest_by(dataset, data_type, tool) %>%
  # Nest by the two subsequent group levels:
  # - Nested named list, indexed by tool
  # <- Named list with data type as keys
  # <-- Named list with data set as keys
  group_by(dataset, data_type) %>%
  list %>% append(1:2) %>%
  reduce(~ summarise(., data = list(deframe(cur_data())))) %>%
  # Convert final data frame to named list
  deframe

saveRDS(predictions_list, "data/hla_predictions.rds")

predictions_tools_ggroup = predictions_list %>%
  map_depth(3, allele_map)

saveRDS(predictions_tools_ggroup, "data/hla_predictions_ggroup.rds")
