library(tidyverse)
future::plan(future::multicore)
options(future.rng.onMisuse="ignore")

# Load functions
source("scripts/functions/extract_dataset.R")
source("scripts/functions/build_group_arg_list.R")

source("scripts/functions/meta_prediction.R")
source("scripts/functions/evaluate_predictions.R")
source("scripts/functions/exact_genotype_check.R")

# Read gold standard data
gs <- readRDS("data/gold_standard_1kg.rds")

# Read predictions for 1000 genomes data
# (on which we optimize the metaclassifier)
predictions_all <- readRDS("data/hla_predictions_ggroup.rds")
pred_ds <- extract_dataset("1kg", predictions_all)

# Convert predictions to long format
pred_long <- pred_ds %>%
  gather("locus", "prediction", -data_type, -tool, -sample_id) %>%
  separate(locus, c("gene", "idx")) %>%
  # Add a column for the MHC class
  mutate(mhc_class = if_else(gene %in% c("A", "B", "C"), "I", "II"))

# List all genes for which we have a gold standard available
# For other genes we cannot train the metaclassifier
supported_genes <- gs %>%
  colnames() %>%
  setdiff("sample_id") %>%
  str_split('\\.') %>%
  map(1) %>%
  unique

# The first two tools to select
initial_tools <- readRDS("data/metaclassifier_dependencies/initial_tools.rds")

# Scores per tool to prefer the output of a particular tool in case of ties
tie_solver <- readRDS("data/metaclassifier_dependencies/tie_solver.rds")

evaluate_selection <- function(predictions, key, tie_solver, tools) {
  print(tools)
  current_meta_pred <- predictions %>%
    # OPTION: choose "gene" or "mhc_class" here to select a different tie solver
    meta_prediction(tie_solver, "gene", tools) %>%
    # Convert to required format for evaluation
    ungroup %>%
    gather("idx", "allele", `1`, `2`) %>%
    mutate(locus=str_glue("{gene}.{idx}")) %>%
    # Remove index and gene to prevent NA while converting to wide format
    # Also remove data_type:
    # <- would cause issues with the evaluation function
    select(-idx, -gene, -data_type) %>%
    spread(locus, allele)
  
  if (nrow(current_meta_pred) != 0) {
    # Calculate the per-allele metrics
    # (Check whether allele is correct)
    e_allele <- evaluate_predictions(list(typing=current_meta_pred), gs) %>%
      chuck(2) %>%
      chuck("typing") %>%
      transmute(gene, undecided_ratio, accuracy, S=(1 - undecided_ratio) * accuracy)
    
    # Calculate the per-genotype metrics
    # (Call is correct if both alleles are correct)
    e_genotype <- exact_genotype_check(current_meta_pred, gs)
    e <- left_join(e_allele, e_genotype, by="gene")
  } else {
    # No valid prediction made -> score is 0 for all genes
    e <- tibble(gene = unique(predictions$gene), undecided_ratio=1, accuracy=0, S=0, H=0)
  }
  e
}

# Helper function: combine selected tools one by one with all other tools
combine_with_others <- function(selected, tools) {
  others <- setdiff(tools, selected)
  map(others, ~union(selected, .))
}

optimize_metaclassifier <- function(predictions, group_keys, initial_tools, tie_solver) {
  print(initial_tools)
  print(group_keys)
  
  # Select the initial tools for the current group
  # (e.g. data type and MHC class)
  selected_tools <- initial_tools %>%
    semi_join(group_keys, by = colnames(group_keys)) %>%
    pull("tools")
  
  # Output dataframe: will be filled in the loop below
  results <- tibble()
  
  # Add additional selected tools one by one until all tools are added
  n = length(unique(predictions$tool))
  for (k in 2:n) {
    # Make and evaluate predictions for all selected tools
    results_new <- selected_tools %>%
      # Make and evaluate the meta-predictions
      # Create data frame with metrics per tuple of tools
      # In each map call (per "tools" tuple) we need to create 1 row data frames
      # Arguments of tibble are columns: therefore nest everything in 1 element lists
      furrr::future_map_dfr( ~ tibble(tools = list(.),
                        val=list(evaluate_selection(predictions, group_keys, tie_solver, .)))) %>%
      # Select the best tool combination until now
      # (take the tool that improves the accuracy the most)
      slice_max(map_dbl(val, ~mean(.$S, na.rm=T))) %>%
      # Ensure there are no ties
      {
        predictions <- .
        if(nrow(predictions) != 1) {
          cat("Tie at level", k, "with:", paste0(predictions$tools, collapse=", "), "\n")
          
          current_tie_solver <- tie_solver %>%
            # Select tie solver for the correct group (e.g. data type and MHC class)
            semi_join(group_keys, by=colnames(group_keys)) %>%
            # For each tool, calculate the mean score
            # for all genes under consideration
            # ignore unsupported genes, but if all are unsupported
            # then give a score of zero to the tool
            
            # Why group_by(across(...))?
            # -> Also add the group keys to the group
            # -> so these columns are still available after summarise
            group_by(across(all_of(colnames(group_keys))), tool) %>%
            summarise(S_mean=coalesce(mean(S, na.rm=T), 0), .groups = "drop")
          
          # Let decision depend on mean score (S_mean) of last added tool
          predictions <- predictions %>%
            # Extract last tool
            mutate(last_tool = map_chr(tools, ~.[[length(.)]])) %>%
            left_join(current_tie_solver,by=c("last_tool" = "tool")) %>%
            slice_max(S_mean, n=1) %>%
            select(tools, val) %>%
            ungroup
        }
        
        # If still not solved choose random selection to continue:
        slice_sample(predictions, n=1)
      } %>%
      # Extract metrics for the current best tool combination
      unnest(val)
    
    results <- bind_rows(results, results_new)
    
    # Combine current best tool combination with all others
    selected_tools <- results_new %>%
      pull("tools") %>%
      pluck(1) %>%
      combine_with_others(unique(predictions$tool))
  }
  # Add group keys to results
  # Remove possible duplicate columns between group_keys and results
  # (e.g. genes for gene-specific optimization)
  # (Analysis preformed for only 1 group_key)
  group_keys %>%
    select(-any_of(colnames(results))) %>% 
    bind_cols(results)
}

# Prepare list of arguments for call_meta_compact function
meta_p_args <- pred_long %>%
  # For optimizing the metaclasifier, we can only consider the supported genes
  # (no gold standard available for DPA1, DPB1)
  filter(gene %in% supported_genes) %>%
  # Drop unsupported gene, tool pairs
  drop_na 

# Optimized per MHC class
# -----------------------
opt_class <- meta_p_args %>%
  group_by(data_type, mhc_class) %>%
  build_group_arg_list("predictions", "group_keys") %>%
  furrr::future_pmap_dfr(optimize_metaclassifier, initial_tools$others$mhc_class, tie_solver)

# Special version for the TCGA data
# HLA-HD was not applied on the TCGA data
# Check what is the best RNA metaclassifier excluding HLA-HD
# -----------------------
opt_class_tcga_rna <- meta_p_args %>%
  # Only re-evaluate the RNA-Seq classifier
  filter(data_type == "RNA-Seq") %>%
  # Exclude HLA-HD for this analysis
  filter(tool != "HLA-HD") %>%
  # Same as before: perform per class optimization
  group_by(data_type, mhc_class) %>%
  build_group_arg_list("predictions", "group_keys") %>%
  # Use the initial tools for TCGA (excludes HLA-HD)
  furrr::future_pmap_dfr(optimize_metaclassifier, initial_tools$tcga$mhc_class, tie_solver)

# Metaclassifier for DNA remains the same
opt_class_tcga <- opt_class %>%
  filter(data_type == "DNA-Seq") %>%
  # Add the TCGA RNA classifier to the data frame
  bind_rows(opt_class_tcga_rna)

# Optimized per gene
# -----------------------
opt_gene <- meta_p_args %>%
  group_by(data_type, gene) %>%
  build_group_arg_list("predictions", "group_keys") %>%
  furrr::future_pmap_dfr(optimize_metaclassifier, initial_tools$others$gene, tie_solver)

# Save the results
list("class" = opt_class,
     "class_tcga" = opt_class_tcga,
     "gene" = opt_gene) %>%
  saveRDS("data/meta_optimized.rds")