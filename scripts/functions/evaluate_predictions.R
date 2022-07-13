library(tidyverse)

tupelize_alleles <- function(df) {
  # Join alleles at different rows to one row = (allele1, allele2)
  # We need to handle 1 special case in case of a homozygous prediction
  
  # In case of a homozygous prediction:
  # - gold standard is also homozygous and equals prediction: 2 times correct
  # - gold standard is heterozygous and one allele matches: 1 time correct, 1 time false
  # -> The only special case
  # - none of the alleles match: 0 correct
  
  # Solution: count alleles "without replacement"
  df %>%
    group_by(sample_id, gene, allele) %>%
    mutate(allele = paste0(allele, ".", row_number())) %>%
    group_by(sample_id, gene) %>%
    summarise(list(allele), .groups="drop")
}

combine_data <- function(pred, gs) {
  predicted_tuples <- pred %>%
    # Do not consider NAs in prediction
    filter(!is.na(allele)) %>%
    # column "prediction": two predicted alleles for sample, gene pair
    tupelize_alleles() %>%
    rename(prediction = `list(allele)`)
  
  gs_tuples <- # column "truth": two gold standard alleles for sample, gene pair
    gs %>%
    # Do not consider NAs in gold standard
    filter(!is.na(allele)) %>%
    # column "prediction": two predicted alleles for sample, gene pair
    tupelize_alleles() %>%
    rename(truth = `list(allele)`)
  
  # Inner join: consider a record when the prediction is made ...
  # ... and the gold standard is available
  inner_join(predicted_tuples, gs_tuples, by = c("sample_id", "gene"))
}

eval_prediction <- function(combined_data) {
  combined_data %>%
    # Check whether predicted alleles are in the gold standard for that (sample, gene)
    mutate(is_correct = map2(prediction, truth, ~ (..1 %in% ..2))) %>%
    # "is_correct": TRUE / FALSE is given in order of predicted alleles
    select(sample_id, gene, allele = prediction, is_correct) %>%
    # Separate alleles over two rows
    unnest(c(allele, is_correct)) %>%
    # Remove counter from allele name
    mutate(allele = str_remove(allele, "\\..*$"))
}

# Calculate accuracy
calc_accuracy <- function(evaluated_data) {
  evaluated_data %>%
    # Calculate accuracy per gene
    # Number of correct predictions / number of predictions
    # For each patient we have two predictions (one per allele)
    group_by(gene) %>%
    summarise(accuracy = mean(is_correct), .groups="drop")
}

source("scripts/functions/calc_undecided.R")

# =========================================

evaluate_predictions <- function(predictions, gs) {
  # Convert gold standard to long table format
  # Actually this only needs to be done once per data set
  # Previously this was done in "process_dataset"
  # Done here to have a clear function syntax:
  # predictions and gs are both in wide format
  gs <- gs %>%
    gather("gene", "allele",-sample_id) %>%
    tidyr::separate(gene, c("gene", "idx"), sep = "\\.")
  
  # Process predictions and store in a data frame
  metrics_df <- tibble()
  evaluated_df <- tibble()
  
  # Assess performance per tool
  for (tool in names(predictions)) {
    pred_wide <- predictions[[tool]]
    
    # Convert to long table format
    pred <- gather(pred_wide, "gene", "allele", -sample_id) %>%
      tidyr::separate(gene, c("gene", "idx"), sep = "\\.")
    
    # Combine prediction and accuracy data
    combined_data <- combine_data(pred, gs)
    
    # When the gold standard is not available for both alleles
    # we cannot sensibly determine true / false
    # e.g. prediction for both alleles, only one allele found
    # What do we do with the rest
    # TODO: This is problematic for the NCI-60, HLA-C
    # Only 17 samples have a ground truth for both alleles!
    combined_data <- combined_data %>%
      filter(map_int(truth, length) == 2)
    
    # Percentage undecided
    # "Undecided": prediction = NA
    # Original idea: calculate undecided as the number of samples (out of the total)
    # for which we have an NA
    # Better calculate this metric on the sample_id, gene pairs that exist in the gold standard
    # -> GS is not NA
    undecided <- pred %>%
      semi_join(gs %>%
                  # Do not consider NAs in gold standard
                  filter(!is.na(allele)), by = c("sample_id", "gene")) %>%
      calc_undecided()
    
    evaluated <- eval_prediction(combined_data)
    evaluated_df <- bind_rows(evaluated_df, tibble(tool,
                                                   iscorrect =
                                                     list(evaluated)
    ))
    
    accuracy <- evaluated %>% calc_accuracy()
    
    # Create one table with all metrics
    # Add the rows to the summarizing table
    metrics <- reduce(
      # Which metrics to include
      list(undecided, accuracy),
      # Join by gene
      full_join,
      by = "gene"
    )
    
    metrics_df <- bind_rows(metrics_df, tibble(tool, metrics = list(metrics)))
  }
  
  # Store as nested list instead of a dataframe
  list(evaluated_df, metrics_df) %>%
    map(deframe)
}

process_dataset <- function(dataset_name, gs_suffix, df_all) {
  # Load corresponding gold standard data
  gs_file <- str_glue("data/gold_standard_{gs_suffix}.rds")
  gs <- readRDS(gs_file)
  
  # Assess performance for all data types in this dataset
  map(df_all[[dataset_name]], ~evaluate_predictions(., gs))
}