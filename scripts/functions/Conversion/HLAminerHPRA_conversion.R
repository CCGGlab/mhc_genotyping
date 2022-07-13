##########################
## HLAminer output to R ##
##########################

toolOutputToR.HLAminer_HPRA <- function(outputFolder) {
  type <- "HPRA"
  types <- c("HPRA", "HPTASR")
  
  if (!type %in% types) {
    stop(paste("Type must be one of", paste(types, collapse = ", ")))
  }
  
  # Get a list of csv file names in the given output folder
  fileList <- list.files(path = outputFolder, pattern = sprintf("^HLAminer_%s_.*\\.csv$", type))
  
  # Define the expected loci
  loci <- c("A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1")
  # Create a data frame to store the results in
  results <- data.frame(matrix(NA, nrow = length(fileList), ncol = length(loci) * 2))
  names(results) <- rep(loci, each = 2)
  IDs <- c()
  
  # For every file in the given folder, extract and store sample ID and genotype output
  for (i in 1:length(fileList)) {
    # Extract sample ID
    # Original code: gsub(".csv", "", strsplit(fileList[i], "_")[[1]][3])
    IDs[i] <- str_match(basename(fileList[i]), 'HLAminer_HPRA_(.*).csv')[,2]
    
    # Parse file line by line because it has no meaningful structure
    fileHandle <- file(paste(outputFolder, fileList[i], sep = "/"), "r")
    
    while (TRUE) {
      line <- readLines(fileHandle, n = 1)
      
      if (length(line) == 0) {
        break
      }
      
      # search and extract allele predictions from lines
      if (grepl("Prediction #[0-9]+", line)) {
        # 1st group: e.g. "1" or "2"
        # 2nd group: e.g. "DRB1" or "C"
        re <- "\\s+Prediction #([0-9]+) - ([A-Z]{1,3}[0-9]*).*"
        idx <- as.integer(gsub(re, "\\1", line))
        loc <- gsub(re, "\\2", line)
        
        # early exit if we're not interested in the locus
        if (!loc %in% loci) {
          next
        }
        
        # Find the corresponding column index for the gene
        colIndex <- which(names(results) %in% loc)
        
        # there can be multiple predictions for each position, only consider first one
        if (!is.na(results[i, colIndex[idx]])) {
          next
        }
        
        # the prediction line only contains single digit if there are ambiguous predictions
        # so always take full result from next line,
        # being the prediction with the highest score
        nextLine <- readLines(fileHandle, n = 1)
        # expect at least 4 digit resolution
        if (grepl("[A-Z]{1,3}[0-9]*\\*[0-9]+:[0-9]+", nextLine)) {
          # get line parts; order: Allele,Score,E-value,Confidence
          lineParts <- strsplit(nextLine, ",")
          # 1st group: e.g. "DRB1" or "C"
          # 2nd group: e.g. "15:01P" or "01:01:01:18"
          allele_re <- "([A-Z]{1,3}[0-9]*)?\\*([0-9]+[:0-9A-Z]+[:0-9A-Z]*[:0-9A-Z]*)"
          prediction <- gsub(allele_re, "\\2", lineParts[[1]][1])
          prediction <- gsub("\\s", "", prediction)
          
          # Store results
          results[i, colIndex[idx]] <- prediction
        }
      }
    }
    close(fileHandle)
  }
  
  # Add the IDs as rownames to the data frame
  row.names(results) <- IDs
  
  # Handle homozygous call correctly:
  # For homozygous calls the first allele is not NA, while the second is NA
  orig_colnames <- colnames(results)
  
  # Make column names unique
  pos_colnames <-
    paste(orig_colnames, rep.int(c(1, 2), length(orig_colnames) / 2), sep =
            ".")
  
  results <- results %>%
    set_names(pos_colnames) %>%
    rownames_to_column("sample_id") %>%
    # Two alleles in different columns, gene in different row
    pivot_longer(-sample_id, names_to=c("gene", "idx"), names_sep = "\\.", values_to = "allele") %>%
    pivot_wider(c(sample_id, gene), names_from="idx", names_prefix="v", values_from = "allele") %>%
    # Replace allele in v2 to v1 if it is NA
    # (Make the call homozygous)
    mutate(v2 = coalesce(v2, v1)) %>%
    # Convert back to wide format: each allele in different column
    pivot_longer(-c(sample_id, gene), names_to="idx", names_prefix="v", values_to = "allele") %>%
    pivot_wider(sample_id, names_from=c("gene", "idx"), names_sep = ".", values_from = "allele") %>%
    # Convert back to original DF format
    column_to_rownames("sample_id") %>%
    set_names(str_remove(colnames(.), '\\.[12]$'))
  
  return(results)
}
