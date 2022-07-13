#######################
## PHLAT output to R ##
#######################

toolOutputToR.PHLAT <- function(outputFolder) {
  
  # Get a list of specific sum files in the given output folder
  files <- list.files(path = outputFolder, pattern = "*_HLA.sum$")
  
  # Define the expected loci
  loci <- c("A","B","C","DPA1","DPB1","DQA1","DQB1","DRB1")
  
  # Create a data frame to store the results in
  results <- data.frame(matrix(NA, nrow = length(files), ncol = length(loci) * 2))
  names(results) <- rep(loci, each = 2)
  IDs <- c()
  
  # For every file in the given folder, extract and store sample ID and genotype output
  for (i in 1:length(files)) {
    # Load in file as R object
    data <- read.delim(paste(outputFolder, files[i], sep = "/"), sep = "\t", header = TRUE)
    
    # Extract sample ID
    IDs[i] <- gsub('\\_HLA\\.sum$', '', basename(files[[i]]))
    
    # Loop through alleles, store the type of gene and allele result
    for (row in 1:nrow(data)) {
      
      # Get locus from first col
      locus <- gsub("HLA_(.+)", "\\1", data[row, 1])
      
      # Early exit loop iteration if we're not interested in the given allele
      if (nrow(data) == 0) {
        next
      }
      if (!(locus %in% loci)) {
        next
      }
      
      
      # Check if predictions are present
      # If so, remove redundant prefix and locus from allele predictions
      if (grepl("[A-Z]{1,3}[0-9]*\\*[0-9]+", data[row, "Allele1"])) {
        # Remove allele name and '*', e.g remove DRB1*
        allele_1 <- gsub("[A-Z]{1,3}[0-9]*\\*", "", data[row, "Allele1"])
        # ensure trimming of whitespace
        allele_1 <- gsub("\\s", "", allele_1)
        # if multiple alleles e.g. "01:22N,01:107" take the first one
        allele_1 <- gsub(",.*", "", allele_1)
      } else {
        allele_1 <- NA
      }
      
      if (grepl("[A-Z]{1,3}[0-9]*\\*[0-9]+", data[row, "Allele2"])) {
        # Remove allele name and '*', e.g remove DRB1*
        allele_2 <- gsub("[A-Z]{1,3}[0-9]*\\*", "", data[row, "Allele2"])
        # ensure trimming of whitespace
        allele_2 <- gsub("\\s", "", allele_2)
        # if multiple alleles e.g. "01:22N,01:107" take the first one
        allele_2 <- gsub(",.*", "", allele_2)
      } else {
        allele_2 <- NA
      }
      
      colIndex <- which(names(results) %in% locus)
      
      # Store results
      results[i, colIndex[1]] <- allele_1
      results[i, colIndex[2]] <- allele_2
      
    }
    
  }
  
  # Add the IDs as rownames to the data frame
  row.names(results) <- IDs
  
  return(results)
}
