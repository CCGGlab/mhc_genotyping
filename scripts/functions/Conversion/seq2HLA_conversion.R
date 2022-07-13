##########################
## seq2HLA  output to R ##
##########################
seq2HLA.getSingleAllele <- function(data_string) {
  
  if (!is.character(data_string)) {
    stop(sprintf("Unexpected data provided as input argument: %s. Expected character.",
                 typeof(data_string))) 
  }
  
  if (grepl("[A-Z]{1,3}[0-9]*\\*[0-9]+", data_string)) {
    # Remove allele name and '*', e.g remove DRB1*
    allele <- gsub("[A-Z]{1,3}[0-9]*\\*", "", data_string)
  } else {
    allele <- NA
  }
  
  return(allele)
}

seq2HLA.writeHLAClassResult <- function(hla_class, outputFolder, loci, results) {
  # note: R is weakly typed, thus 1 == "1"
  valid_class_aliases <- c(1, "I", 2, "II")
  
  if (!is.element(hla_class, valid_class_aliases)) {
    stop(sprintf("Unexpected HLA class provided as input argument: %s. Expected one of %s", 
                 hla_class, 
                 paste(valid_class_aliases, sep = ", ", collapse = TRUE)))
  }
  
  # get the character representation of the class
  class_str <- switch(hla_class, "1" = "I", "2" = "II", hla_class)
  
  pattern <- switch(
    class_str, 
    "I" = "ClassI-class.HLAgenotype4digits$",
    "II" = "ClassII.HLAgenotype4digits$")
  
  files <- list.files(path = outputFolder, pattern = pattern)
  
  for (file in files) {
    # Generalized for general sample name
    data <- read.delim(paste(outputFolder, file, sep = "/"), sep = "\t", header = TRUE)
    sample_name <- str_remove(file, '(-ClassI|-ClassII)?(-class|-nonclass)?\\.HLAgenotype4digits$')
    
    for (row in 1:nrow(data)) {
      # get locus from first col
      locus <- data[row, 1]
      
      # Early exit loop iteration if we're not interested in the given allele
      if (!locus %in% loci) {
        next
      }
      
      colIndex <- which(names(results) %in% locus)
      
      allele_1 <- gsub("'", "", seq2HLA.getSingleAllele(data[row, "Allele.1"]))
      allele_2 <- gsub("'", "", seq2HLA.getSingleAllele(data[row, "Allele.2"]))
      
      # ensure trimming of whitespace
      allele_1 <- gsub("\\s", "", allele_1)
      allele_2 <- gsub("\\s", "", allele_2)
      
      # Store result
      results[sample_name, colIndex[1]] <- allele_1
      results[sample_name, colIndex[2]] <- allele_2
    }
  }
  
  return(results)
}

toolOutputToR.seq2HLA <- function(outputFolder) {
  
  # define the expected loci
  loci <- c("A","B","C","DPA1","DPB1","DQA1","DQB1","DRB1")

  # Adapted: Generalized version of code students.
  sample_names <-
    list.files(path = outputFolder, pattern = "HLAgenotype4digits$") %>%
    map(~str_remove(., '(-ClassI|-ClassII)?(-class|-nonclass)?\\.HLAgenotype4digits$')) %>%
    unique
  
  # create a data frame to store the results in
  results <- data.frame(matrix(NA, nrow = length(sample_names), ncol = length(loci) * 2))
  colnames(results) <- rep(loci, each = 2)
  rownames(results) <- sample_names
  
  # store results
  results <- seq2HLA.writeHLAClassResult("I", outputFolder, loci, results)
  results <- seq2HLA.writeHLAClassResult("II", outputFolder, loci, results)
  return(results)
}
