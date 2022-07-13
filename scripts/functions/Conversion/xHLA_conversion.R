
# Read in the output of a xHLA and store as R dataframe

toolOutputToR.xHLA <- function(outputFolder){
  
  # Load needed packages
  library("rjson")
  
  # Get a list of .json file names in the given outputFolder
  fileList <- list.files(outputFolder, pattern = "\\.json$", full.names = T)
  
  # Define the expected gene types
  geneTypes <- c("A","B","C","DPA1","DPB1","DQA1","DQB1","DRB1")
  
  # Create a data frame to store the results in
  results <- data.frame(matrix(NA, nrow = length(fileList), ncol = length(geneTypes)*2))
  names(results) <- rep(geneTypes, each = 2)
  IDs <- c()
  
  # For every file in the given folder, extract and store sample ID and genotype output
  for(i in 1:length(fileList)){
    # Load in file as R object
    file <- fromJSON(file = fileList[i])
    
    # Extract sample ID and a list of alleles
    IDs[i] <- file$sample_id
    alleleList <- file$hla$alleles
    
    # For every allele in the list, store the type of gene and allele result
    for (allele in alleleList) {
      gene <- strsplit(allele, "[*]")[[1]][1]
      allele <- strsplit(allele, "[*]")[[1]][2]
      
      # Find the corresponding column index for the gene
      colIndex <- which(names(results) %in% gene)
      
      # Store allele in the first of the two gene columns when available, otherwise in the second
      if(is.na(results[i, colIndex[1]])){
        results[i, colIndex[1]] <- allele
      } else if(is.na(results[i, colIndex[2]])){
        results[i, colIndex[2]] <- allele
      } else {
        print("PROBLEM: no empty spots available to store more genotypes")
      }
    }
  }
  
  # Add the IDs as rownames to the data frame
  row.names(results) <- IDs
  
  return(results)
}
