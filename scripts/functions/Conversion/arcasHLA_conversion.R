
toolOutputToR.arcasHLA <- function(outputFolder){
  
  # Load needed packages
  library("rjson")
  
  # Get a list of .json file names in the given outputFolder
  fileList <- list.files(outputFolder, pattern = "\\.genotype.json$", full.names = T)
  
  # Define the expected gene types
  geneTypes <- c("A","B","C","DPA1","DPB1","DQA1","DQB1","DRB1")
  
  # Create a data frame to store the results in
  results <- data.frame(matrix(NA, nrow = length(fileList), ncol = length(geneTypes)*2))
  names(results) <- rep(geneTypes, each = 2)
  IDs <- c()
  
  # For every file in the given folder, extract and store sample ID and genotype output
  for(i in 1:length(fileList)){
    # Load in file as R object
    alleleList <- fromJSON(file = fileList[i])
    
    # Extract sample ID
    IDs[i] <- basename(fileList[i]) %>% str_replace('\\.genotype.json$', '')
    
    # For every allele in the list, store the type of gene and allele result
    for (allele in alleleList) {
      gene <- strsplit(allele[1], "[*]")[[1]][1]
      allele1 <- strsplit(allele[1], "[*]")[[1]][2]
      if(length(allele)==2){
        allele2 <- strsplit(allele[2], "[*]")[[1]][2]
      } else {
        allele2 <- allele1
      }
      
      # Find the corresponding column index for the gene
      colIndex <- which(names(results) %in% gene)
      
      # Store allele in results
      results[i, colIndex[1]] <- allele1
      results[i, colIndex[2]] <- allele2
    }
  }
  
  # Add the IDs as rownames to the data frame
  row.names(results) <- IDs
  
  return(results)
}
