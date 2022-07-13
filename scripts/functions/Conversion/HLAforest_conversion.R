
# Read in the output of hlaforest and store as R dataframe
toolOutputToR.hlaforest <- function(outputFolder){
  
  # Get a list of file names in the given outputFolder
  fileList <- list.files(path = outputFolder, pattern = "\\.txt$")
  
  # Define the needed loci
  geneTypes <- c("A","B","C","DPA1","DPB1","DQA1","DQB1","DRB1")
  
  # Create a data frame to store the results in
  results <- data.frame(matrix(NA, nrow = length(fileList), ncol = length(geneTypes)*2))
  names(results) <- rep(geneTypes, each = 2)
  IDs <- c()
  
  # For every file in the given folder, extract and store sample ID and genotype output
  for(i in 1:length(fileList)){
    # Load in file as R object
    file <- readLines(paste0(outputFolder, "/", fileList[i]), warn=FALSE)
    
    # Get sample ID
    IDs[i] <- strsplit(fileList[i], "_")[[1]][1]
    
    # Extract predicted alleles
    for (line in file) {
      resolutionLevels <- strsplit(line, "\t")[[1]]
      if(strsplit(resolutionLevels[1], ":")[[1]][2] %in% geneTypes){
        gene <- strsplit(resolutionLevels[1],":")[[1]][2]
        allele <- gsub(paste0("ROOT:",gene,":"), "", resolutionLevels[length(resolutionLevels)])
        
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
  }
  
  # Add the IDs as rownames to the data frame
  row.names(results) <- IDs
  
  return(results)
}