
# Read in the output of hlaforest and store as R dataframe
toolOutputToR.HLA_HD <- function(outputFolder){
  
  # Get a list of file names in the given outputFolder
  fileList <- list.files(path = outputFolder, pattern = "final\\.result\\.txt$")
  
  # Define the needed loci
  geneTypes <- c("A","B","C","DPA1","DPB1","DQA1","DQB1","DRB1")
  
  # Create a data frame to store the results in
  results <- data.frame(matrix(NA, nrow = length(fileList), ncol = length(geneTypes)*2))
  names(results) <- rep(geneTypes, each = 2)
  IDs <- c()
  
  # For every file in the given folder, extract and store sample ID and genotype output
  for(i in 1:length(fileList)){
    # Load in file as R object
    file <- readLines(paste0(outputFolder,"/",fileList[i]))
    
    # Get sample ID
    IDs[i] <- str_remove(basename(fileList[i]), '_final\\.result\\.txt$')
    
    # Extract predicted alleles
    for (line in file) {
      typing <- strsplit(line, "\t")[[1]]
      if(typing[1] %in% geneTypes){
        gene <- typing[1]
        alleles <- c(typing[2],typing[3])
        
        # Find the corresponding column index for the gene
        colIndex <- which(names(results) %in% gene)
        
        if(alleles[2] == "-"){
          alleles[2] <- alleles[1]
        }
        
        # Store allele in result matrix
        for (j in 1:2) {
          if(alleles[j] != "Not typed"){
            results[i, colIndex[j]] <- strsplit(alleles[j],"[*]")[[1]][2]
          }
        }
      }
    }
  }
  
  # Add the IDs as rownames to the data frame
  row.names(results) <- IDs
  
  return(results)
}