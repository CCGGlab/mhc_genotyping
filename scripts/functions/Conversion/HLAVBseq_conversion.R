
toolOutputToR.HLA_VBSeq <- function(outputFolder){
  
  # Get a list of .json file names in the given outputFolder
  fileList <- list.files(outputFolder, pattern = "d4.txt$", full.names = T)
  
  # Define the expected loci
  geneTypes <- c("A","B","C","DPA1","DPB1","DQA1","DQB1","DRB1")
  
  # Create a data frame to store the results in
  results <- data.frame(matrix(NA, nrow = length(fileList), ncol = length(geneTypes)*2))
  names(results) <- rep(geneTypes, each = 2)
  IDs <- c()
  
  # For every file in the given folder, extract and store sample ID and genotype output
  for(i in 1:length(fileList)){
    
    # Get sample ID
    IDs[i] <- strsplit(basename(fileList[i]), ".report")[[1]][1]
    
    if(file.info(fileList[i])$size != 0){ # If file is empty, leave all NAs in results as they are
      # Load in file as R object
      file <- read.table(fileList[i], sep = "\t", header = T)
      
      # For every allele in the list, store the type of gene and allele result
      # When multiple alleles are predicted, take the first one (!!!)
      all = c()
      for (geneType in geneTypes) {
        
        if(geneType %in% file$Gene){
          # extracting alleles for each genetype
          alleles <- file[which(file$Gene == geneType),]
          # Adapted from student version to prevent the code from crashing on NAs
          allele_1 <- str_remove(alleles$Allele1, "^.*\\*")
          allele_2 <- str_remove(alleles$Allele2, "^.*\\*")
          
          allele_1 = if_else(allele_1 == "", NA_character_, allele_1)
          allele_2 = if_else(allele_2 == "", NA_character_, allele_2)
          
          all <- c(all, allele_1, allele_2)
        } else {
          all <- c(all, NA, NA)
        }
        
      }
      # Add to result data frame
      results[i,] <- all
    }
  }
  
  # Add the IDs as rownames to the data frame
  row.names(results) <- IDs
  
  return(results)
}
