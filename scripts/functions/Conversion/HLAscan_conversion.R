
toolOutputToR.HLAscan <- function(outputFolder){
  
  # Get a list of .json file names in the given outputFolder
  fileList <- list.files(outputFolder, pattern = "\\_result.txt$", full.names = T)
  
  # Define the expected loci
  geneTypes <- c("A","B","C","DPA1","DPB1","DQA1","DQB1","DRB1")
  
  # Create a data frame to store the results in
  # Not necessary to define the number of rows???
  # length(unique(substr(fileList,1,7)))
  results <- data.frame(matrix(NA,
                               nrow = 0,
                               ncol = length(geneTypes)*2)
                        )
  names(results) <- rep(geneTypes, each = 2)
  IDs <- c()
  
  # For every file in the given folder, extract and store sample ID and genotype output
  for(i in 1:length(fileList)){
    # Load in file as R object
    lines <- readLines(fileList[i])
    
    # Extract sample ID and gene
    ID <- str_remove(basename(fileList[i]), '_[^_]*_result.txt$')
    # ID <- strsplit(basename(fileList[i]), "_")[[1]][1]
    
    if (! ID %in% IDs) {IDs[length(IDs)+1] <- ID}
    gene <- str_match(basename(fileList[i]), '_([^_]*)_result.txt$')[,2]
    
    if ((length(lines) == 0) | (substr(lines[1],1,5)=="ERROR")){
      allele1 <- NA
      allele2 <- NA
    } else if (substr(lines[length(lines)],1,14)=="HLAscan cannot"){
      allele1 <- NA
      allele2 <- NA
    } else {
      allele1 <- strsplit(lines[length(lines)-1], "\t")[[1]][2]
      allele2 <- strsplit(lines[length(lines)], "\t")[[1]][2]
    }
    
    # Find the corresponding column index for the gene
    colIndex <- which(names(results) %in% gene)
    
    # Store allele in results
    results[length(IDs), colIndex[1]] <- allele1
    results[length(IDs), colIndex[2]] <- allele2
  }
  
  # Add the IDs as rownames to the data frame
  row.names(results) <- IDs
  
  fileInfo <- strsplit(fileList, "_")
  DQA1_IDs <- c()
  for (i in 1:length(fileInfo)) {
    if(fileInfo[[i]][2]=="DQA1"){
      DQA1_IDs <- c(DQA1_IDs, fileInfo[[i]][1])
    }
  }
  DQA1_absent <- IDs[-which(IDs %in% DQA1_IDs)]
  results[which(IDs %in% DQA1_absent), which(names(results) == "DQA1")] <- "00:00"
  
  return(results)
}
