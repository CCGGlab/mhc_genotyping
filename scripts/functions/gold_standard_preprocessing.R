# Author: Jasper Staut

prepareGoldStandard <- function(goldDataLocation, digitResolution){
  
  # Read in the gold standard data set
  data <- read.table(goldDataLocation, header = T, sep = "\t", quote = "")
  
  # Extract typing information of the relevant alleles
  alleles.data <- data.frame(data$HLA.A.1,data$HLA.A.2,data$HLA.B.1,data$HLA.B.2,data$HLA.C.1,
                             data$HLA.C.2,data$HLA.DQA1.1,data$HLA.DQA1.2,data$HLA.DRB1.1,
                             data$HLA.DRB1.2,data$HLA.DQB1.1,data$HLA.DQB1.2)
  # Explicitly indicate the allele names extracted above (in the right order)
  extractedLoci <- rep(c("A","B","C","DQA1","DRB1","DQB1"), each = 2)
  
  # Store sample names (and clean them up, for example NA18579_1 -> remove underscore)
  sampleNames <- c()
  for (i in 1:nrow(data)) {sampleNames[i] <- strsplit(data$Subject[i], "_")[[1]][1]}
  
  # Clean up the alleles and store in nested lists format, with the different ambiguous allele
  # possibilities stored as a vector
  alleles <- list()
  for(locus in 1:ncol(alleles.data)){
    
    # Split ambiguous allele possibilities into vector with the alleles
    alleles[[locus]] <- strsplit(alleles.data[,locus],"/")
    
    for(sample in 1:nrow(alleles.data)) {
      
      # Fill in NA if only an empty character vector is present for an allele
      if(length(alleles[[locus]][[sample]]) == 0){
        alleles[[locus]][[sample]] <- NA
        
        # Otherwise, go over all possible ambiguous allele options and clean up allele types
      } else {
        for(i in 1:length(alleles[[locus]][[sample]])){
          # Remove first part A* in A*01:01:01 for example
          alleles[[locus]][[sample]][i] <- gsub("[a-zA-Z0-9]*[*]","",alleles[[locus]][[sample]][i])
          # Remove letters, for example in 02:01N
          alleles[[locus]][[sample]][i] <- gsub("[A-Z]","",alleles[[locus]][[sample]][i])
        }
      }
    }
  }
  
  # Adapt the gold standard data to a given resolution level
  for(locus in 1:length(alleles)){
    for(sample in 1:length(alleles[[1]])) {
      if(!is.na(alleles[[locus]][[sample]][1])){
        
        # Vector to store indices (of ambiguous possibilities vector) with insufficient resolution
        lowResolution <- c()
        
        for(i in 1:length(alleles[[locus]][[sample]])){

          # split a single allele into a vector with the digits of each resolution level 
          digitsVector <- strsplit(alleles[[locus]][[sample]][i], ":")[[1]]
          
          # If the resolution of a given allele is insufficient, store its index in lowResolution
          if(length(digitsVector) < digitResolution/2){
            lowResolution[length(lowResolution)+1] <- i
            
          # If resolution is higher than required, retain only the needed digits and join 
          # these into one string again
          } else if(length(digitsVector) > digitResolution/2){
            alleles[[locus]][[sample]][i] <- paste(digitsVector[1:(digitResolution/2)], collapse = ":")
          }
          
        }
        
        # Remove alleles with insufficient resolution from the vector
        if(length(lowResolution) > 0){
          alleles[[locus]][[sample]] <- alleles[[locus]][[sample]][-lowResolution]
        }
        
        # Check if vector of allele possibilities is empty. If so, set to NA
        if(length(alleles[[locus]][[sample]]) == 0){
          alleles[[locus]][[sample]][1] <- NA
          
          # Otherwise, deduplicate the elements in the ambiguous allele possibilities vector
        } else {
          alleles[[locus]][[sample]] <- unique(alleles[[locus]][[sample]])
        }
      }
    }
  }
  return(list(extractedLoci, sampleNames, alleles))
}
