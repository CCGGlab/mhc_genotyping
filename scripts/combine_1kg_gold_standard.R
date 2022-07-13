library(tidyverse)
# Function to rename gene columns
column_rename <- function(x) {
  str_remove(x, '\\..*') %>% as_tibble_col("gene") %>% group_by(gene) %>% mutate(idx = 1:2) %>% ungroup %>% transmute(paste0('HLA.', gene, '.', idx)) %>% deframe
}

# Standardize the output data format
uniform_columns <- c("Col", "Subject", "Sample", "Population", "Pedigree", "HLA.A.1", 
                     "HLA.A.2", "HLA.B.1", "HLA.B.2", "HLA.C.1", "HLA.C.2", "HLA.DRB1.1", 
                     "HLA.DRB1.2", "HLA.DRB3.1", "HLA.DRB3.2", "HLA.DRB4.1", "HLA.DRB4.2", 
                     "HLA.DRB5.1", "HLA.DRB5.2", "HLA.DQA1.1", "HLA.DQA1.2", "HLA.DQB1.1", 
                     "HLA.DQB1.2", "HLA.DPA1.1", "HLA.DPA1.2", "HLA.DPB1.1", "HLA.DPB1.2", 
                     "Lab", "STRMICA.1", "STRMICA.2", "D6S276.1", "D6S276.2", "D6S2239.1", 
                     "D6S2239.2", "D6S2223.1", "D6S2223.2", "D6S105.1", "D6S105.2", 
                     "D6S2222.1", "D6S2222.2", "D6S2972.1", "D6S2972.2", "D6S510.1", 
                     "D6S510.2", "D6S265.1", "D6S265.2", "D6S2685.1", "D6S2685.2", 
                     "D6S2674.1", "D6S2674.2", "D6S2673.1", "D6S2673.2", "TNFa.1", 
                     "TNFa.2", "D6S2791.1", "D6S2791.2", "D6S2789.1", "D6S2789.2", 
                     "D6S2671.1", "D6S2671.2", "D6S273.1", "D6S273.2", "D6S1666.1", 
                     "D6S1666.2", "D6S2447.1", "D6S2447.2", "D6S2663.1", "D6S2663.2", 
                     "DQIV.1", "DQIV.2", "D6S2874.1", "D6S2874.2", "D6S2656.1", "D6S2656.2", 
                     "D6S291.1", "D6S291.2")

# For 10 samples 2 rows are present
# one with higher and one with lower resolution
# or unknown vs known
# prefer the rows with high resolution over low resolution
# prefer known over unknown
preferred_rows <-
  c(
    "NA18577_2",
    "NA18579_1",
    "NA18612_1",
    "NA18949_2",
    "NA18968_1",
    "NA18997_1",
    "NA19119_1",
    "NA19210_2",
    "NA19221_1",
    "NA19223_2"
  )

# Load the data
goldStandardLocation <- "downloads/1kg/gold standard/"

gourraud_df <- read_delim(paste0(goldStandardLocation, "20140702_hla_diversity.txt"), delim = ' ') %>%
  group_by(id) %>%
  mutate(num_present = length(id)) %>%
  mutate(id_idx = paste0(id, '_', seq_along(id))) %>%
  ungroup %>%
  filter((num_present == 1) | (id_idx %in% preferred_rows)) %>%
  select(-id_idx, -num_present) %>%
  rename(Subject = id) %>%
  select(-sbgroup) %>%
  rename_with(column_rename, .cols = -Subject)

missing_cols <- setdiff(uniform_columns, colnames(gourraud_df))
empty_cols <- as_tibble_col(missing_cols) %>% mutate(v = NA) %>% deframe

gourraud_df <- gourraud_df %>%
  add_column(!!!empty_cols) %>%
  mutate(Col = row_number()) %>%
  select(!!uniform_columns)

# Read the two supplementary tables of the De Bakker data
## Supplementary tables 3 and 7 of seq2HLA paper (doi:10.1186/gm403)
debakker1 = read.csv2(paste0(goldStandardLocation,"/DeBakker1.csv"))
# Table I of doi:10.5281/zenodo.1338172
# https://doi.org/10.5281/zenodo.1338172
debakker2 = read.table(paste0(goldStandardLocation,"/DeBakker2.txt"))

# Standardize debakker1 format
## Column names
names(debakker1) = c("Sample",
                     "A",
                     "A",
                     "B",
                     "B",
                     "C",
                     "C",
                     "DRB1",
                     "DRB1",
                     "DQB1",
                     "DQB1",
                     "DQA1",
                     "DQA1")
for (i in 1:nrow(debakker1)) {
  for(j in 2:ncol(debakker1)){
    # "-" -> no call present: replace by empty string
    if(debakker1[i,j] == "-"){debakker1[i,j] = ""}
    # - Only keep allele name (split on star, take second part)
    # - Remove :xx from allele.
    #   (This means we don't have this sample typed at the requested resolution)
    else {debakker1[i,j] = gsub(":xx","",strsplit(debakker1[i,j],"[*]")[[1]][2])}
  }
}

# Standardize debakker2 format
## Remove column 2 and 3: we are not interested in the SRR IDs
debakker2 = debakker2[,-c(2,3)]
## Column names
names(debakker2) = c("Sample",
                     "A",
                     "A",
                     "B",
                     "B",
                     "C",
                     "C",
                     "DRB1",
                     "DRB1",
                     "DQB1",
                     "DQB1",
                     "DQA1",
                     "DQA1")
for (i in 1:nrow(debakker2)) {
  for(j in 2:ncol(debakker2)){
    typing = toString(debakker2[i,j])
    # "-" -> no call present: replace by empty string
    if(typing == "-"){debakker2[i,j] = ""}
    # Convert allele name from old to new format
    # (insert colon between each field)
    else if(nchar(typing) == 4){debakker2[i,j] = paste0(substr(debakker2[i,j],1,2),":",substr(debakker2[i,j],3,4))}
    # Prepend zeros if necessary
    else if(nchar(typing) == 3){debakker2[i,j] = paste0("0",substr(debakker2[i,j],1,1),":",substr(debakker2[i,j],2,3))}
    else {print(paste0("problem"))} # => does not occur
    # - Remove :xx from allele.
    #   (:xx means we don't have this sample typed at the requested resolution)
    debakker2[i,j] = gsub(":xx","",debakker2[i,j])
  }
}

# Check if overlapping samples of both De Bakker subsets have the same typing result
for (i1 in 1:nrow(debakker1)) {
  for (i2 in 1:nrow(debakker2)) {
    if(debakker1$Sample[i1] == debakker2$Sample[i2]){
      print(debakker1$Sample[i1])
      print(sum(c(debakker1[i1,2]!=debakker2[i2,2],debakker1[i1,3]!=debakker2[i2,3])))
      print(sum(c(debakker1[i1,4]!=debakker2[i2,4],debakker1[i1,5]!=debakker2[i2,5])))
      print(sum(c(debakker1[i1,6]!=debakker2[i2,6],debakker1[i1,7]!=debakker2[i2,7])))
      print(sum(c(debakker1[i1,8]!=debakker2[i2,8],debakker1[i1,9]!=debakker2[i2,9])))
      print(sum(c(debakker1[i1,10]!=debakker2[i2,10],debakker1[i1,11]!=debakker2[i2,11])))
      print(sum(c(debakker1[i1,12]!=debakker2[i2,12],debakker1[i1,13]!=debakker2[i2,13])))
    }
  }
}

# NA12156: HLA-C (column 6 or 7) does not correspond
# -> according to De Bakker 1: HLA-C = 15:02, 06:06
# -> according to De Bakker 2 + Gourraud: HLA-C = 15:02, 06:02

# NA12873: HLA-B (column 4 or 5) does not correspond
# -> according to De Bakker 1: HLA-B = 39:06, 07:02
# -> according to De Bakker 2: HLA-B = 39:05, 07:02
# --> Not in Gourraud dataset

# Fill in all the DQA1 info from De Bakker in the gold data
for (i in 1:nrow(gourraud_df)){
  if(gourraud_df$Subject[i] %in% debakker1$Sample){
    gourraud_df$HLA.DQA1.1[i] = debakker1[which(debakker1$Sample == gourraud_df$Subject[i]),12]
    gourraud_df$HLA.DQA1.2[i] = debakker1[which(debakker1$Sample == gourraud_df$Subject[i]),13]
  } else if(gourraud_df$Subject[i] %in% debakker2$Sample){
    gourraud_df$HLA.DQA1.1[i] = debakker2[which(debakker2$Sample == gourraud_df$Subject[i]),12]
    gourraud_df$HLA.DQA1.2[i] = debakker2[which(debakker2$Sample == gourraud_df$Subject[i]),13]
  }
}

# Append additional samples to the gold standard
allDebakkerSamples = unique(c(debakker1$Sample, debakker2$Sample))
# Get sample IDs
# -> NA18968_1 and NA18997_1 -> remove underscore, hence the regex search
goldSamples = gsub("_[0-9]$","",gourraud_df$Subject)

# Samples available in De Bakker data that were not available in
# the Gourraud data
extraSamples = allDebakkerSamples[which(!(allDebakkerSamples%in%goldSamples))]
newRow = c()
for(i in 1:length(extraSamples)){
  newRow = c(# Row number
    as.numeric(gourraud_df$Col[nrow(gourraud_df)]) + 1,
    # Subject
    extraSamples[i],
    # Sample (= Subject identifier)
    extraSamples[i],
    # Use "De Bakker, NA" as population
    "De Bakker, NA",
    # Pedigree not defined
    NA
  )
  # ... Next columns in "newRow" are for the HLA alleles
  
  # Fill in from De Bakker 1 if possible, otherwise fill in from De Bakker 2
  if(extraSamples[i]%in%debakker1$Sample){
    # First 8 columns (A,B,C,DRB1)
    newRow = c(newRow, unname(unlist(debakker1[which(debakker1$Sample==extraSamples[i]),c(2:9)])))
    newRow = c(newRow,
               # After (A,B,C,DRB1) there are 3 genes (6 columns)
               # that we don't have information about
               # (DRB3,4,5)
               rep("",6),
               # Fill in other genes in correct order
               # (order of DQA1 and DQB1 differs between Gourraud and De Bakker)
               unname(unlist(debakker1[which(debakker1$Sample==extraSamples[i]),c(12,13,10,11)])))
  } else {
    # Same reasoning for De Bakker 2 data
    newRow = c(newRow, unname(unlist(debakker2[which(debakker2$Sample==extraSamples[i]),c(2:9)])))
    newRow = c(newRow, rep("",6),unname(unlist(debakker2[which(debakker2$Sample==extraSamples[i]),c(12,13,10,11)])))
  }
  # Add this information to the data frame
  # The next 27 (= 5 + 22) we do not have information about
  gourraud_df = rbind(gourraud_df, c(newRow,rep("",5),rep(NA,22)))
}

# Write new gold standard to file
write.table(gourraud_df, "temp/GourroudAndDeBakker_gold_standard.txt",
            quote = FALSE, append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
