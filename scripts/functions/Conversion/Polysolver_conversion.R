# Read in the output and store as R dataframe
toolOutputToR.Polysolver <- function(outputFolder) {
  # Get a list of file names in the given outputFolder
  fileList <-
    list.files(path = outputFolder,
               pattern = "winners.hla.nofreq.txt$",
               full.names = T)
  
  res = fileList %>%
    vroom::vroom(col_names=c("gene", "1", "2"),
                 id = "filename",
                 show_col_types = FALSE) %>%
    mutate(sample_id=str_remove(basename(filename), "_winners.hla.nofreq.txt$"),
           gene=str_remove(gene, "^HLA-"),
           filename=NULL) %>%
    gather("key","value", -gene,-sample_id) %>%
    mutate(value=str_replace_all(str_remove(value, "^hla_[^_]+_"), "_", ":")) %>%
    mutate(key=paste0(gene,key)) %>%
    select(sample_id, key, value) %>%
    spread("key", "value")
  
  # Force non-unique column names as this is required by the scripts of the students
  colnames(res) <- str_remove(colnames(res), "[0-9]+$")
  # Use classical dataframe with rownames
  res = column_to_rownames(res, "sample_id")
  res
}