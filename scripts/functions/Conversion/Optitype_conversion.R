# Read in the output and store as R dataframe
toolOutputToR.Optitype <- function(outputFolder) {
  # Get a list of file names in the given outputFolder
  fileList <- list.files(path = outputFolder,
                         pattern = "_result.tsv$",
                         full.names = T)
  res = vroom::vroom(fileList, id="filename", show_col_types = FALSE, .name_repair = "minimal") %>%
    select(-2, -Reads, -Objective) %>%
    mutate(sample_id=str_remove(basename(filename), "_([^_]+_){6}result\\.tsv$"),
           filename=NULL) %>%
    select(sample_id, everything()) %>%
    mutate(across(-sample_id, ~str_remove(., "^[^*]+\\*")))
  
  colnames(res) <- str_remove(colnames(res), "[0-9]+$")
  # Use classical dataframe with rownames
  res = column_to_rownames(res, "sample_id")
  res
}