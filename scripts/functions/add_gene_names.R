add_gene_names <- function(df, cols=sample_id) {
  cols = substitute(cols)
  df %>%
    dplyr::ungroup() %>%
    # Apply the following function to all columns except for the column "sample_id"
    dplyr::mutate(dplyr::across(!!cols, function(val) {
      # If value is NA: return NA, else add the gene name to the allele name, separated by star
      # The gene name is derived from the current column name, by removing everything after the last dot
      dplyr::if_else(is.na(val), NA_character_, paste0(stringr::str_remove(dplyr::cur_column(), '\\..*$'), '*', val))
    })) %>%
    # Restore original grouping
    group_by(!!!groups(df))
}
