library(tidyverse)

# HLA naming at the time of publication of Adams 2005
old_naming = read_delim(
  "downloads/HLA_nomenclature/Allelelist.2090.txt",
  delim = " ",
  col_names = c("AlleleID", "Allele"),
  col_types = cols(AlleleID = col_character(),
                   Allele = col_character())
)

new_naming = read_delim(
  "downloads/HLA_nomenclature/Allelelist.3440.txt",
  delim = ",",
  comment = "#",
  col_types = cols(
    AlleleID = col_character(),
    Allele = col_character()
  )
)

allele_map_fullres = left_join(old_naming, new_naming, by="AlleleID", suffix=c(".old", ".new")) %>%
  select(Allele.old, Allele.new) %>%
  tidyr::separate("Allele.old", c("gene.old", "allele.old"), sep="\\*") %>%
  tidyr::separate("Allele.new", c("gene.new", "allele.new"), sep="\\*")

# Manually add custom mappings between alleles based on information from IMGT-HLA
# Alleles that could not be mapped accurately are listed in "problematic_mappings" below
# When such an allele appears in the data, a solution is provided here.

# 020204 was deleted from the database as the sequence corresponds to C*02:10:01:01
# Source: https://www.ebi.ac.uk/cgi-bin/ipd/imgt/hla/get_allele.cgi?Cw*020204
# We also need to remap alleles at a higher resolution than 4 digits,
#  sometimes these alleles are mapped to a different allele at a lower resolution
#  e.g. 02:02:04 -> 02:10:01
custom_mappings = tribble(~ gene.old, ~ allele.old, ~ gene.new, ~ allele.new,
                          "Cw", "020204", "C", "02:10:01:01",
                          "A", "021701", "A", "02:17:02:01",
                          "B", "47010101", "B", "47:01:01:03")

# Replace the (gene.old, allele.old) pairs for which a manual mapping was defined
allele_map_fullres = bind_rows(anti_join(
  allele_map_fullres,
  custom_mappings,
  by = c("gene.old", "allele.old")
),
custom_mappings)

list_not_one_to_one <- function(df) {
  df %>%
    distinct %>%
    group_by(gene.old, allele.old) %>%
    add_count() %>%
    filter(n != 1 | any(is.na(allele.new))) %>%
    distinct(gene.old, allele.old)
}

# Ensure the mapping is always 1 -> 1
# List all (old) alleles that do not exist anymore in the database
# + old alleles that are mapped to multiple new alleles
problematic_mapping_fullres = list_not_one_to_one(allele_map_fullres)

allele_map_4d = allele_map_fullres %>%
  mutate(allele.old=str_replace(allele.old, "([0-9][0-9])([0-9][0-9]).*", "\\1\\2"),
         allele.new=str_replace(allele.new, "([0-9]+):([0-9]+).*", "\\1:\\2")) %>%
  distinct %>%
  # Filter for the genes we need (MICA and MICB cause problems)
  filter(gene.old %in% c("A", "B", "Cw", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"))

# List for which alleles there is still no 1 -> 1 mapping at 4 digit level
problematic_mapping_4d = list_not_one_to_one(allele_map_4d)

# The mapping at 4 digit level is not unique for HLA-C
#   Allele.old gene.old allele.old Allele.new    gene.new allele.new
# 1 Cw*020201  Cw       0202       C*02:02:01    C        02:02     
# 2 Cw*020205  Cw       0202       C*02:10:06    C        02:10     
# 3 Cw*070101  Cw       0701       C*07:01:01:01 C        07:01     
# 4 Cw*070103  Cw       0701       C*07:18:04    C        07:18

read_html_table <- function(path) {
  xml2::read_html(path) %>%
    rvest::html_element(css = ".data") %>%
    rvest::html_table()
}

split_allele_col <- function(col) {
  # Get current column name
  colname = colnames(col)
  # Remove "Locus" from column name
  gene = str_remove(colname, regex(' Locus$', ignore_case=TRUE))
  
  # Rename DRβ1 to DRB1
  if (gene == "DRβ1") {
    gene = "DRB1"
  }

  # Split into two allele columns
  tidyr::separate(col, colname, paste0("c", c(1, 2)), ', *', fill='right') %>%
    # If no second allele is specified, assume that the sample is homozygous for the first allele
    mutate(c2=if_else(is.na(c2), c1, c2)) %>%
    # Add gene name to the column name
    rename_with(~paste0(gene, '.', 1:2))
}

clean_allele_name <- function(col) {
  gather(col, "key", "allele.old") %>%
    tidyr::separate("key", c("gene.old", "idx"), sep="\\.", remove=F) %>%
    # Ensure at least 4 digit resolution
    # Replace all alleles that have less than 4 digit resolution by NA
    mutate(allele.old=if_else(str_detect(allele.old, '^([0-9]{4}).*'), allele.old, NA_character_)) %>%
    # Convert pre-2010 naming to current naming
    # Some alleles at 6 digit level (in 2005) correspond to an allele with a different 4 digit part in 2021
    # (This is the case e.g. for HLA-C)
    # First try to map based on the full resolution
    # TODO: Check whether mapping is one-to-one for full resolution conversion
    left_join(allele_map_fullres, by=c("gene.old", "allele.old")) %>%
    mutate(rownum=row_number()) %>%
    group_by(!is.na(allele.new)) %>%
    group_modify(function(df, state) {
      if (state == FALSE) {
        df_4d_trimmed = df %>%
          # Remove previous (failed) mapping
          select(-ends_with('.new')) %>%
          # Perform 4 digit mapping
          mutate(allele.old=str_replace(allele.old, '^([0-9]{4}).*', '\\1'))
        
        # Check whether the mapping at 4 digit level is problematic
        df_4d_trimmed %>%
          semi_join(problematic_mapping_4d, by = c("gene.old", "allele.old")) %>%
          nrow %>%
          walk(~ stopifnot(. == 0))
        
        # Perform the actual mapping
        df_4d_trimmed %>%
          left_join(allele_map_4d, by = c("gene.old", "allele.old"))
      } else {
        df
      }
    }) %>%
    # Remove grouping
    ungroup %>%
    select(-`!is.na(allele.new)`) %>%
    # Restore original order
    arrange(rownum) %>%
    # Limit resolution to 4 digits (ignoring trailing letters)
    mutate(allele.new=str_replace(allele.new, '^([0-9]+:[0-9]+).*', '\\1')) %>%
    transmute(rownum, key=paste0(gene.new, '.', idx), value=allele.new) %>%
    spread('key', 'value') %>%
    select(-rownum, -starts_with('NA'))
}

tbl.orig = Sys.glob("downloads/pub/adams_2005/MHC*.html") %>%
  as_tibble_col("path") %>%
  # Parse the HTML pages
  mutate(html_obj=map(path, read_html_table)) %>%
  # Join the two tables (MHC-I and MHC-II) based on their cell line ID 
  pull(html_obj) %>%
  reduce(full_join, by=c("ID", "Cell Line", "Tissue")) %>%
  # Split allele columns
  transmute(
    # Keep all columns without allele information
    select(cur_data(), !ends_with("Locus")),
    # Split all columns with allele information into two columns
    alleles = lmap(select(cur_data(), ends_with("Locus")), split_allele_col)
  )

tbl = tbl.orig %>%
  mutate(lmap(alleles, clean_allele_name)) %>%
  select(-alleles)

saveRDS(tbl, "downloads/pub/adams_2005/hla_calls.rds")