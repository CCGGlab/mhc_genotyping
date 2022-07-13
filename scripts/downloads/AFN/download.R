library(tidyverse)

# Download population frequencies and metadata for the following list of populations
# Note that the San Francisco HLA types (population 3098) are based on NGS data <-> Sequence Primer
# Software: Connexio Genomics ATF
# This dataset is not used for the consensus frequencies
# Liberia Bong County = 2522

# Population 2779 uses "PCR with sequence-specific primers, performed using ARMS technology, was utilized for high-resolution typing of HLA-DRB1, HLA-DRB3, and HLA-DPB1."
# (was not mentioned on HLA frequency net)
hla_populations = c(1279, 1359, 1401, 1619, 1620, 2419, 2531, 2779, 2780, 1479, 1586, 
                    1588, 1895, 1513, 2570, 1480, 2223, 1511)
mhc1_genes = c("A", "B", "C")
mhc2_genes = c("DPA1", "DPB1", "DQA1", "DQB1", "DRB1")
all_mhc_genes = c(mhc1_genes, mhc2_genes)

# Function was corrected to query per locus.
# This was needed to prevent issues when > 1000 alleles are reported.
get_table_url = function(population_ids, hla_locus) {
  paste0(
    "http://allelefrequencies.net/hla6006b.asp?",
    str_glue(
      "hla_locus={hla_locus}",
      "hla_locus_type=Classical",
      "hla_allele1=",
      "hla_allele2=",
      "hla_selection=",
      "hla_pop_selection=",
      "hla_population={population_ids}",
      "hla_country=",
      "hla_dataset=",
      "hla_region=",
      "hla_ethnic=",
      "hla_study=",
      "hla_sample_size=",
      "hla_sample_size_pattern=",
      "hla_sample_year=",
      "hla_sample_year_pattern=",
      "hla_level=",
      "hla_level_pattern=",
      "hla_show=",
      "hla_order=order_1",
      "standard=",
      .sep = "&"
    )
  )
}

get_info_url <- function(pop_id) {
  str_glue("http://allelefrequencies.net/pop6001c.asp?pop_id={pop_id}")
}

extract_quality <- function(html_source) {
  # The quality code is given in a poorly parseable format
  # This is a best-effort to automatically extract this information.
  
  # Select the second column of the row that contains a column with "Loci typed:"
  xml2::xml_find_first(html_source, "//td[b/text() = 'Loci typed:']/following-sibling::td") %>%
    as.character %>%
    # Delete surrounding <td> and </td> tags
    str_sub(5,-6) %>%
    # Trim string
    str_trim %>%
    # Replace images by a label
    str_replace_all(r"{<img src=\"images/icon_(.+?).gif\">}", " : \\1") %>%
    # Split listed alleles
    str_split(" *, *") %>%
    map_depth(0, 1) %>%
    # Split into locus and label columns
    as_tibble_col %>%
    tidyr::separate("value", c("locus", "quality_code"), " : ", fill = "right") %>%
    spread(locus, quality_code)
}

parse_info_table <- function(html_source) {
  rvest::html_table(html_source, fill=TRUE) %>%
    # Remove first rows from the table (containing wrongly parsed  table headers)
    tail(-2) %>%
    select(key=X2, value=X3) %>%
    # Remove the colon at the end of the key column
    mutate(key=str_replace(key, ":$", ""))
}

parse_info_page <- function(html_source) {
  xml2::xml_find_all(html_source, xpath = "//table[@class = 'table04']") %>%
    # Convert to data frame with a column HTML
    map(identity) %>% as_tibble_col("html") %>%
    # Extract property category from the HTML table
    mutate(title=map(html, compose(xml2::as_list, xml2::xml_find_first), ".//h2/text()[2]")) %>%
    # Remove all tables without a title
    filter(modify(title, length) > 0) %>%
    # Cleanup the title if it exists
    mutate(across(title, str_trim)) %>%
    # Convert all HTML tables to data frames
    mutate(table=map(html, parse_info_table)) %>%
    # Stack all properties
    unnest(table) %>%
    drop_na %>%
    # Replace "Loci typed" by table containing color code
    # mutate(across(value, as.list)) %>%
    # group_by(key) %>%
    # group_modify(function(df, key) {
    #   if (key == "Loci typed") {
    #     df$value = extract_quality(df$html[[1]]) %>% list
    #   }
    #   df
    # }) %>%
    # ungroup %>%
    mutate(map_dfr(html, extract_quality) %>% rename_all(~paste0("quality_", .))) %>%
    # Remove HTML nodes
    select(-html) %>%
    # Put properties into separate columns
    select(-title) %>% spread(key, value)
    # Convert the columns to the most appropriate data type
    # (is only needed if the value column is converted to a list above)
    # mutate(across(everything(), simplify))
}

parse_frequency_page <- function(html_source) {
  node = rvest::html_node(html_source, css = ".tblNormal")
  # If this element is not present on the page,
  # no information for this (dataset, gene) pair was found
  # Return a NULL object
  if (class(node) == "xml_missing") {
    return(NULL)
  }
  node %>%
    rvest::html_table(fill = TRUE) %>%
    select(1:4, frequency := 3) %>%
    # 4th column header contains the name of the population.
    # Remove trailing spaces from the population name.
    mutate(population = colnames(.)[4]) %>%
    # Otherwise column 4 is blank. Remove it.
    select(-4) %>%
    # Population column contains the sample size as well
    tidyr::extract(
      population,
      c("population", "sample_size"),
      r"{(.*?) +\(n=([0-9]+)\)}"
    ) %>%
    # Only keep rows containing a line number
    filter(str_detect(Line, "^[0-9]+$")) %>%
    select(-Line) %>%
    # Split gene and allele name
    tidyr::separate(Allele, c("gene", "allele"), sep = "\\*") %>%
    # Add a column with the resolution
    mutate(resolution = (str_count(allele, ":") + 1) * 2) %>%
    # Cast all columns to the most appropriate type
    type_convert(
      cols(
        gene = col_character(),
        allele = col_character(),
        frequency = col_double(),
        population = col_character(),
        sample_size = col_integer(),
        resolution = col_integer()
      )
    )
}

download_html_cached <- function(url) {
  print(paste("Download", url))
  # Goal: prevent overloading the server by repeatedly requesting the same page
  # (e.g. while debugging this script)
  # TODO: cache does not work
  # cache = paste0(cache_dir, "/", digest::digest(url))
  xml2::read_html(url)
  # Conversion appears to be extremely slow and it still doesn't work.
  # => external pointer is not valid
  
  # use_cached = file.exists(cache)
  # if (use_cached) {
  #   print(paste("Read from cache", url))
  #   html = xml2::as_xml_document(readRDS(cache))
  # } else {
  #   print(paste("Read from web", url))
  #   html = xml2::read_html(url)
  #   saveRDS(xml2::as_list(html), file = cache)
  # }
  # html
}

download_population_metadata <- function(population_ids) {
    tibble(population_id=population_ids) %>%
    mutate(info_url = get_info_url(population_ids),
           info_html = map(info_url, download_html_cached),
           map_dfr(info_html, parse_info_page)
           ) %>%
    # Remove column with HTML code
    select(-info_html)
}

download_population_data_raw <- function(info_df) {
  info_df %>%
    # Repeat dataset per gene and determine URL per (dataset, gene) pair
    # (Only needed as a workaround during download)
    group_by(population_id) %>%
    mutate(gene_download=list(all_mhc_genes)) %>%
    ungroup %>%
    unnest(gene_download) %>%
    # Calculate URLs for (dataset, gene)
    mutate(
      frequency_tbl_url = get_table_url(population_id, gene_download),
      frequency_tbl_html = map(frequency_tbl_url, download_html_cached)
    )
}

process_population_data_raw <- function(freq_df) {
  freq_df %>%
    # Parse downloaded page and replace current dataframe
    transmute(population_id,
              data = map(frequency_tbl_html, parse_frequency_page)) %>%
    # Unpack the data
    unnest(data)
}

sum_to_4d_alleles <- function(freq_df, population_id) {
  res = freq_df %>%
    mutate(allele_4d = if_else(resolution >= 4,
                               str_replace(allele, "^([0-9]+:[0-9]+).*", "\\1"),
                               NA_character_),
           allele_6d = if_else(resolution >= 6,
                               str_replace(allele, "^([0-9]+:[0-9]+:[0-9]+).*", "\\1"),
                               NA_character_)
    )
  
  exact4_4 = res %>%
    filter(resolution == 4) %>%
    pull(allele_4d) %>%
    unique
  
  exact6_4 = res %>%
    filter(resolution == 6) %>%
    pull(allele_4d) %>%
    unique
  
  # Information at 6 digit level, not at 4 digit level
  only_higher = setdiff(exact6_4, exact4_4)
  
  # Limit to extension by 6 digit level to avoid complications
  res2 <- res %>%
    filter(allele_4d %in% only_higher, resolution == 6) %>%
    group_by(allele_4d) %>%
    mutate(allele=allele_4d, allele_6d = NA_character_, frequency=sum(frequency)) %>%
    distinct
  
  bind_rows(res, res2)
}

population_metadata = download_population_metadata(hla_populations)
population_data_raw = download_population_data_raw(population_metadata)
population_data = process_population_data_raw(population_data_raw)

# Post-processing:

results = population_data %>%
  # Note: combination below was not performed eventually.
  # (We are not sure that sum of 6 digit alleles = all 4 digit alleles)
  # # Combine higher resolutions to 4 digit resolution if not explicitly present
  # # Calculation must be performed per dataset (population_id)
  # group_by(population_id) %>%
  # group_modify(sum_to_4d_alleles) %>%
  # ungroup %>%
  # Combine with metadata in a single object
  nest_by(population_id, population, .key="hla_frequencies") %>%
  inner_join(population_metadata, by="population_id")

# TODO: test whether some (population, gene) tuples return more than 1000 records
# -> currently not the case

population_data %>%
  saveRDS("downloads/AFN/hla_allele_frequencies_raw_tmp.rds")

results %>%
  saveRDS("downloads/AFN/hla_allele_frequencies.rds")
