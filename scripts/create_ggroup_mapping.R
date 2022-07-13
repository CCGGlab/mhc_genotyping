library(tidyverse)

trim_genes <- function(x) {
  str_replace(x, "^([0-9]+:[0-9]+).*", "\\1") %>%
    # Ensure the same allele is always chosen as the group designation
    # Sort numerically on digits
    # (Prevents problems when multiple 3 digits are present at a single position)
    as_tibble_col("val") %>%
    tidyr::separate(sep=":", val, into=c("v1", "v2"), remove=F) %>%
    mutate(across(c(v1,v2), as.numeric)) %>%
    arrange(v1, v2) %>%
    pull(val) %>%
    unique
}

# G-groups
g_groups = readLines("downloads/HLA_nomenclature/hla_nom_g.txt") %>%
  discard(~str_detect(.,'^ *#')) %>%
  as_tibble_col("col") %>%
  tidyr::separate("col", c("locus", "alleles", "col3"), sep=';') %>%
  # Very few alleles are in the third column: ignore this
  # Difficult to know what this represents without a header ...
  select(-col3) %>%
  mutate(alleles=str_split(alleles, '/')) %>%
  # Ignore "G-groups" containing a single allele.
  # These cause problems when these alleles are simultaneously mentioned in another G-group.
  filter(map_lgl(alleles, ~length(.) > 1)) %>%
  # Trim to 4 digits, remove additional info at the end
  mutate(alleles=map(alleles, trim_genes)) %>%
  # Name the G-groups according to the first allele that is member of that group
  mutate(key=map_chr(alleles, ~.[[1]])) %>%
  unnest(alleles) %>%
  transmute(gene=str_replace(locus, '\\*$', ''), group=paste0(locus, key), allele=paste0(locus, alleles)) %>%
  distinct

# Ensure the mapping allele -> group is unique
g_groups %>%
  group_by(allele) %>%
  summarise(l=list(unique(group))) %>%
  filter(map_lgl(l, ~length(.) > 1))

g_groups %>%
  filter(allele == "DPA1*02:12")

# Remark:
# DPA1*02:12 cannot be uniquely mapped to a G-group
# DPA1*02:12:02 is assigned to DPA1*02:02:02G
# DPA1*02:12:01 is assigned to DPA1*02:07:01G
# => Discard the mapping to DPA1*02:07 (which further contains no other alleles than DPA1*02:07 itself)
g_groups = g_groups %>%
  filter(!(allele == "DPA1*02:12" & group == "DPA1*02:07"))

# Ensure the mapping allele -> group is now unique
g_groups %>%
  group_by(allele) %>%
  summarise(l=list(unique(group))) %>%
  filter(map_lgl(l, ~length(.) > 1))

g_groups %>%
  saveRDS("data/ggroup_mapping.rds")