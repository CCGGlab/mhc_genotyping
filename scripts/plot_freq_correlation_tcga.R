library(tidyverse)
# stop("Rename tools according to author names")
# Required: ggtext, latex2exp, Metrics
# Load functions
source("scripts/functions/extract_dataset.R")
source("scripts/functions/calc_undecided.R")
source("scripts/functions/add_gene_names.R")

formatted_tool_name = c(HLAminerHPRA = "HLAminer",
                        HLAVBseq = "HLA-VBSeq",
                        kourami = "Kourami",
                        `HLA-LA` = "HLA*LA")

# Import custom plot elements
source("scripts/functions/custom_plot_functions.R")

# Read in AFN: use consensus.rds data
# Both prediction and consensus frequencies were mapped to G-groups
# (Allele frequencies often reported at G-group resolution due to PCR-SSOP)
consensus_frequencies <- readRDS("downloads/AFN/consensus_ggroup.rds")

# Add clinical data: see check_tcga_hla_hd (necessary for ethnicity)
metadata <- read.table("downloads/TCGA/HLA/metadata.tsv",
  sep = "\t",
  header = T
)

ethnicity_info <- metadata %>%
  distinct(submitter_id = cases.0.submitter_id, race = cases.0.demographic.race)

# Check whether submitter_id -> race mapping is unique
ethnicity_info %>%
  group_by(submitter_id) %>%
  count %>%
  filter(n > 1)
# => OK

results <- tibble::tibble()
results_undecided <- tibble::tibble()

# Read in result dataframe for TCGA (mapped to G-groups)
predictions_all_datasets <- readRDS("data/hla_predictions_ggroup.rds")
predictions_all_meta <- readRDS("data/hla_predictions_ggroup_meta.rds") %>%
  chuck("compact")

predictions_tcga <- extract_dataset("TCGA", predictions_all_datasets)
# Only results of DNA metaclassifier are included for now
predictions_meta <- extract_dataset("TCGA", predictions_all_meta) %>%
  filter(tool == "meta") %>%
  mutate(tool = "Meta") %>%
  filter(data_type == "DNA-Seq")

predictions_tcga <- predictions_tcga %>%
  bind_rows(predictions_meta)

predictions <- predictions_tcga %>%
  mutate(tool = if_else(tool %in% names(formatted_tool_name), unname(formatted_tool_name[tool]), tool)) %>%
  # Convert to long table format, keep metadata as columns
  gather("gene", "allele", -sample_id, -tool, -data_type) %>%
  # Separate gene and index
  tidyr::separate(gene, c("gene", "idx"), sep = "\\.") %>%
  group_by(data_type, tool, sample_id, gene) %>%
  # Sort the alleles (make sure the order of the alleles per sample is fixed)
  arrange(allele, .by_group=T) %>%
  mutate(idx = seq_along(allele)) %>%
  # Shorten sample_id to case submitter ID
  mutate(submitter_id = str_sub(sample_id, 1, 12)) %>%
  # Subject submitter_id is not unique
  # Only consider 1 sample per subject
  pivot_wider(names_from = "idx", names_prefix = "a", values_from = "allele") %>%
  group_by(across(-c(sample_id, a1, a2))) %>%
  # Remove records where one or both alleles are not predicted
  drop_na %>%
  # Remove duplicate rows: select first portion / vial
  slice_max(sample_id) %>%
  # Convert back to long table format
  pivot_longer(cols = c(a1, a2), names_prefix = "a", names_to = "idx", values_to = "allele")

nrow(predictions)

predictions <- predictions %>%
  # Add column for ethnicity
  # (Removing rows where ethnicity information is not available)
  inner_join(ethnicity_info, by = "submitter_id")

nrow(predictions)

# How many samples (aliquots) remaining?
predictions %>%
  ungroup %>%
  distinct(data_type, sample_id) %>%
  group_by(data_type) %>%
  count
# => Samples with no prediction are deleted now (drop_na in "predictions")

# How many subjects remaining?
predictions %>%
  ungroup %>%
  distinct(data_type, submitter_id) %>%
  group_by(data_type) %>%
  count

# How many samples left per ethnicity
samples_per_ethnicity <- predictions %>%
  ungroup %>%
  distinct(data_type, submitter_id, race) %>%
  group_by(data_type, race) %>%
  count

print("Samples per ethnicity:")
samples_per_ethnicity
predictions %>%
  ungroup %>%
  distinct(submitter_id, race) %>%
  group_by(race) %>%
  count

predictions %>%
  ungroup %>%
  distinct(data_type, tool, submitter_id, race) %>%
  group_by(data_type, tool, race) %>%
  count

# Check whether every sample is in a category (including "not reported" or "")
samples_per_ethnicity %>%
  group_by(data_type) %>%
  summarise(sum(n))

# # The code block below is not valid:
# # - Needs to be executed before NA predictions are removed
# # - What about tools that were not run on all samples?
# predictions_undecided <- predictions %>%
#   group_by(tool, data_type) %>%
#   calc_undecided() %>%
#   mutate(decided_ratio = 1 - undecided_ratio, undecided_ratio = NULL)

# Calculate predicted frequency per ethnicity for each tool
predictions_frequencies <- predictions %>%
  # We require information about ethnicity
  filter(race != "not reported", race != "") %>%
  # Rename factors for race to terminology used for consensus population frequencies
  mutate(
    race = recode(
      race,
      "white" = "Caucasian American",
      "black or african american" = "African American"
    ) %>%
      as.character()
  ) %>%
  # Only consider made predictions
  # Predictions that cannot be made do not count for the total for that tool, gene
  filter(!is.na(allele)) %>%
  # Calculate frequency of HLA genes per ethnicity
  # 1) Calculate number of times each allele was predicted (per race, gene)
  group_by(tool, data_type, race, gene, allele) %>%
  count() %>%
  # 2) Count total number of predictions for each gene (per race)
  group_by(tool, data_type, race, gene) %>%
  mutate(tot = sum(n)) %>%
  # 3) Relative frequency = (1) / (2)
  mutate(freq = n / tot) %>%
  ungroup()

predictions_frequencies %>%
  saveRDS(str_glue("data/predicted_freqs_tcga_ggroup.rds"))

# Create dataframe containing both predicted and expected frequencies
df_joined <- predictions_frequencies %>%
  # Rename the frequency column to avoid confusion
  rename("predicted.frequency" = freq) %>%
  # Add rownumber: for technical purposes:
  # Will be NA in join if an allele was not included.
  tibble::rowid_to_column("rowid.pred") %>%
  # Filter for ethnicities where we have consensus frequencies
  filter(race %in% c("African American", "Caucasian American")) %>%
  # Add column for population consensus frequencies
  full_join(
    gather(
      consensus_frequencies,
      "race",
      "population.frequency",
      -gene,
      -allele
    ),
    by = c("gene", "allele", "race")
  ) %>%
  # It makes sense to set the frequency of an allele that was never predicted on 9000 samples
  # to zero for the NGS methods
  # Contrary, for the consensus frequencies we cannot do this!
  # With PCR-SSOP it is possible certain alleles were never tested for.
  # => If allele found in population data was never predicted, set 0 as predicted frequency
  mutate(predicted.frequency = if_else(is.na(rowid.pred), 0, predicted.frequency)) %>%
  select(-rowid.pred) %>%
  # Tool column can be NA if an allele was predicted by no tool for a given ethnicity
  # Remove these records as this is nonsensical for our purpose.
  # TODO: reconsider. Instead include these records for all tools with freq = 0?
  filter(!is.na(tool))

calc_metrics <- function(df, key) {
  # TODO: handle case where prediction or population variance = 0

  # If less than 3 records are available with known consensus frequency,
  # we cannot calculate a correlation
  if (sum(!is.na(df$population.frequency)) < 3) {
    correlation <- tibble::tibble(correlation.pearson = NA_real_)
  } else {
    # TODO: explicitly handle NAs for calculation of population frequency
    # -> discard population NAs?
    correlation <- broom::tidy(cor.test(
      df$predicted.frequency,
      df$population.frequency
    )) %>%
      select(correlation.pearson = estimate, p.pearson = p.value)
  }

  root_mean_squared_error <- df %>%
    filter(!is.na(population.frequency)) %>%
    with(Metrics::rmse(population.frequency, predicted.frequency))

  ff_plot <- ggpubr::ggscatter(
    df,
    "predicted.frequency",
    "population.frequency",
    label = "allele",
    repel = T,
    title = paste(unique(key$gene))
  ) +
    geom_abline(slope = 1, intercept = 0) +
    xlab("predicted frequency") +
    ylab("expected frequency")

  tibble::tibble(correlation, RMSE = root_mean_squared_error, plt = list(ff_plot))
}

# Create heatmaps
# normal rectangle
# rows: tools, reference
# cols: genes (duplicated: next to each other for CA and AA)
# value: correlation for
# triangular shape
tcga_metrics <- df_joined %>%
  group_by(tool, data_type, race, gene) %>%
  group_modify(calc_metrics)

saveRDS(tcga_metrics, "data/freq_based_metrics.rds")

heatmap_df <- tcga_metrics %>%
  select(tool, data_type, race, gene, p.pearson, correlation.pearson) %>%
  # Add NA for all missing combinations
  # cross_df: data frame with all combinations
  # map(select(...), unique): for tool, data_type, race, gene
  # -> give unique values
  right_join(cross_df(map(select(., tool, data_type, race, gene), unique))) %>%
  # Filter out what should not be displayed
  # Discard for which we have all NAs for a given population
  group_by(race, gene) %>%
  filter(!all(is.na(correlation.pearson))) %>%
  # Discard tools for the data types they do not support
  group_by(tool, data_type) %>%
  filter(!all(is.na(correlation.pearson))) %>%
  ungroup() %>%
  # Formatting
  mutate(data_type = recode(data_type, `RNA-Seq` = "RNA", `DNA-Seq` = "DNA"))

# Store interesting scatter plots as a supplementary figure
supp_plots <- list(
  tcga_metrics %>% filter(tool == "PHLAT", data_type == "DNA-Seq", gene == "A", race == "African American"),
  tcga_metrics %>% filter(tool == "PHLAT", data_type == "DNA-Seq", gene == "A", race == "Caucasian American"),
  tcga_metrics %>% filter(tool == "PHLAT", data_type == "RNA-Seq", gene == "A", race == "African American"),
  tcga_metrics %>% filter(tool == "PHLAT", data_type == "RNA-Seq", gene == "A", race == "Caucasian American"),

  tcga_metrics %>% filter(tool == "arcasHLA", data_type == "RNA-Seq", gene == "DRB1", race == "African American"),
  tcga_metrics %>% filter(tool == "arcasHLA", data_type == "RNA-Seq", gene == "DRB1", race == "Caucasian American")
) %>% map("plt") %>% flatten

supp_panel = cowplot::plot_grid(plotlist = supp_plots,
                                align='v', ncol=2, labels = "AUTO", greedy=F, label_size=24)
ggsave(str_glue("results/figs/supplementary_tcga_panel.pdf"), supp_panel, width=21, scale=2, units = "cm")

heatmap_plt <- ggplot(
  heatmap_df,
  aes(
    x = gene,
    # Row labels are ordered by mean Pearson correlation, ...
    # ... ignoring NA values.
    y = reorder(fct_cross(tool, data_type), correlation.pearson, mean, na.rm = T),
    ylabel = tool,
    ycolour = data_type,
    colour = correlation.pearson,
    size = p.pearson
  )
) +
  facet_custom(~ race, scales = "free_x") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.text.y = ggtext::element_markdown(size = 18),
    strip.text = element_text(size = 18),
    panel.spacing = unit(1, "lines"),
    panel.grid = element_blank(),
    strip.background = element_bottom_line(),
    legend.text = element_text(size = 18),
    legend.title.align = 0.1,
    # legend.position = "bottom",
    # legend.box = "horizontal",
    text = element_text(size = 18),
    axis.title = element_blank()
  ) +
  labs(
    #colour = latex2exp::TeX("$\\rho$"),
    colour = "Pearson r",
    size = latex2exp::TeX("$\\mathit{p}\\mathrm{-value}$"),
    ycolour = "data type"
  ) +
  scale_x_discrete(labels = ~ paste0("HLA-", .)) +
  colorspace::scale_colour_discrete_diverging(palette = "Blue-Red", aesthetics = "ycolour", guide = guide_legend(order=2)) +
  geom_custom_legend(custom_aes="ycolour") +
  geom_point(shape = 16) +
  scale_size(
    range = c(2, 11),
    trans = reverse_log_trans(10),
    guide = guide_legend(
      override.aes = aes(shape = 21),
      label.theme = element_text(
        size = 18, # Was too small by default. Make sure this equals legend.text
        hjust = 0,
        margin = margin(l = -12)
      )
    ),
    labels = size_labels
  ) +
  colorspace::scale_colour_continuous_sequential(
    # Blues 3: few contrast, better candidate for linear scale
    # Hawaii: without reversing and exponential scale (base 10^5)
    # -> good contrast between good and bad tools,
    # also contrast within the good tools
    palette = "OrYel",
    rev = T,
    na.value = "grey50",
    trans = scales::exp_trans(base = 100),
    limits = c(0, 1),
    breaks = c(0.5, 0.9, 0.95, 1),
    guide = guide_colorbar(order = 1)
  )

# res_plt <- heatmap_plt +
  # theme(plot.margin = margin(0, 4, 0, 4, 'cm'))

ggsave(
  "results/figs/tcga_heatmap_with_meta.pdf",
  heatmap_plt + theme(plot.margin=margin(l=2,r=2,unit="cm")),
  scale = 2,
  width = 13.86, # space between left and right margin for IJMS
  height = 11,
  units = "cm"
)
# Also save as SVG (Inkscape conversion was not perfect)
ggsave(
  "results/figs/tcga_heatmap_with_meta.svg",
  heatmap_plt + theme(plot.margin=margin(l=2,r=2,unit="cm")),
  scale = 2,
  width = 13.86, # space between left and right margin for IJMS
  height = 11,
  units = "cm"
)

# TODO: rename tools according to author names

# Conclusions:
# DQA1 05:05 vs 05:01 issue seems to be solved by
# switch to G-groups (both are in same G-group)
#   + use of consensus frequencies
# - All well performing tools have a very high correlation
# - Poor performing tools also perform poor here.
