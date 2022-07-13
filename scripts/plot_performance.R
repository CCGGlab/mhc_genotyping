library(tidyverse)
library(patchwork)

source("scripts/functions/radar_plot.R")

# Functions
source("scripts/functions/extract_dataset.R")

arrange_in_panel <- function(plts) {
  # Set top margin to 0 for plots on second row
  for (i in 3:4) {
    plts[[i]]$theme$plot.margin[[1]] = unit(0, 'cm')
  }
  
  pnl <- (
    (plts[[1]] | plts[[2]]) /
      (plts[[3]] | plts[[4]])
  )
  
  grb <- patchworkGrob(pnl)
  grb2 <- gtable::gtable_add_grob(grb, grid::textGrob("DNA", rot=90,
                                                      gp = grid::gpar(fontsize = 24),
                                                      x = unit(1, 'npc'),
                                                      hjust = 0.5, vjust=0),
                                  # Center of row 1 panel
                                  t = 2, b = 19, l=1, r=1) %>%
    gtable::gtable_add_grob(grid::textGrob("RNA", rot=90,
                                           gp = grid::gpar(fontsize = 24),
                                           x = unit(1, 'npc'),
                                           hjust = 0.5, vjust=0),
                            # Panel object of row 2
                            # (somehow needed to center the label correctly)
                            t = 29, b = 29, l=1, r=1)
  # For first row only
  grb2$grobs[[2]] <- grb2$grobs[[2]] %>%
    # Add row of tag, but column of panel
    gtable::gtable_add_grob(grid::textGrob("MHC-I",
                                           gp = grid::gpar(fontsize = 24),
                                           x = unit(0.5, 'npc'),
                                           y = unit(1, 'npc'),
                                           hjust = 0.5, vjust=1),
                            t = 2, b = 2, l=8, r=8) %>%
    gtable::gtable_add_grob(grid::textGrob("MHC-II",
                                           gp = grid::gpar(fontsize = 24),
                                           x = unit(0.5, 'npc'),
                                           y = unit(1, 'npc'),
                                           hjust = 0.5, vjust=1),
                            t = 2, b = 2, l=23, r=23)
  
  grb2$widths[[1]] = unit(48, 'pt')
  grb2$grobs[[2]]$heights[[2]] = unit(48, 'pt')
  grb2
}

get_mhc_class <- function(df, gene_col) {
  mutate(df, mhc_class = case_when(
    .data[[gene_col]] %in% c("A", "B", "C") ~ "I",
    .data[[gene_col]] %in% c("DPA1", "DPB1", "DQA1", "DQB1", "DRB1") ~ "II"
  ))
}

get_sample_count <- function(evaluated_df, gs, data_type) {
  # Note: we have many more samples in the gold standard than were actually downloaded
  available_samples <- evaluated_df %>%
    #unnest(iscorrect) %>%
    filter(.data[["data_type"]] == .env[["data_type"]]) %>%
    distinct(sample_id)

  sample_count_available <- gs %>%
    filter(!is.na(allele)) %>%
    semi_join(available_samples, by = "sample_id") %>%
    group_by(sample_id, gene) %>%
    count() %>%
    filter(n == 2) %>%
    group_by(gene) %>%
    count()

  sample_count_available
}

evaluated_real <- readRDS("data/iscorrect_tools.rds")
metrics_real <- readRDS("data/performance_tools.rds")

evaluated_meta <- readRDS("data/iscorrect_meta_compact.rds")
metrics_meta <- readRDS("data/performance_meta_compact.rds")

for (dataset_name in c("1kg", "NCI-60")) {
  evaluated_df <- extract_dataset(dataset_name, evaluated_real)
  metrics_df <- extract_dataset(dataset_name, metrics_real)
  
  if (dataset_name %in% names(metrics_meta)) {
    evaluated_dfm <- extract_dataset(dataset_name, evaluated_meta) %>%
      filter(data_type == "DNA-Seq") %>%
      filter(tool == "meta") %>%
      mutate(tool = "Meta")
    metrics_dfm <- extract_dataset(dataset_name, metrics_meta) %>%
      filter(data_type == "DNA-Seq") %>%
      filter(tool == "meta") %>%
      mutate(tool = "Meta")
    
    evaluated_df <- bind_rows(evaluated_df, evaluated_dfm)
    metrics_df <- bind_rows(metrics_df, metrics_dfm)
  }
  
  # ===
  # Determine GS sample size per gene
  # Read gold standard data for 1000 genomes
  # Convert to long table format
  gs_name <- list("1kg"="1kg", "NCI-60"="nci60")[[dataset_name]]
  gs <- readRDS(str_glue("data/gold_standard_{gs_name}.rds")) %>%
    gather("gene", "allele", -sample_id) %>%
    tidyr::separate(gene, c("gene", "idx"), sep = "\\.")
  
  # Use exact tool names and case styling used by authors
  formatted_tool_name = c(HLAminerHPRA = "HLAminer",
                          HLAVBseq = "HLA-VBSeq",
                          kourami = "Kourami",
                          `HLA-LA` = "HLA*LA")
  
  # Format table and write to results folder
  metrics_formatted <- metrics_df %>%
    get_mhc_class("gene") %>%
    # HLA typing on NCI-60 data makes no sense for MHC-II
    filter(!(dataset_name == "NCI-60" & data_type == "RNA-Seq" & mhc_class == "II")) %>%
    # For some genes no information was present in the gold standard
    filter(gene %in% unique(gs$gene)) %>%
    # Polysolver and Optitype are not compatible with MHC-II
    # Do not display a dot for these tools
    filter((mhc_class != "II") | !(tool %in% c("Polysolver", "Optitype"))) %>%
    # Proportion correct relative to the total number of predictions
    mutate(decided_ratio = (1 - undecided_ratio)) %>%
    mutate(decided_correct = decided_ratio * accuracy) %>%
    mutate(tool = if_else(tool %in% names(formatted_tool_name), unname(formatted_tool_name[tool]), tool))
  
  gene_label <- function(data_type) {
    # Use official gene names (with HLA prefix) and add sample count
    data_type <- unique(data_type)
    sample_counts <- get_sample_count(evaluated_df, gs, data_type)
    function(x) {
      left_join(as_tibble_col(x),
                sample_counts,
                by = c("value" = "gene")) %>%
        str_glue_data("HLA-{value}\n(n = {n})")
    }
  }
  
  plots <- metrics_formatted %>%
    group_by(data_type, mhc_class) %>%
    summarise(plt = list(radar_plot(
      cur_data(),
      x_var = "tool",
      y_var = "decided_correct",
      color_var = "gene",
      labeller = gene_label(data_type) ) +
        theme(plot.margin = margin(2,3.5,2,2, unit = "cm"),
              legend.text = element_text(hjust = 0.5, margin = margin(r = 15, unit = "pt"))
        )
        )
    )
  
  # Supplementary figures
  supp_plots_decided <- metrics_formatted %>%
    # Do not display unsupported genes as "undecided"
    drop_na %>%
    group_by(data_type, mhc_class) %>%
    summarise(plt = list(
      radar_plot(
        mutate(cur_data(), decided_ratio = 1 - undecided_ratio),
        x_var = "tool",
        y_var = "decided_ratio",
        color_var = "gene",
        labeller = gene_label(data_type)
      ) +
        theme(plot.margin = margin(2,3,2,2, unit = "cm"),
              legend.text = element_text(hjust = 0.5, margin = margin(r = 15, unit = "pt"))
        )
    ), .groups = "drop")
  
  supp_plots_correct <- metrics_formatted %>%
    group_by(data_type, mhc_class) %>%
    summarise(plt = list(
      radar_plot(
        mutate(cur_data(), decided_ratio = 1 - undecided_ratio),
        x_var = "tool",
        y_var = "accuracy",
        color_var = "gene",
        labeller = gene_label(data_type)
      ) +
        theme(plot.margin = margin(2,3,2,2, unit = "cm"),
              legend.text = element_text(hjust = 0.5, margin = margin(r = 15, unit = "pt"))
              )
    ), .groups = "drop")

  # Arrange plots in panel for publication 
  
  # Quick fix for missing panel with NCI-60
  plts <- plots$plt
  if (length(plts) == 3) {
    plts[[4]] <- zeroGrob()
  } 
  pnl_grb <- arrange_in_panel(plts)
  # Save plot
  ggsave(str_glue("results/figs/accuracy_radar_{dataset_name}_panel.pdf"), pnl_grb, height=21, width=21, scale=2, units = "cm")
  
  # Arrange plots in panel for publication with cowplot
  plts <- supp_plots_decided$plt
  if (length(plts) == 3) {
    plts[[4]] <- zeroGrob()
  } 
  pnl_grb <- arrange_in_panel(plts)
  ggsave(str_glue("results/figs/supplementary_{dataset_name}_decided.pdf"), pnl_grb, height=21, width=21, scale=2, units = "cm")
  
  # Arrange plots in panel for publication with cowplot
  plts <- supp_plots_correct$plt
  if (length(plts) == 3) {
    plts[[4]] <- zeroGrob()
  } 
  pnl_grb <- arrange_in_panel(plts)
  ggsave(str_glue("results/figs/supplementary_{dataset_name}_correct.pdf"), pnl_grb, height=21, width=21, scale=2, units = "cm")
}
