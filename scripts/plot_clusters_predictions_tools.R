library(tidyverse)
library(patchwork)

# Load functions
source("scripts/functions/extract_dataset.R")
source("scripts/functions/dendrogram_coord.R")

# Mapping from tool names as used in predictions to
# names used on plot
formatted_tool_name = c(
  HLAminerHPRA = "HLAminer",
  HLAVBseq = "HLA-VBSeq",
  kourami = "Kourami",
  `HLA-LA` = "HLA*LA"
)

# Load data
is_correct_1kg <-
  extract_dataset("1kg", readRDS("data/iscorrect_tools.rds")) %>%
  mutate(tool = if_else(tool %in% names(formatted_tool_name), unname(formatted_tool_name[tool]), tool))

hierarchical_clustering <- function(long_df) {
  # Convert "is correct" to matrix format
  # (needed for "hclust" cluster method)
  m = spread(long_df, sample_id, correct) %>%
    mutate(across(!tool, as.numeric)) %>%
    # If NA -> replace by incorrect
    mutate(across(!tool, ~ replace_na(., 0))) %>%
    column_to_rownames("tool") %>%
    as.matrix()
  
  # Hierarchical clustering on the tools and samples
  clusters_tools <-
    hclust(dist(m, method = "euclidean"), method = "ward.D")
  clusters_samples <-
    hclust(dist(t(m), method = "euclidean"), method = "ward.D")
  
  # Return the order of the clustered tools / samples
  # and the hclust objects used by ggdendrogram
  list(
    tools = clusters_tools,
    samples = clusters_samples,
    xorder = rownames(m)[clusters_tools$order],
    yorder = colnames(m)[clusters_samples$order]
  )
}

build_plot <- function(is_correct = cur_data()) {
  # Executed per data type and gene
  # Replace "not predicted" by FALSE
  is_correct <- is_correct %>%
    # Convert to wide format
    pivot_wider(id_cols = "sample_id", names_from="tool", values_from = "correct") %>%
    # Replace NA by FALSE
    mutate(across(-sample_id, ~replace_na(., FALSE))) %>%
    # Go back to long format
    pivot_longer(cols = -"sample_id", names_to="tool", values_to = "correct")
  
  # Perform hierachical clustering
  clusters <- hierarchical_clustering(is_correct)
  
  res <- ggplot(is_correct,
                aes(x = tool,
                    y = sample_id,
                    fill = correct)) +
    # Heatmap
    geom_tile() +
    # Order x and y axis according to clusters
    scale_x_discrete(limits = clusters$xorder) +
    scale_y_discrete("samples", limits = clusters$yorder,
                     guide = guide_none()) +
    scale_fill_manual("correct",
      limits = c(FALSE, TRUE),
      values = c("#DE319C", "#4DEB8A"),
      labels = as_mapper( ~ if_else(., "yes", "no"))
    ) +
    # Add a dendrogram for the tools at the top of the plot
    dendrogram_coord(clusters$tools, clusters$samples) +
    # Set plot options
    theme_minimal() +
    theme(
      # Do not display a label for the x-axis
      axis.title.x = element_blank(),
      # Rotate the labels of the y-axis
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      # Left alignment of the legend title
      legend.title.align = 0,
      # Hide the grid
      panel.grid = element_blank(),
      # Add a margin to leave sufficient space between the plots on the panel
      plot.margin = margin(t = 12, l = 24)
    )
  
  # Return plot as single object for summarise function
  list(res)
}

genotype_correct_1kg <- is_correct_1kg %>%
  # Genotype is correct for sample if both alleles are correct
  group_by(data_type, tool, sample_id, gene) %>%
  summarise(correct = all(is_correct))

# Amount of samples wrongly typed by all tools simultaneously
genotype_correct_1kg %>%
  ungroup %>%
  nest_by(data_type) %>%
  deframe %>%
  map(function(df) {
    # To wide format
    pivot_wider(
      df,
      id_cols = c(sample_id, gene),
      names_from = tool,
      values_from = correct
    ) %>%
      # Replace NA by FALSE
      mutate(across(-c(sample_id, gene), ~ replace_na(., FALSE))) %>%
      group_by(sample_id, gene) %>%
      summarise(all_wrong=if_all(everything(), ~!.)) %>%
      group_by(gene) %>%
      # Calculate proportion of samples wrongly predicted by all
      summarise(mean(all_wrong))
  }) %>%
  enframe %>%
  unnest(value) %>%
  group_by(name) %>%
  # Percentage
  summarise(median(`mean(all_wrong)`) * 100)
  

plts <- genotype_correct_1kg %>%
  # Make plot per data type and gene
  group_by(data_type, gene) %>%
  summarise(plt = build_plot(), .groups = "drop")

plts %>%
  # Create a panel per data type
  group_by(data_type) %>%
  group_walk(function(d, data_type) {
    # Arrange the plots in a grid
    # panel <- cowplot::plot_grid(
    #   plotlist = d$plt,
    #   align = 'v',
    #   ncol = 3,
    #   labels = paste0("HLA-", d$gene),
    #   greedy = F,
    #   label_size = 24
    # ) + theme(plot.margin = margin(r = 12, b = 12))
    dplts <- map2(d$plt, d$gene, ~.x + ggtitle(paste0("HLA-", .y)))
    panel <- (wrap_plots(dplts) +
                plot_layout(ncol=3) +
                plot_layout(guides = "collect") &
                theme(legend.position = 'bottom',
                      plot.title = element_text(hjust = 0.5, size = 24)))
    
    # Write (per data type) the panel to a separate figure
    ggsave(
      str_glue("results/figs/supplementary_clustering_{data_type}.pdf"),
      panel,
      width = 21,
      scale = 2,
      units = "cm"
    )
  })
