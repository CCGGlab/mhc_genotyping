library(tidyverse)
source("scripts/functions/clone_axis.R")
source("scripts/functions/facet_add_table.R")
source("scripts/functions/stat_mark_final.R")

# Functions
# Mapping between tool names as used in the data frame to the names on the plot
formatted_tool_name = c(HLAminerHPRA = "HLAminer",
                        HLAVBseq = "HLA-VBSeq",
                        kourami = "Kourami",
                        `HLA-LA` = "HLA*LA")

reformat_tool <- function(tool) {
  if_else(tool %in% names(formatted_tool_name), unname(formatted_tool_name[tool]), tool)
}

build_plt_dfs <- function(results_all_levels) {
  plt_df <- results_all_levels %>%
    # Reformat the tool names
    mutate(tools = modify(tools, reformat_tool)) %>%
    # Calculate the number of selected tools
    mutate(num_tools = map_int(tools, length))
  
  # Add the MHC class if it is missing
  if ("mhc_class" %in% names(results_all_levels)) {
    selected_tools_df <- plt_df %>%
      # Extract list of (unique) selected tools
      distinct(data_type, mhc_class, tools, num_tools) %>%
      # Determine the earliest n where a tool was selected
      unnest(tools) %>%
      group_by(across(!num_tools)) %>%
      mutate(min(num_tools)) %>%
      # How to separate the table over facets
      group_by(data_type, mhc_class)
  } else {
    plt_df <- plt_df %>%
      mutate(mhc_class = if_else(gene %in% c("A", "B", "C"), "I", "II"))
    selected_tools_df <- NULL
  }
  
  list(plt_df = plt_df, selected_tools_df = selected_tools_df)
}

# Plot the table below the accuracy plots
plot_selected_tools <- function(df, key, x_scales) {
  # We want to keep the same scale for the x-axis as the main panel:
  # x-scale of the main panel
  x_scale_table <- clone_axis(x_scales[[key$PANEL]])
  
  # Plot a table with selected tools below the main panel
  tbl_plt <- ggplot(df, aes(x = num_tools,
                            # order y-axis in increasing order of num_tools
                            y = reorder(tools,-`min(num_tools)`))) +
    theme_minimal() +
    # Only display major grid lines
    theme(panel.grid = element_blank(),
          panel.grid.major.y = element_line(colour = "grey70"),
          axis.text.y = element_text(size = 16)) +
    # Cloned x-axis of main panel
    x_scale_table +
    # Indicate selected tools with black circles of size 4
    geom_point(shape = 16, size = 4) +
    # Set colour of dots to black and hide legend
    scale_fill_manual(values = "black",
                      guide = guide_none())
  
  # Return result: y-axis and panel
  list(y_axis = cowplot::get_y_axis(tbl_plt),
       panel = cowplot::get_panel(tbl_plt))
}

# Plot the main panel (the accuaracy of the metaclassifier)
plt_compact_metaclassifier <- function(plt_df, selected_tools_df) {
  # With or without a table below the main plot:
  # split into a panel per data type
  applied_layout <- facet_wrap(~ data_type, scales = "free_x")
  
  if (!is.null(selected_tools_df)) {
    # Add table below the plot if selected_tools_df is provided
    applied_layout <- facet_add_table(applied_layout,
                                      selected_tools_df,
                                      plot_selected_tools,
                                      c("MHC-I", "MHC-II"),
                                      # margin, table height
                                      list(unit(0.5, 'cm'), unit(6, 'cm')))
    # Change aspect ratio and margin depending on whether we put a table below or not
    aspect_ratio_theme <- theme(# aspect.ratio = 1,
      # Use less margin to left because of MHC-I, MHC-II and accuracy labels
      plot.margin = margin(l = 0.5, r = 1.5, t = 1, b = 1, unit = 'cm'))
      # plot.margin = margin(l = 5.7, r = 5.7, unit = 'cm'))
  } else {
    aspect_ratio_theme <- theme(#aspect.ratio = 0.8,
                                plot.margin = margin(l = 1, r = 1, t = 1, b = 1, unit = 'cm'))
  }
  
  # Main plot:
  ggplot(plt_df, aes(x = num_tools, y = S, colour = gene)) +
    theme_minimal() +
    theme(
      # Put the legend on top instead of right to the plot
      legend.position = "top",
      # legend.box = "vertical",
      # legend.box.just = "left",
      legend.text = element_text(size = 18),
      # Font size of the panel labels (data type)
      strip.text = element_text(size = 18),
      # Hide minor grid lines
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 18),
      text = element_text(size = 18)
    ) +
    # Add a table below the plot if requested
    # (and adjust aspect ratio and margin accordingly)
    applied_layout +
    aspect_ratio_theme +
    scale_x_continuous(
      name = "number of tools",
      # Add a margin of 0.2 to both sides of the x-axis
      # Place the ticks at all integers
      limits = lift_dl(as_mapper( ~ c(.x - 0.2, .y + 0.2))),
      breaks = lift_dl(as_mapper( ~ seq(.x + 0.2, .y - 0.2))),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      "accuracy",
      # Display accuracy on y-axis as percentage
      label = function(x) {
        sprintf("%0.0f%%", 100 * x)
      }
    ) +
    # Plot accuracy per gene
    # Use an open circle for first final point per gene
    geom_point(size = 4, stat = stat_mark_final) +
    # Shape: closed circle and open circle
    scale_shape_manual(values = c(19, 21)) +
    # Also connect the dots
    geom_line() +
    # Add HLA as a prefix in the colour legend
    scale_colour_hue(
      label = function(x) {
        paste0("HLA-", x)
      },
      guide = guide_legend(nrow = 2, order = 1)
    ) +
    # Add MHC class average as black lines with different dash style
    geom_point(aes(group = mhc_class), stat = "summary", fun = "mean") +
    geom_line(
      aes(linetype = mhc_class),
      stat = "summary",
      fun = "mean",
      colour = "black",
      size = 1
    ) +
    # Different dash style for MHC-I and MHC-II
    scale_linetype("MHC class",
                   limits = c("I", "II"),
                   guide = guide_legend(nrow=1, order = 2))
}

# Read the results of the optimization and plots
res <- readRDS("data/meta_optimized.rds") %>%
  map(build_plt_dfs) %>%
  map(lift_dl(plt_compact_metaclassifier))

# Plot for compact per class metaclassifier as main figure
ggsave(
  str_glue("results/figs/metaclassifier.pdf"),
  res$class,
  width = 13.86, # 15 cm for IJMS: 21 - 2*(1.27) - 4.6
  height = 15, # Good looking aspect ratio
  scale = 2,
  units = "cm"
)

# ggsave(
#   str_glue("results/figs/metaclassifier_tcga.pdf"),
#   res$class_tcga,
#   width = 21,
#   height = 15,
#   scale = 2,
#   units = "cm"
# )
# 
# # Plot for per gene classifier as supplementary figure
# ggsave(
#   str_glue("results/figs/metaclassifier_per_gene.pdf"),
#   res$gene,
#   width = 21,
#   height = 10,
#   scale = 2,
#   units = "cm"
# )
