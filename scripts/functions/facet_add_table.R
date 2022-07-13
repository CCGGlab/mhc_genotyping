library(tidyverse)
source("scripts/functions/facet_draw_wrapper.R")

assemble_table_facet <- function(plt, df) {
  plt %>%
    # Add y axis to appropriate column
    gtable::gtable_add_grob(
      df$y_axis,
      t = -1,
      l = df$y_axis_col,
      name = str_glue_data(df, "ytable-b-{COL}-{ROW}-{table_row}")
    ) %>%
    # Add panel to column right to the y axis
    gtable::gtable_add_grob(
      df$panel,
      t = -1,
      l = df$y_axis_col + 1,
      name = str_glue_data(df, "ptable-b-{COL}-{ROW}-{table_row}")
    )
}

append_table_axis <- function(plt, selected_tools, heights) {
  plt <- plt %>%
    # Add a separator
    gtable::gtable_add_rows(heights[[1]], pos = -1) %>%
    # Add a new row for the table
    gtable::gtable_add_rows(heights[[2]], pos = -1)
  
  # Add the objects that should be on that row
  selected_tools %>%
    mutate(y_axis_col = 1 + 4 * (COL - 1)) %>%
    group_by(COL) %>%
    group_split %>%
    reduce(assemble_table_facet, .init=plt)
}

assemble_plot_table <- function(self, facet_params, parent, tbl_df, tbl_plot, table_height) {
  parent_panel <- lift_dl(parent$draw_panels)(facet_params)
  row_separator <- map_chr(groups(tbl_df), rlang::as_string) %>%
    setdiff(names(facet_params$params$facets))
  
  tbl_df %>%
    # Separate the data for the tables per panel
    inner_join(facet_params$layout, by = names(facet_params$params$facets)) %>%
    mutate(across(PANEL, as.numeric)) %>%
    group_by(across(all_of(row_separator)), PANEL, ROW, COL) %>%
    summarise(plt = list(
      tbl_plot(cur_data(), cur_group(), facet_params$x_scales)
    ), .groups = "drop") %>%
    # unnest_wider(plt) %>%
    
    # For some bizarre reason the line above does not work anymore
    # after a conda update
    # Bug?
    # Workaround:
    # - First encapsulate the plots in a list
    # - Unnest the list to columns
    # - Remove the encapsulation
    mutate(plt = map_depth(plt, 2, list)) %>%
    unnest_wider(plt) %>%
    mutate(across(c(y_axis, panel), flatten)) %>%
    
    # Add the table rows (one per MHC class)
    group_by(across(all_of(row_separator))) %>%
    mutate(table_row=cur_group_id()) %>%
    group_split %>%
    reduce(append_table_axis, table_height, .init = parent_panel) %>%
    # Increase the size of the y-axis labels
    # -> Set widths of columns 1 and 5 to 3 cm
    modify_in('widths', modify_at, c(1,5), ~unit(3.5, 'cm'))
}

construct_label <- function(panel, lbl, posy, posx) {
  panel %>%
    gtable::gtable_add_grob(grid::textGrob(
      lbl,
      x = unit(0.5, 'npc') + unit('0.1', 'cm'),
      y = unit(0.5, 'npc'),
      rot = 90,
      gp = grid::gpar(fontsize = 18)
    ),
    t = posy,
    l = posx)
}

draw_labels <- function(self, facet_params, parent, table_labels) {
  # Draw axis labels of actual drawing area
  res <- lift_dl(parent$draw_labels)(facet_params)
  
  # Calculate the label positions
  positions <- nest_by(res$layout, name) %>% deframe
  table_ycol <- positions[["ylab-l"]]$l
  table_yrow <- positions[["ytable-b-1-1-1"]]$t
  x_axis_row <- positions[["axis-b-1-1"]]$t
  x_label_height <- res$heights[[positions$`xlab-b`$t]]
  
  # Add labels to the bottom table
  reduce2(table_labels,
          table_yrow + 2 * (seq_along(table_labels) - 1),
          construct_label,
          table_ycol,
          .init = res) %>%
    # Move x label under the x-axis
    gtable::gtable_add_rows(x_label_height, pos = x_axis_row) %>%
    modify_in("layout", function(x) {
      x[x$name == "xlab-b", c("t", "b")] <- x_axis_row + 1
      x
    })
}

facet_add_table <- function(parent, tbl_df, tbl_plot, table_labels, table_height) {
  ggproto(
    NULL,
    parent,
    draw_panels = facet_draw_wrapper(assemble_plot_table, "panels", parent, tbl_df, tbl_plot, table_height),
    draw_labels = facet_draw_wrapper(draw_labels, "labels", parent, table_labels)
  )
}