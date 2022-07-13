dendrogram_coord <- function(column_clusters, row_clusters) {
  CoordDendro <- ggproto(
    NULL,
    coord_cartesian(),
    render_axis_h = function(self, panel_params, theme) {
      # Calculate original horizontal axis
      res <- ggproto_parent(CoordCartesian, self)$render_axis_h(panel_params, theme)
      
      # Calculate dendrogram based on column clusters
      dendro <- ggdendro::ggdendrogram(column_clusters) %>%
        cowplot::get_panel()
      
      dendro_height <- unit(48, 'pt')
      gtdendro <- gtable::gtable(
        widths = unit(1, 'npc'),
        heights = dendro_height,
        # Calculate viewport
        # Align the vp-bottom to the bottom of the axis drawing area
        vp = grid::viewport(
          y = unit(0, 'npc'),
          just = "bottom",
          height = dendro_height
        )
      ) %>%
        gtable::gtable_add_grob(dendro, t=1,l=1,b=1,r=1)
      
      res$top <- gtdendro
      res
    },
    render_axis_v = function(self, panel_params, theme) {
      # Calculate original horizontal axis
      res <- ggproto_parent(CoordCartesian, self)$render_axis_v(panel_params, theme)
      
      # Calculate dendrogram based on row clusters
      # dendro <-
      #   ggdendro::ggdendrogram(
      #     row_clusters,
      #     labels = F,
      #     leaf_labels = F,
      #     rotate = T
      #   ) %>%
      #   cowplot::get_panel()
      ddata <- ggdendro::dendro_data(row_clusters, type = "rectangle")
      p <- ggplot(ggdendro::segment(ddata)) + 
        geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
        coord_flip() + 
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        ggdendro::theme_dendro()
      dendro <- cowplot::get_panel(p)
      
      # TODO: find other solution
      # segments run from ca. y = 0.05 to 0.95
      # -> 5% added to use entire y-axis...
      dendro_width<- unit(48, 'pt')
      gtdendro <- gtable::gtable(
        # heights = unit(1.05, 'npc'),
        heights = unit(1, 'npc'),
        widths = dendro_width,
        # Calculate viewport
        # Align the vp-left to the left of the axis drawing area
        vp = grid::viewport(
          x = unit(0, 'npc'),
          just = "left",
          width = dendro_width
        )
      ) %>%
        gtable::gtable_add_grob(dendro, t=1,l=1,b=1,r=1)

      res$right <- gtdendro
      res
    }
  )
}