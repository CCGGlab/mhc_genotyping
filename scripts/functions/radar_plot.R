# coord_radar from package "ggigraphExtra"
# The radar coordinate system is a modification of polar coordinate system, commonly used for radar chart
coord_radar <- function() {
  ggproto(
    "CoordRadar",
    ggplot2::CoordPolar,
    theta = "x",
    r = "y",
    start = 0,
    direction = 1,
    is_linear = function(coord) TRUE,
    # Prevent clipping text of tool labels
    clip = "off",

    # CoordPolar ignores hjust and vjust parameters
    # Fix here
    render_fg = function(self, panel_params, theme) {
      # Call parent render function
      # res <- ggproto_parent(CoordPolar, self)$render_fg(panel_params, theme)

      theta <- ggplot2:::theta_rescale(self, panel_params$theta.major, panel_params)
      labels <- panel_params$theta.labels
      theta <- theta[!is.na(theta)]

      res <- grid::grobTree(
        if (length(labels) > 0) {
          element_render(
            theme, "axis.text.x",
            labels,
            unit(0.44 * sin(theta) + 0.5, "native"),
            unit(0.44 * cos(theta) + 0.5, "native"),
            # - Change justification of labels
            #   Justification set such that innermost point is aligned on outer circle
            # - Move labels closer to the drawing area
            #   (Drawing area until r = 0.4. (Hidden) outer circle at r = 0.45)
            hjust = (1 - sin(theta)) / 2, vjust = (1 - cos(theta)) / 2
            # hjust = 0.5, vjust = 0.5
          )
        },
        element_render(theme, "panel.border")
      )

      # Rotate text along the circle: convert to radians and rotate clockwise
      # res$children[[1]]$children[[1]]$rot = -((theta * (180 / pi)))
      res
    },
    render_bg = function(self, panel_params, theme) {
        # Call parent render function
        res <- ggproto_parent(CoordPolar, self)$render_bg(panel_params, theme)
        
        # These variables need to be recalculated here for r.minor.
        # Do not add 0.45 to rfine. We do not want an outer border.
        # Note: r_rescale is implemented uses scales::rescale
        # -> scales r.minor to the 0 → 0.4 range (drawing area)
        thetafine <- seq(0, 2 * pi, length.out = 100)
        rfine <- ggplot2:::r_rescale(self, panel_params$r.minor, panel_params$r.range)
        minorr <- paste("panel.grid.minor.", self$r, sep = "")
        
        grid_minor_y <- element_render(
            theme, minorr, name = "radius",
            x = rep(rfine, each = length(thetafine)) * sin(thetafine) + 0.5,
            y = rep(rfine, each = length(thetafine)) * cos(thetafine) + 0.5,
            id.lengths = rep(length(thetafine), length(rfine)),
            default.units = "native"
        )
        
        grid::addGrob(res, grid_minor_y)
    }
  )
}

# Customized grid structure
element_custom_grid <- function(...) {
  structure(ggplot2::element_line(...),
    class = c("custom_grid", "element_line", "element")
  )
}

# Render function for grid
# Remove the outer border
element_grob.custom_grid <- function(...) {
  params <- list(...)
  # Remove outer border
  outer_length <- params$id.lengths[[length(params$id.lengths)]]
  params$id.lengths[[length(params$id.lengths)]] <- 0

  params$x <- params$x[seq(1, length(params$x) - outer_length)]
  params$y <- params$y[seq(1, length(params$y) - outer_length)]

  do.call(ggplot2:::element_grob.element_line, params)
}

# Close the radar plot: connect first and last element.
StatConnect <- ggproto(
  "StatConnect",
  StatIdentity,
  compute_layer = function(data, params, layout) {
    data %>%
      group_by(group) %>%
      # Copy first row to the end of the dataframe
      group_modify(~ bind_rows(.x, .x[1, ]))
  }
)

radar_plot <-
  function(df,
           x_var = "tool",
           y_var = "decided_correct",
           color_var = "gene",
           labeller = identity) {

    # Sorting labels is obligatory.
    # Conversion character -> factor -> number for group index
    # Which dots are connected with polyline depends on this order.
    ggplot(
      arrange(df, .data[[x_var]], .data[[color_var]]),
      aes_string(
        x = x_var,
        y = y_var,
        colour = color_var,
        group = color_var
      )
    ) +
      coord_radar() +
      # To consider:
      # - 10: circle with cross
      # - 8: star
      # - 19: filled circle (default), e.g. with alpha = 0.7?
      # - 1: circle
      # Below: use unicode character of circle with cross (2A02)
      # We might consider using geom_text for better customization.
      # geom_text(label = "\u25D2", size=10, family = "Lucida Sans Unicode") +
      geom_point(shape = 10, size = 6, stroke = 1.5) +
      geom_path(size = 1.5, stat = StatConnect, na.rm = T) +
      scale_color_hue(label = as_mapper(labeller)) +
      scale_y_continuous(
        breaks = seq(0, 1, 0.1),
        minor_breaks = seq(0.05, 0.95, 0.1),
        expand = c(0, 0)
      ) +
      theme_bw() +
      theme(
        # Change font size
        axis.text.x = element_text(size = 24),
        # Hide labels for radial variable at the side
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        # Hide x / y axis titles. These are confusing for a radar plot.
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        # Tangential lines
        panel.grid.major.y = element_custom_grid(colour = "#0000003f", linetype = "solid"),
        panel.grid.minor.y = element_line(
          colour = "#0000003f", linetype =
            "dashed"
        ),
        # No radial lines
        panel.grid.major.x = element_blank(),
        # No border around panel
        panel.border = element_blank(),
        # Move legend to bottom
        legend.position = "bottom",
        legend.text = element_text(size = 18),
        legend.title = element_blank(), #element_text(size = 18),
        legend.key.width = unit(0, 'cm')
      ) +
      # Add label at grid lines
      geom_text(
        aes(x = 0, y = tick, label = lbl),
        data = tibble(
          tick = seq(0, 1, 0.1),
          lbl = seq(0, 100, 10)
        ),
        show.legend = F,
        colour = "black",
        family = "serif",
        size = 6,
        inherit.aes = FALSE
      )
  }
