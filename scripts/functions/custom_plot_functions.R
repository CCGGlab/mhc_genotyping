# Custom transformations
reverse_log_trans <- function (base = exp(1)) 
{
  force(base)
  trans <- function(x) {x = pmax(x, 1e-10); -log(x, base)}
  inv <- function(x) (base^(-x))
  scales::trans_new(paste0("rlog-", format(base)), trans, inv, scales::log_breaks(base = base), 
                    domain = c(1e-100, Inf))
}

size_labels <- function(labels) {
  scientific_notation <- format(labels, scientific = T)
  mantissa_exponent <- str_split(scientific_notation, "e")
  mantissa <- as.numeric(map_chr(mantissa_exponent, 1))
  exponents <- as.numeric(map_chr(mantissa_exponent, 2))
  
  prefix_str <- if_else(mantissa == "1", "", paste0(mantissa, " *"))
  labels_str <- str_glue("{prefix_str}10^{{{exponents}}}")
  labels_str[[1]] = paste0("\\leq ", labels_str[[1]])
  latex2exp::TeX(str_glue("${labels_str}"))
}

# Custom theme elements
## Facet strip ("category title") with only a line at the bottom
element_bottom_line <- function() {
    structure(element_rect(),
        class = c("line_strip", "element_rect", "element")
    )
}

element_grob.line_strip <- function(element,
                                    x = 0.5,
                                    y = 0.5,
                                    width = 1,
                                    height = 1,
                                    fill = NULL,
                                    colour = NULL,
                                    size = NULL,
                                    linetype = NULL,
                                    ...) {
  element_gp <-
    grid::gpar(
      lwd = NULL,
      col = colour,
      fill = fill,
      lty = linetype
    )

  grid::linesGrob(unit(c(0, 1), "npc"),
    unit(c(0, 0), "npc"),
    gp = element_gp
  )
}

geom_custom_legend <- function(custom_aes,
                                 ...) {
  # Create a blank layer with custom key_glyph function
  # and customized version of GeomBlank to set the required_aes
  ggproto(
    NULL,
    geom_blank(...),
    geom = ggproto(
      "GeomCustomLegend",
      ggplot2::GeomBlank,
      # Define which aesthetics should be included in the legend
      # → Scales with aesthetic "ycolour" will be displayed
      default_aes = ggplot2::GeomTile$default_aes,
      required_aes = custom_aes,
      draw_key = function(data, params, size) {
        if (custom_aes %in% names(data)) {
          data$fill <- data[[custom_aes]]
          ggplot2::draw_key_polygon(data, params, size)
        } else {
          # Ignore if the scale used for rendering does not contain the aes
          zeroGrob()
        }
      }
    ),
    show.legend = TRUE,
    inherit.aes = TRUE
  )
}

facet_custom = function(...) {
  ggproto(
    NULL,
    facet_wrap(...),
    # Change axis labels based on calculated aesthetics
    draw_panels = function(self, panels, layout, x_scales, y_scales, panel_params, coord, data, theme, params) {
      # Build y label here based on current value in data$ycolour (is already transformed by scale)
      # Guides are contained in the "ranges" parameter.
      mapper <- distinct(data[[1]], y, ycolour, ylabel)
      
      for (i in seq_along(panel_params)) {
        panel_params[[i]]$guides$y$key = panel_params[[i]]$guides$y$key %>%
          rownames_to_column("idx") %>%
          inner_join(mutate(mapper, idx = as.character(y), y = NULL), by = "idx") %>%
          transmute(.value, .label = str_glue('<span style="color: {ycolour}">{ylabel}</span>'), y, x)
      }
      
      ggproto_parent(ggplot2::FacetWrap, self)$draw_panels(panels, layout, x_scales, y_scales, panel_params, coord, data, theme, params)
    }
  )
}