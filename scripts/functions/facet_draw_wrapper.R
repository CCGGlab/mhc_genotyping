# Wrapper function for ggplot facet functions
# Used to prevent passing a million parameters
facet_draw_wrapper <- function(f, fun_type, parent, ...) {
  extra_pars <- list(...)
  
  function(self, ...) {
    facet_draw_params <-
      list(
        panels = c(
          "panels",
          "layout",
          "x_scales",
          "y_scales",
          "ranges",
          "coord",
          "data",
          "theme",
          "params"
        ),
        labels = c(
          "panels",
          "layout",
          "x_scales",
          "y_scales",
          "ranges",
          "coord",
          "data",
          "theme",
          "labels",
          "params"
        )
      )
    facet_params <- rlang::set_names(list(...), facet_draw_params[[fun_type]])
    args <- list(self, facet_params, ggplot2::ggproto_parent(parent, self))
    do.call(f, c(args, extra_pars))
  }
}