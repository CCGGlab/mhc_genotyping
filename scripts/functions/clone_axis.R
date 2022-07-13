clone_axis <- function(x_scale_main) {
  # - Do not re-train the range of the scale
  # - Set the range equal to the range of the main panel
  x_range_main <- x_scale_main$range$range # the range of this scale
  x_scale_table <-
    ggproto(
      NULL,
      x_scale_main,
      train = function(self, x) {
        # When the train function is called always set the range equal
        # to the range of the main panel x-axis
        self$range$range <- x_range_main
      }
    )
}