library(tidyverse)

stat_mark_final <- ggproto(
  NULL,
  StatIdentity,
  compute_layer = function(self, data, scales, ...) {
    data %>%
      group_by(colour,PANEL) %>%
      mutate(final=(y >= max(y))) %>%
      group_by(colour, PANEL, final) %>%
      mutate(shape=final & x == min(x)) %>%
      ungroup %>%
      select(-final)
  }
)