library(tidyverse)
extract_dataset <- function(dataset_name, all_predictions) {
  all_predictions[[dataset_name]] %>%
    enframe(name = "data_type") %>%
    mutate(value = map(value, enframe, name = "tool")) %>%
    unnest(value) %>%
    unnest(value)
}