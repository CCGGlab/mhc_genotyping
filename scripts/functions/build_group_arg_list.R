# Convert a grouped data frame to a list with one element per group
# A common use case is to apply pmap on this output: pmap(output, f, extra_args)
# -> f(df_arg = df, group_arg = key, extra_args) will by applied on all groups
# This function is also compatible with furrr::future_pmap
# -> Allows parallellization of the calculations per group
build_group_arg_list <- function(df, df_arg, group_arg) {
  list(
    group_split(df, .keep = T),
    group_keys(df) %>% rowwise %>% group_split
  ) %>%
    set_names(c(df_arg, group_arg))
}
