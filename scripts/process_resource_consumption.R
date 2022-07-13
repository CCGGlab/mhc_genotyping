library(tidyverse)
library(patchwork)

# future::plan(future::multicore)
# options(future.rng.onMisuse="ignore")

units <- c('KiB', 'MiB', 'GiB') %>%
  imap_dfr(~tibble(.x,  2^(10 * .y))) %>%
  deframe()

predictions <- readRDS("temp/hla_predictions_ggroup_resmon.rds")
mapped_tools <- c(
  "HLAminer" = "HLAminerHPRA",
  "HLAScan" = "HLAscan",
  "Kourami" = "kourami",
  "HLA-VBSeq" = "HLAVBseq"
)

# List all files in resource monitoring folder
resource_folders <- Sys.glob("temp/resmon_single_tcga3/*/*/*")

resource_folders %>%
  str_split('/') %>%
  map(~ .[3:5]) %>%
  transpose %>%
  set_names(c("data_type", "tool", "sample_id")) %>%
  as_tibble %>%
  mutate_all(simplify) %>%
  group_by(data_type, tool) %>%
  count

metrics <- resource_folders %>%
  purrr::map_dfr(function (f) {
    start <- str_glue("{f}/start.txt")
    end <- str_glue("{f}/end.txt")

    path_parts <- str_split(f, '/', simplify = T)[, 3:5]

    output_exists <- path_parts %>%
      set_names(c("data_type", "tool", "sample_id")) %>%
      str_glue_data(paste0(
        "temp/output_resmon_single_tcga3/",
        "{data_type}/{tool}/*{sample_id}*"
      )) %>%
      Sys.glob() %>%
      {
        length(.) > 0
      }

    predictions_made <- 0
    if (output_exists) {
      data_type <- path_parts[[1]]
      tool <- path_parts[[2]]

      if (tool %in% names(mapped_tools)) {
        tool <- mapped_tools[[tool]]
      }

      predictions_made <- predictions %>%
        chuck("single_tcga", data_type, tool) %>%
        filter(str_starts(path_parts[[3]], sample_id)) %>%
        filter(!all(is.na(c_across(!sample_id)))) %>%
        nrow
    }

    if ((predictions_made == 0) | !file.exists(start) | !file.exists(end)) {
      print("WARNING: no prediction made for:")
      print(f)
      # Sample not completed successfully
      return(tibble(path = f))
    }

    duration <-
      (as.numeric(read_lines(end)) - as.numeric(read_lines(start))) / 3600

    p <- str_glue("{f}/resources.txt") %>%
      read_lines %>%
      str_remove('^.{7}') %>%
      map(possibly(jsonlite::fromJSON, otherwise = NULL))

    mem_human <- p %>%
      map("MemUsage") %>%
      str_remove(' */ *.*$') %>%
      str_match("([0-9.]*)([KMG]iB)")
    mem_bytes <-
      as.numeric(mem_human[, 2]) * unname(units[mem_human[, 3]])
    mem_gib_max <- max(mem_bytes, na.rm = T) / (2 ^ 30)

    cpu <- p %>%
      map("CPUPerc") %>%
      str_remove('%$') %>%
      as.numeric()
    cpu_mean <- mean(cpu, na.rm=T)
    cpu_max <- max(cpu, na.rm=T)

    res <- tibble(path = f, duration, mem_gib_max, cpu_mean, cpu_max)
    res
  }) %>%
  mutate(across(path, str_remove, 'temp/resmon_single_tcga3/')) %>%
  tidyr::separate(path, c("data_type", "tool", "sample_id"), '/')

saveRDS(metrics, "data/resource_consumption.rds")