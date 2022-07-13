library(tidyverse)
library(patchwork)
future::plan(future::multicore)
#fo <- furrr::furrr_options(seed=TRUE)
options(future.rng.onMisuse="ignore")

source("scripts/functions/extract_dataset.R")
source("scripts/functions/evaluate_predictions.R")

predictions_all <- readRDS("data/hla_predictions_ggroup.rds")
meta_predictions_compact_all <- readRDS("data/hla_predictions_ggroup_meta.rds")$compact

for (dataset in c("1kg", "TCGA")) {
  pred_ds <- extract_dataset(dataset, predictions_all)
  predictions_meta <- extract_dataset(dataset, meta_predictions_compact_all) %>%
    filter(tool == "meta") %>%
    mutate(tool = "Meta")
  pred_ds <- pred_ds %>%
    bind_rows(predictions_meta)
  
  # Use exact tool names and case styling used by authors
  formatted_tool_name = c(HLAminerHPRA = "HLAminer",
                          HLAVBseq = "HLA-VBSeq",
                          kourami = "Kourami",
                          `HLA-LA` = "HLA*LA")
  
  evaluate_combinations <- function(d) {
    # Create list of (tool, data) elements
    # - Group by tool
    # - Separate all tools to separate element of list
    # - Extract the one and only data frame from the data column
    g <- nest_by(d, tool) %>%
      group_split %>%
      map(map_at, "data", pluck, 1)
    
    # Cross product of all tools
    cross_df(list(x = g, y = g)) %>%
      # filter(map_chr(x, "tool") != map_chr(y, "tool")) %>%
      # Evaluate the predictions by tool "x" using "y" as gold standard
      furrr::future_pmap_dfr(function(x, y) {
        # Exact genotype matching
        # Sample corresponds if both alleles equal
        a <- list(x$data, y$data) %>%
          map(function(.) {
            pivot_longer(., !sample_id) %>%
              separate(name, c("gene", "idx")) %>%
              group_by(sample_id, gene) %>%
              arrange(value, .by_group = T) %>%
              mutate(idx = seq_along(value)) %>%
              # Drop sample_id, gene pairs where the prediction is NA
              # for any of the 2 alleles
              group_by(sample_id, gene) %>%
              filter(!any(is.na(c_across(value))))
          })
        
        r <- inner_join(a[[1]], a[[2]], by=  c("sample_id", "gene", "idx")) %>%
          group_by(sample_id, gene) %>%
          mutate(is_correct = all(value.x == value.y)) %>%
          group_by(gene) %>%
          summarise(accuracy = mean(is_correct)) %>%
          bind_cols(tool.x = x$tool, tool.y = y$tool)
        r
      })
  }
  
  metrics <- pred_ds %>%
    mutate(tool = if_else(tool %in% names(formatted_tool_name), unname(formatted_tool_name[tool]), tool)) %>%
    group_by(data_type) %>%
    summarise(evaluate_combinations(cur_data()), .groups="drop")

  saveRDS(metrics, str_glue("data/concordance_tools_{dataset}_exact.rds"))
  # metrics <- readRDS(str_glue("data/concordance_tools_{dataset}_exact.rds"))
  
  # Heatmap of concordance when both tools were able to make a prediction
  plot_conco <- function(df) {
    m = df %>%
      select(tool.x, tool.y, accuracy) %>%
      spread(tool.y, accuracy) %>%
      mutate(across(!tool.x, as.numeric)) %>%
      mutate(across(!tool.x, ~ replace_na(., 0))) %>%
      column_to_rownames("tool.x") %>%
      as.matrix()
    
    clusters_x <-
      hclust(dist(m, method = "euclidean"), method = "ward.D")
    clusters_y <-
      hclust(dist(t(m), method = "euclidean"), method = "ward.D")
    
    ggplot(df, aes(x=tool.x, y=tool.y, fill=accuracy)) +
      geom_raster() +
      theme_minimal() +
      theme(strip.text = element_text(size=14),
            axis.text.x = element_text(angle=45, hjust=1, vjust=1),
            axis.title.x = element_blank(),
            axis.title.y = element_blank()) +
      colorspace::scale_fill_continuous_sequential(limits=c(0,1),
                                                   palette = "OrYel",
                                                   na.value="grey50",
                                                   rev = T,
                                                   trans = scales::exp_trans(base = 5),
                                                   breaks = c(0.25, 0.5, 0.8, 0.9, 1),
                                                   label = function(x) {sprintf("%0.0f%%", 100*x)}) +
      scale_x_discrete(limits = rownames(m)[clusters_x$order]) +
      scale_y_discrete(limits = rownames(m)[clusters_x$order]) +
      #scale_y_discrete(limits = colnames(m)[clusters_y$order]) +
      labs(fill="concordance")
  }
  
  plots <- metrics %>%
    drop_na %>%
    group_by(data_type, gene) %>%
    summarise(plt=list(plot_conco(cur_data()) ))
  
  # + theme(plot.margin=margin(t=12,l=72))
  
  group_by(plots, data_type) %>%
    group_walk(function(d, data_type) {
      # panel <- cowplot::plot_grid(
      #   plotlist = d$plt,
      #   align = 'v',
      #   ncol = 3,
      #   #labels = d$gene,
      #   labels = "AUTO",
      #   label_x = 0.05, label_y = 0.95,
      #   hjust = 0, vjust = 1,
      #   greedy = F,
      #   label_size = 24
      # ) + theme(plot.margin = margin(r=12, b=12))
      
      dplts <- map2(d$plt, d$gene, ~.x + ggtitle(paste0("HLA-", .y)))
      panel <- (wrap_plots(dplts) +
                  plot_layout(ncol=3) +
                  plot_layout(guides = "collect") &
                  theme(legend.position = 'bottom',
                        legend.key.width = unit(2, 'cm'),
                        legend.title = element_text(hjust = 0.5, vjust = 0.5),
                        plot.title = element_text(hjust = 0.5, size = 24)))
      
      
      ggsave(
        str_glue("results/figs/supplementary_concordance_{dataset}_{data_type}_exact.pdf"),
        panel,
        width = 21,
        height = 15,
        scale = 2,
        units = "cm"
      )
    })
}