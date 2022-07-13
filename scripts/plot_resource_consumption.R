library(tidyverse)
metrics <- readRDS("data/resource_consumption.rds") %>%
  mutate(tool = recode(tool, "HLA-LA" = "HLA*LA", "HLAScan" = "HLAscan"))

# Colour palette for 13 tools
colour_palette13 <-
  c(
    "#d2b798",
    "#815ace",
    "#6ccb5e",
    "#c14ca4",
    "#c6ce50",
    "#cc943c",
    "#8dbad7",
    "#cc5334",
    "#7cd4ae",
    "#ca4766",
    "#61803d",
    "#537371",
    "#8b5c42"
  )

resource_plot <- function(property, label) {
  property <- substitute(property)
  ggplot(metrics, aes(
    x = reorder(fct_cross(data_type, tool, sep = '\31'),
                !!property, FUN = median, na.rm=T),
    y = !!property,
    colour = tool
  )) +
    geom_boxplot(outlier.shape=NA, lwd=1) +
    ggpubr::theme_pubclean() +
    scale_y_log10(label,
                  # Almost perfect equidistant:
                  # (rounded to 1 significant digit)
                  breaks = ~signif(10^seq(log10(.[[1]]), log10(.[[2]]), length.out = 10), 1),
                  # n.breaks = 20,
                  labels = prettyNum) +
    # scale_y_continuous(label, trans="log10", n.breaks = 10, labels = prettyNum) +
    scale_x_discrete('tool', labels = ~map(str_split(., pattern = '\31'), 2)) +
    scale_colour_manual(values = colour_palette13) +
    #scale_colour_hue(guide = guide_none()) +
    theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
          axis.title.x = element_blank()) +
    facet_wrap(~data_type, scales = "free_x")
}

time_plt <- resource_plot(duration, 'time (hours)')
# cpu_plt <- resource_plot(cpu_mean, 'average CPU usage (%)')
mem_plt <- resource_plot(mem_gib_max, 'maximum memory usage (GiB)')

plt_panel <- ((time_plt / mem_plt) +
                plot_layout(guides = "collect") +
                plot_annotation(tag_levels = 'A')) &
  theme(
    legend.position = "bottom",
    legend.box.margin = margin(0, 0, 0, 12),
    plot.tag = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 18),
    # Font size of the panel and axis labels
    strip.text = element_text(size = 18),
    axis.text = element_text(size = 18),
    text = element_text(size = 18)
  )

ggsave(
  str_glue("results/figs/resources_panel.pdf"),
  plt_panel,
  width = 13.86,
  height = 20,
  scale = 2,
  units = "cm"
)
