#!/usr/bin/env Rscript

library(tidyverse)
library(ggtext)

set.seed(19760620)

add_es <- function(x) {

  paste0("ES = ", format(as.numeric(x), nsmall = 2L))

}

add_nl <- function(x) {

  paste0("N~L~ = ", format(as.numeric(x), big.mark = "\\,"))

}

pretty_distances <- c(bray = "Bray-Curtis",
                      euclidean = "Euclidean",
                      poisson = "Poisson",
                      logFC = "Mean Squared Difference",
                      uunifrac = "Unweighted UniFrac",
                      wunifrac = "Weighted UniFrac")

plot_five <- function(sim, cluster_method, transformation){

  clusters <- read_tsv("data/simulation_clusters.tsv.gz") %>%
    filter(simulation == paste0("sim_", sim)) %>%
    filter(filter == "filter") %>%
    filter(method == cluster_method) %>%
    filter(distance != "bcv") %>%
    filter(str_detect(transform, transformation)) %>%
    filter(fraction != "s1" & fraction != "1") %>%
    select(fraction, n_seqs, rep, transform, distance, method, fracCorrect) %>%
    group_by(fraction, n_seqs, transform, distance, method) %>%
    summarize(median = median(fracCorrect),
              lci = quantile(fracCorrect, 0.025),
              uci = quantile(fracCorrect, 0.975), .groups = "drop") %>%
    mutate(quant = as.numeric(str_replace(transform, transformation, "")),
          distance = pretty_distances[distance])

  clusters %>%
    ggplot(aes(x = quant, y = median, group = distance,
              color = distance, fill = distance, shape = distance)) +
    geom_segment(x = 0, y = 1,
                xend = 100, yend = 0, color = "gray") +
    geom_linerange(aes(ymin = lci, ymax = uci),
                  alpha = 0.6, position = position_dodge(width = 1.5),
                  show.legend = FALSE) +
    geom_line(position = position_dodge(width = 1.5)) +
    geom_point(position = position_dodge(width = 1.5),
              size = 2) +
    facet_grid(n_seqs ~ fraction,
              labeller = labeller(fraction = add_es, n_seqs = add_nl)) +
    scale_fill_manual(
      values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#000000")
    ) +
    scale_color_manual(
      values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#000000")
    ) +
    scale_shape_manual(
      values = c(21, 22, 23, 24, 25, 21)
    ) +
    coord_cartesian(ylim = c(0.4, 1.05),
                    xlim = c(0, 40)) +
    labs(x = "Library Size Minimum Quantile (%)",
        y = "Accuracy",
        color = "Distance Calculation:",
        fill = "Distance Calculation:",
        shape = "Distance Calculation:") +
    theme_light() +
    theme(
      legend.position = "top",
      legend.box = "horizontal",
      legend.margin = margin(t = 0, r = -0.5, b = 0, l = -0.5, unit = "cm"),
      strip.text.y = element_markdown()
    ) +
    guides(color = guide_legend(nrow = 1),
          fill = guide_legend(nrow = 1),
          shape = guide_legend(nrow = 1))

  ggsave(paste0("results/figures/",
              "fig_5_", cluster_method, "_", transformation, "_", sim, ".pdf"),
          width = 11, height = 10)

}

args <- commandArgs(trailingOnly = TRUE)

plot_five(args[1], args[2], args[3])
