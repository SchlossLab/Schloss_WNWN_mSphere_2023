#!/usr/bin/env Rscript

library(tidyverse)
library(ggtext)

set.seed(19760620)

add_nl <- function(x) {

  paste0("N~L~ = ", format(as.numeric(x), big.mark = "\\\\,"))

}

pretty_distances <- c(bray = "Bray-Curtis",
                      euclidean = "Euclidean",
                      poisson = "Poisson",
                      logFC = "Mean Squared Difference",
                      uunifrac = "Unweighted UniFrac",
                      wunifrac = "Weighted UniFrac")

pretty_transform <- c(deseq = "DeSeq VS",
                    none = "Raw counts",
                    proportion = "Relative abundance",
                    subsample15 = "Sub-sampled (15%)",
                    rarefaction00 = "Rarefaction (0%)",
                    upperquartile = "Upper Quartile Log Fold Change")


plot_four <- function(cluster_method, transformation){

  read_tsv("data/simulation_clusters.tsv.gz") %>%
    rename(subset = fracCorrectPred, all = fracCorrect) %>%
    filter(simulation == "sim_a") %>%
    filter(filter == "filter" &
            transform %in% names(pretty_transform) &
            distance %in% names(pretty_distances)) %>%
    filter(method == cluster_method) %>%
    filter(fraction != "s1") %>%
    mutate(fraction = as.numeric(fraction)) %>%
    filter((distance == "bray" &
              transform %in% c("proportion", transformation, "none", "deseq")) |
          (distance == "euclidean" &
              transform %in% c("none", "deseq", transformation)) |
          (distance == "poisson" &
              transform %in% c(transformation, "none")) |
          (distance == "logFC" &
              transform %in% c("upperquartile", transformation)) |
          (distance == "uunifrac" &
              transform %in% c(transformation, "proportion")) |
          (distance == "wunifrac" &
              transform %in% c(transformation, "none", "proportion", "deseq"))
          ) %>%
    group_by(fraction, n_seqs, transform, distance) %>%
    summarize(mean = mean(all, na.rm = TRUE),
              lci = quantile(all, 0.025, na.rm = TRUE),
              uci = quantile(all, 0.975, na.rm = TRUE), .groups = "drop") %>%
    mutate(distance = factor(pretty_distances[distance],
                            levels = pretty_distances),
          transform = factor(pretty_transform[transform],
                              levels = pretty_transform)) %>%
    ggplot(aes(x = fraction, y = mean, group = transform,
              color = transform, shape = transform, fill = transform)) +
    geom_hline(yintercept = 0.5, color = "darkgray") +
    geom_line(position = position_dodge(width = 0.075)) +
    geom_linerange(aes(ymin = lci, ymax = uci),
                  alpha = 0.6, position = position_dodge(width = 0.075),
                  show.legend = FALSE) +
    geom_point(position = position_dodge(width = 0.075),
              size = 2) +
    facet_grid(n_seqs ~ distance,
              labeller = labeller(n_seqs = add_nl)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(
      values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#000000")
    ) +
    scale_color_manual(
      values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#000000")
    ) +
    scale_shape_manual(
      values = c(21, 22, 23, 24, 25, 21)
    ) +
    coord_cartesian(ylim = c(0.4, 1.05)) +
    labs(x = "Effect Size",
        y = "Accuracy",
        color = "Normalization Method:",
        fill = "Normalization Method:",
        shape = "Normalization Method:") +
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
                cluster_method, "_", transformation, "_fig_4.pdf"),
          width = 11, height = 10)

}

args <- commandArgs(trailingOnly = TRUE)

plot_four(args[1], args[2])
