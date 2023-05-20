#!/usr/bin/env Rscript

library(tidyverse)
library(ggtext)

set.seed(19760620)

add_nl <- function(x) {

  paste0("&Ntilde;~L~ = ", format(as.numeric(x), big.mark = "\\\\,"))

}

pretty_distances <- c(bray = "Bray-Curtis",
                      euclidean = "Euclidean",
                      poisson = "Poisson",
                      logFC = "Mean Squared Difference",
                      uunifrac = "Unweighted UniFrac",
                      wunifrac = "Weighted UniFrac")

pretty_transform <- c(deseq = "DESeq VS",
                    none = "Raw counts",
                    proportion = "Relative abundance",
                    subsample15 = "Sub-sampled (15%)",
                    rarefaction00 = "Rarefaction (0%)",
                    upperquartile = "Upper Quartile Log Fold Change")


plot_four <- function(gp_log, cluster_method, transformation, deseq) {

  df <- readr::read_tsv("data/simulation_clusters.tsv.gz") %>%
    dplyr::rename(subset = fracCorrectPred,
                  all = fracCorrect) %>%
    dplyr::filter(model == gp_log & distribution == "random") %>%
    dplyr::filter(filter == "filter" &
            transform %in% names(pretty_transform) &
            distance %in% names(pretty_distances)) %>%
    dplyr::filter(method == cluster_method) %>%
    dplyr::mutate(fraction = as.numeric(fraction)) %>%
    dplyr::filter((distance == "bray" &
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
    dplyr::group_by(fraction, n_seqs, transform, distance) %>%
    dplyr:: summarize(median = median(all, na.rm = TRUE),
                      lci = quantile(all, 0.025, na.rm = TRUE),
                      uci = quantile(all, 0.975, na.rm = TRUE),
                      .groups = "drop")

  if (deseq == "nds") {

    df <- df %>%
      filter(transform != "deseq" |
              (transform == "deseq" & distance == "euclidean"))
  }

  df %>%
    dplyr::mutate(distance = factor(pretty_distances[distance],
                            levels = pretty_distances),
          transform = factor(pretty_transform[transform],
                              levels = pretty_transform)) %>%
    ggplot2::ggplot(ggplot2::aes(x = fraction, y = median,
                                group = transform, color = transform,
                                shape = transform, fill = transform)) +
    ggplot2::geom_line(position = position_dodge(width = 0.075)) +
    ggplot2::geom_linerange(aes(ymin = lci, ymax = uci),
                  alpha = 0.6, position = position_dodge(width = 0.075),
                  show.legend = FALSE) +
    ggplot2::geom_point(position = position_dodge(width = 0.075),
                        size = 2) +
    ggplot2::facet_grid(n_seqs ~ distance,
              labeller = ggplot2::labeller(n_seqs = add_nl)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::scale_fill_manual(
      values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")
    ) +
    ggplot2::scale_color_manual(
      values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")
    ) +
    ggplot2::scale_shape_manual(
      values = c(21, 22, 23, 24, 25)
    ) +
    ggplot2::coord_cartesian(ylim = c(0.4, 1.05)) +
    ggplot2::labs(x = "Effect Size",
        y = "Accuracy",
        color = NULL, #"Normalization Method:",
        fill = NULL, #"Normalization Method:",
        shape = NULL #"Normalization Method:"
        ) +
    ggplot2::theme_light() +
    ggplot2::theme(
      legend.position = "top",
      legend.box = "horizontal",
      legend.margin = ggplot2::margin(t = 0, r = -0.5, b = 0, l = -0.5,
                                      unit = "cm"),
      strip.text.y = ggtext::element_markdown()
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 1),
                    fill = ggplot2::guide_legend(nrow = 1),
                    shape = ggplot2::guide_legend(nrow = 1))

  ggplot2::ggsave(
      paste0("results/figures/",
            "fig4_", cluster_method, "_", transformation, "_", gp_log,
            "_", deseq, ".tiff"),
      width = 11, height = 10,
      compression = "lzw+p"
    )

}

args <- commandArgs(trailingOnly = TRUE)

plot_four(args[1], args[2], args[3], args[4])
