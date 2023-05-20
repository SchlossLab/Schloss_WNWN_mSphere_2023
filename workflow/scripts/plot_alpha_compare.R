#!/usr/bin/env Rscript

library(tidyverse)
library(ggh4x)
library(ggtext)
library(here)


add_nl <- function(x) {
  paste0("&Ntilde;~L~ = ", format(as.numeric(x), big.mark = "\\\\,"))
}

pretty_alpha <- c(sobs_p = "Richness",
                 shannon_p = "Shannon diversity")

pretty_transform <- c(deseq = "DESeq VS",
                    none = "Raw counts",
                    proportion = "Relative abundance",
                    subsample15 = "Sub-sampled (15%)",
                    rarefaction00 = "Rarefaction",
                    upperquartile = "Upper Quartile Log Fold Change")

pretty_distro <- c(random = "Randomized",
                  skew = "Confounded")


plot_alpha_compare <- function(gp_log) {

  alpha <- read_tsv(here("data/simulation_alpha.tsv.gz")) %>%
    filter(model == gp_log & filter == "filter") %>%
    filter(transform %in%
        c("proportion", "none", "rarefaction00", "deseq", "upperquartile")) %>%
    group_by(distribution, fraction, n_seqs, transform) %>%
    summarize(sobs_p = mean(sobs_p <= 0.05),
              shannon_p = mean(shannon_p <= 0.05), .groups = "drop") %>%
    pivot_longer(cols = c(sobs_p, shannon_p)) %>%
    drop_na() %>%
    filter(!(name == "shannon_p" & transform == "deseq")) %>%
    mutate(name = factor(name, levels = c("sobs_p", "shannon_p")),
          transform = factor(pretty_transform[transform],
                                levels = pretty_transform))

  alpha %>%
    ggplot(aes(x = as.numeric(fraction), y = value, shape = transform,
              group = transform, color = transform, fill = transform)) +
    geom_line(position = position_dodge(width = 0.075)) +
    geom_point(position = position_dodge(width = 0.075),
              size = 2) +
    facet_nested(n_seqs ~ name + distribution,
                labeller = labeller(name = pretty_alpha, n_seqs = add_nl,
                            distribution = pretty_distro),
                strip = strip_nested(
                  text_x = elem_list_text(colour = c("#000000", "#FFFFFF")),
                  background_x = elem_list_rect(fill = c("#FFFFFF", "grey70")),
                  by_layer_x = TRUE
                )) +
      scale_fill_manual(
      values = c("#d95f02", "#7570b3", "#e7298a", "#66a61e")
    ) +
    scale_color_manual(
      values = c("#d95f02", "#7570b3", "#e7298a", "#66a61e")
    ) +
    scale_shape_manual(
      values = c(22, 23, 24, 25)
    ) +
    scale_x_continuous(breaks = c(1, 2, 3)) +
    labs(x = "Effect Size",
        y = "Fraction of significant tests",
        color = NULL, #"Normalization Method:",
        fill = NULL, #"Normalization Method:",
        shape = NULL #"Normalization Method:"
        ) +
    theme_light() +
    theme(
      legend.position = "top",
      legend.box = "horizontal",
      legend.margin = margin(t = 0, r = -0.5, b = 0, l = -0.5, unit = "cm"),
      strip.text.y = element_markdown()
    ) +
    guides(color = guide_legend(nrow = 2),
          fill = guide_legend(nrow = 2),
          shape = guide_legend(nrow = 2))

  ggsave(paste0("results/figures/alpha_compare_", gp_log, ".tiff"),
        width = 6, height = 8,
        compression = "lzw+p")
}

args <- commandArgs(trailingOnly = TRUE)
plot_alpha_compare(args[1])