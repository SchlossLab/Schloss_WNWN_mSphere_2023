#!/usr/bin/env Rscript

library(tidyverse)
library(here)
library(ggtext)

add_nl <- function(x) {

  paste0("N~L~ = ", format(as.numeric(x), big.mark = "\\\\,"))

}

y_axes <- c(
  median_diff = "Difference between\nrarefaction and subsampling",
  iqr_diff = "Difference between the IQR\nof subsampling and rarefaction"
  )

pretty_distances <- c(bray = "Bray-Curtis",
                      euclidean = "Euclidean",
                      poisson = "Poisson",
                      logFC = "Mean Squared Difference",
                      uunifrac = "Unweighted UniFrac",
                      wunifrac = "Weighted UniFrac")

rare_sub <- read_tsv(here("data/simulation_clusters.tsv.gz")) %>%
  filter(filter == "filter") %>%
  filter(fraction != "s1" & fraction != "1") %>%
  filter(transform %in% c("subsample15", "rarefaction15")) %>%
  filter(distance %in%
            c("bray", "euclidean", "logFC", "poisson",
            "uunifrac", "wunifrac")) %>%
  filter(method == "pam") %>%
  select(fraction, n_seqs, rep, transform, distance, fracCorrect) %>%
  mutate(fraction = as.numeric(fraction),
            distance = pretty_distances[distance])

#the results are qualitatively the same if I do the median accuracy by replicate
#versus by median of replicates
rare_sub %>%
  group_by(fraction, n_seqs, distance, transform) %>%
  summarize(iqr = IQR(fracCorrect),
            median = median(fracCorrect), .groups = "drop") %>%
  pivot_wider(names_from = transform, values_from = c("median", "iqr")) %>%
  mutate(iqr_diff = iqr_subsample15 - iqr_rarefaction15,
            median_diff = median_rarefaction15 - median_subsample15) %>%
  select(fraction, n_seqs, distance, iqr_diff, median_diff) %>%
  pivot_longer(cols = c(iqr_diff, median_diff), names_to = "metric") %>%
  mutate(metric = factor(metric, levels = c("median_diff", "iqr_diff"))) %>%
  ggplot(aes(x = fraction, y = value, color = distance, group = distance,
            shape = distance, fill = distance)) +
  geom_line() +
  geom_point() +
  facet_grid(metric ~ n_seqs, scales = "free_y",
            labeller = labeller(n_seqs = add_nl, metric = y_axes),
            switch = "y") +
  scale_fill_manual(
    values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#000000")
  ) +
  scale_color_manual(
    values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#000000")
  ) +
  scale_shape_manual(
    values = c(21, 22, 23, 24, 25, 21)
  ) +
  coord_cartesian(xlim = c(1, 3.5)) +
  labs(x = "Effect Size",
      y = NULL,
      color = "Distance Calculation:",
      fill = "Distance Calculation:",
      shape = "Distance Calculation:") +
  theme_light() +
  theme(
    legend.position = "top",
    legend.box = "horizontal",
    legend.margin = margin(t = 0, r = -0.5, b = 0, l = -0.5, unit = "cm"),
    strip.text.x = element_markdown(),
    strip.placement = "outside",
    strip.background.y = element_blank(),
    strip.text.y = element_text(color = "black", size = 11)
  ) +
  guides(color = guide_legend(nrow = 1),
        fill = guide_legend(nrow = 1),
        shape = guide_legend(nrow = 1))

ggsave("results/figures/subsample_rarefaction_compare.pdf",
      width = 11, height = 6)
