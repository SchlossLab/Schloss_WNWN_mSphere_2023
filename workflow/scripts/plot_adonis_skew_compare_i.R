#!/usr/bin/env Rscript

library(tidyverse)
library(here)

pretty_distances <- c(bray = "Bray-Curtis",
                      euclidean = "Euclidean",
                      poisson = "Poisson",
                      logFC = "Mean Squared Difference",
                      uunifrac = "Unweighted UniFrac",
                      wunifrac = "Weighted UniFrac")

pretty_transform <- c(deseq = "DeSeq VS",
                    none = "Raw counts",
                    proportion = "Relative abundance",
                    rarefaction00 = "Rarefaction",
                    upperquartile = "Upper Quartile Log Fold Change")

pretty_simulation <- c(a = "GlobalPatterns",
                      log = "Log-distributed")

pretty_fraction <- c("1" = "No",
                    "s1" = "Yes")

adonis_data <- read_tsv("old_data/simulation_adonis.tsv.gz")

p <- adonis_data %>%
  mutate(simulation = str_replace(simulation, ".*_", "")) %>%
  filter(filter == "filter") %>%
  filter(n_seqs == 10000) %>%
  filter((distance == "bray" &
            transform %in% c("proportion", "none", "rarefaction00")) |
          (distance == "euclidean" &
            transform %in% c("none", "deseq", "rarefaction00")) |
          (distance == "poisson" &
            transform %in% c("none", "rarefaction00")) |
          (distance == "logFC" &
            transform %in% c("upperquartile", "rarefaction00")) |
          (distance == "uunifrac" &
            transform %in% c("proportion", "rarefaction00")) |
          (distance == "wunifrac" &
            transform %in% c("none", "proportion", "rarefaction00"))
          ) %>%
  group_by(simulation, fraction, n_seqs, transform, distance) %>%
  summarize(p = mean(p_value <= 0.05), .groups = "drop") %>%
  mutate(distance = factor(pretty_distances[distance],
                        levels = pretty_distances),
      transform = factor(pretty_transform[transform],
                          levels = pretty_transform),
      simulation = factor(pretty_simulation[simulation],
                          levels = pretty_simulation)
  )

p %>%
  filter(fraction == "s1" | fraction == "1") %>%
  mutate(fraction = factor(pretty_fraction[fraction],
                          levels = pretty_fraction)
                          ) %>%
  ggplot(aes(x = fraction, y = p, color = transform, shape = transform,
              fill = transform)) +
  geom_hline(yintercept = 0.05, color = "darkgray") +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  facet_grid(simulation ~ distance) +
#  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(
    values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")
  ) +
  scale_color_manual(
    values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")
  ) +
  scale_shape_manual(
    values = c(21, 22, 23, 24, 25, 21)
  ) +
  coord_cartesian(ylim = c(0, 1.05)) +
  labs(x = "Confounded by sample size",
      y = "Fraction of significant tests",
      color = "Normalization Method:",
      fill = "Normalization Method:",
      shape = "Normalization Method:") +
  theme_light() +
  theme(
    legend.position = "top",
    legend.box = "horizontal",
    legend.margin = margin(t = 0, r = -0.5, b = 0, l = -0.5, unit = "cm")
  ) +
  guides(color = guide_legend(nrow = 1),
        fill = guide_legend(nrow = 1),
        shape = guide_legend(nrow = 1))

ggsave("results/figures/adonis_skew_compare_i.pdf", width = 11, height = 5)
