#!/usr/bin/env Rscript

library(tidyverse)
library(ggh4x)
library(here)

pretty_distances <- c(bray = "Bray-Curtis",
                      euclidean = "Euclidean",
                      poisson = "Poisson",
                      logFC = "Mean Squared Difference",
                      uunifrac = "Unweighted UniFrac",
                      wunifrac = "Weighted UniFrac")

pretty_transform <- c(deseq = "DESeq VS",
                    none = "Raw counts",
                    proportion = "Relative abundance",
                    rarefaction00 = "Rarefaction",
                    upperquartile = "Upper Quartile Log Fold Change")

pretty_model <- c(gp = "GlobalPatterns",
                  log = "Log-scaled")

pretty_distribution <- c("random" = "No",
                        "skew" = "Yes")

adonis_data <- read_tsv("data/simulation_adonis.tsv.gz") %>%
  filter(filter == "filter" & n_seqs == 10000 &
          (fraction == 1 | fraction == 1.15)) %>%
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
  group_by(model, fraction, distribution, transform, distance) %>%
  summarize(p = mean(p_value <= 0.05), .groups = "drop") %>%
  mutate(distance = factor(pretty_distances[distance],
                        levels = pretty_distances),
      transform = factor(pretty_transform[transform],
                          levels = pretty_transform),
      model = factor(pretty_model[model],
                          levels = pretty_model),
      distribution = factor(pretty_distribution[distribution],
                          levels = pretty_distribution),
      fraction = paste0("ES = ", format(fraction))
  )

adonis_data %>%
  ggplot(aes(x = distribution, y = p, color = transform, shape = transform,
            fill = transform)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  facet_nested(fraction + model ~ distance,
              strip = strip_nested(
              text_y = elem_list_text(colour = c("#000000", "#FFFFFF")),
              background_y = elem_list_rect(fill = c("#FFFFFF", "grey70")),
              by_layer_y = TRUE
            )) +
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
      color = NULL, #"Normalization Method:",
      fill = NULL, #"Normalization Method:",
      shape = NULL #"Normalization Method:"
      ) +
  theme_light() +
  theme(
    legend.position = "top",
    legend.box = "horizontal",
    legend.margin = margin(t = 0, r = -0.5, b = 0, l = -0.5, unit = "cm")
  ) +
  guides(color = guide_legend(nrow = 1),
        fill = guide_legend(nrow = 1),
        shape = guide_legend(nrow = 1))

ggsave("results/figures/adonis_skew_compare.tiff",
      width = 11, height = 7,
      compression = "lzw+p")
