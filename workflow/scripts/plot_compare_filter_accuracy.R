#!/usr/bin/env Rscript

library(tidyverse)
library(ggtext)
library(here)

pretty_distances <- c(bray = "Bray-Curtis",
                      euclidean = "Euclidean",
                      poisson = "Poisson",
                      logFC = "Mean Squared<br>Difference",
                      uunifrac = "Unweighted<br>UniFrac",
                      wunifrac = "Weighted<br>UniFrac")

pretty_transform <- c(deseq = "DeSeq VS",
                    none = "Raw counts",
                    proportion = "Relative abundance",
                    rarefaction00 = "Rarefaction",
                    upperquartile = "UQ Log Fold Change")

pretty_models <- c(gp = "Global Patterns",
                  log = "Log-distributed")


cluster_data <- read_tsv(here("data/simulation_clusters.tsv.gz")) %>%
  filter(n_seqs == 10000 & fraction != 1 & method == "kmeans" &
          distribution == "random") %>%
  filter((distance == "bray" &
            transform %in% c("proportion", "rarefaction00", "none")) |
        (distance == "euclidean" &
            transform %in% c("none", "rarefaction00", "deseq")) |
        (distance == "poisson" &
            transform %in% c("rarefaction00", "none")) |
        (distance == "logFC" &
            transform %in% c("upperquartile", "rarefaction00")) |
        (distance == "uunifrac" &
            transform %in% c("rarefaction00", "proportion")) |
        (distance == "wunifrac" &
            transform %in% c("rarefaction00", "none", "proportion"))
        ) %>%
  select(-fracCorrectPred) %>%
  pivot_wider(names_from = filter, values_from = fracCorrect) %>%
  mutate(f_nf = abs(filter - nofilter),
        fraction = as.numeric(fraction)) %>%
  group_by(model, distribution, fraction, n_seqs, transform, distance) %>%
  summarize(mean = mean(f_nf),
            lci = quantile(f_nf, prob = 0.025),
            uci = quantile(f_nf, prob = 0.975), .groups = "drop") %>%
  mutate(distance = factor(pretty_distances[distance],
              levels = pretty_distances),
          transform = factor(pretty_transform[transform],
                levels = pretty_transform),
          model = factor(pretty_models[model],
                levels = pretty_models))

cluster_data %>%
  ggplot(aes(x = fraction, y = mean,
            group = transform, color = transform,
            shape = transform, fill = transform)) +
  geom_line(position = position_dodge(width = 0.075), show.legend = FALSE) +
  geom_point(position = position_dodge(width = 0.075), size = 2) +
  geom_linerange(aes(ymin = lci, ymax = uci),
                position = position_dodge(width=0.075), alpha = 0.6,
                show.legend = FALSE) +
  facet_grid(distance ~ model) +
  scale_fill_manual(
    values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#000000")
  ) +
  scale_color_manual(
    values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#000000")
  ) +
  scale_shape_manual(
    values = c(21, 22, 23, 24, 25, 21)
  ) +
  scale_x_continuous(
    breaks = seq(0, 3.5, 0.5)
  ) +
  labs(x = "Effect size",
      y = "Absolute value of difference in accuracy",
      color = NULL,
      fill = NULL,
      shape = NULL) +
  theme_light() +
  theme(
    legend.position = "top",
    legend.box = "horizontal",
    legend.spacing.x = unit(0.05, "cm"),
    legend.margin = margin(t = 0, r = -0.5, b = 0, l = -0.5, unit = "cm"),
    strip.text.y = element_markdown()
  )

ggsave("results/figures/compare_filter_accuracy.pdf", width = 6, height = 8)
