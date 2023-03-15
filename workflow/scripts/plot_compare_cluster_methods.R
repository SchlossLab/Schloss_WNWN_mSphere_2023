#!/usr/bin/env Rscript

library(tidyverse)
library(ggtext)

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
                    upperquartile = "Upper Quartile Log Fold Change",
                    rarefaction15 = "Rarefaction (15%)")

pretty_method <- c(pam = "PAM",
                  kmeans = "K-means",
                  hclust = "HClust")

best_method <- read_tsv("data/simulation_clusters.tsv.gz") %>%
  filter(filter == "filter") %>%
  filter(fraction != "s1" & fraction != "1") %>%
  filter((distance == "bray" &
            transform %in% c("proportion", "subsample15", "rarefaction15", "none", "deseq")) |
        (distance == "euclidean" &
            transform %in% c("none", "deseq", "subsample15", "rarefaction15")) |
        (distance == "poisson" &
            transform %in% c("subsample15", "none", "rarefaction15")) |
        (distance == "logFC" &
            transform %in% c("upperquartile", "subsample15", "rarefaction15")) |
        (distance == "uunifrac" &
            transform %in% c("subsample15", "proportion", "rarefaction15")) |
        (distance == "wunifrac" &
            transform %in% c("subsample15", "none", "proportion", "deseq", "rarefaction15"))
        ) %>%
  filter(fraction == 1.15) %>%
  select(n_seqs, rep, transform, distance, method, fracCorrect) %>%
    group_by(n_seqs, transform, distance, rep) %>%
    slice_max(fracCorrect) %>%
    ungroup() %>%
    mutate(fracCorrect = 1) %>%
    pivot_wider(names_from = method, values_from = fracCorrect, values_fill = 0)

# overall <- best_method %>%
#               summarize(pam = 100 * mean(pam),
#                         kmeans = 100 * mean(kmeans),
#                         hclust = 100 * mean(hclust))

# by_n_seqs <- best_method %>%
#                 group_by(n_seqs) %>%
#                 summarize(pam = 100 * mean(pam),
#                           kmeans = 100 * mean(kmeans),
#                           hclust = 100 * mean(hclust))

best_method %>%
  group_by(n_seqs, transform, distance) %>%
  summarize(pam = 100 * mean(pam),
          kmeans = 100 * mean(kmeans),
          hclust = 100 * mean(hclust), .groups = "drop") %>%
  pivot_longer(c(pam, kmeans, hclust),
              names_to = "method", values_to = "freq") %>%
  mutate(method = factor(pretty_method[method],
                          levels = pretty_method),
        distance = factor(pretty_distances[distance],
                          levels = pretty_distances),
        transform = factor(pretty_transform[transform],
                          levels = pretty_transform)) %>%

  ggplot(aes(x = method, y = freq,
            color = transform, fill = transform,
            shape = transform, group = transform)) +
  geom_point(position = position_dodge(0.25), size = 2) +
  geom_line(position = position_dodge(0.25)) +
  facet_grid(n_seqs ~ distance,
            labeller = labeller(n_seqs = add_nl)) +
  # scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(
    values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#000000")
  ) +
  scale_color_manual(
    values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#000000")
  ) +
  scale_shape_manual(
    values = c(21, 22, 23, 24, 25, 21)
  ) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "Clustering Method",
      y = "Percentage of Simulations (N=100)",
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

ggsave("results/figures/compare_cluster_methods.pdf", height = 10, width = 11)