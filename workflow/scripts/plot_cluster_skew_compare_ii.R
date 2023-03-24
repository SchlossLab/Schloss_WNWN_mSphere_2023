#!/usr/bin/env Rscript

library(tidyverse)
library(here)
library(ggtext)


set.seed(19760620)

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

pretty_simulation <- c(sim_a = "GlobalPatterns",
                      sim_log = "Log-distributed")

# need to get global patterns skew for all effect sizes to compare
# type i/type ii errors

cluster_data <- read_tsv(here("data/simulation_clusters.tsv.gz"))

cluster_data %>%
  filter(filter == "filter") %>%
  filter(method == "kmeans") %>%
  filter(fraction != "s1" & fraction != "1") %>%
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
  group_by(simulation, fraction, transform, distance) %>%
  summarize(mean = mean(fracCorrect, na.rm = TRUE),
            lci = quantile(fracCorrect, 0.025, na.rm = TRUE),
            uci = quantile(fracCorrect, 0.975, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(distance = factor(pretty_distances[distance],
                          levels = pretty_distances),
        transform = factor(pretty_transform[transform],
                            levels = pretty_transform),
        simulation = factor(pretty_simulation[simulation],
                            levels = pretty_simulation)) %>%
  ggplot(aes(x = fraction, y = mean,
            color = transform, fill = transform, group = transform)) +
    geom_point() +
    geom_line() +
    facet_grid(simulation ~ distance) +
    theme_light()

ggsave("results/figures/cluster_skew_compare_ii.pdf", width = 11, height = 5)
