#!/usr/bin/env Rscript

library(tidyverse)
library(here)
library(ggtext)

read_tsv(here("data/simulated_nseqs_distros.tsv")) %>%
  mutate(simulation = factor(simulation,
                             levels = c("Log-scaled", "GlobalPatterns"))) %>%
  ggplot(aes(x = n_seqs, y = simulation)) +
  geom_jitter(height = 0.25, width = 0.0, alpha = 0.25) +
  scale_y_discrete(name = NULL) +
  scale_x_continuous(limits = c(0, 2.5e6),
                    breaks = seq(0, 2.5e6, 0.5e6),
                    labels = seq(0, 2.5, 0.5),
                    name = "Number of sequences per sample (x10^6^)") +
  theme_classic() +
  theme(
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x = element_markdown()
        )

ggsave("results/figures/distribution_shape.pdf",
        width = 4, height = 2)
