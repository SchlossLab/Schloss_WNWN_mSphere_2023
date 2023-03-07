#!/usr/bin/env Rscript

library(tidyverse)
library(here)
library(ggtext)

read_tsv(here("data/globalpatterns_distribution.tsv")) %>%
  ggplot(aes(x = n_seqs, y=1)) +
  geom_jitter(height = 0.25, width = 0.0) +
  scale_y_continuous(limits = c(0, 2), name = NULL) +
  scale_x_continuous(limits = c(0, 2.5e6), 
                    breaks = seq(0, 2.5e6, 0.5e6),
                    labels = seq(0, 2.5, 0.5),
                    name = "Number of sequences per sample (x10^6^)") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x = element_markdown()
        )

ggsave("results/figures/globalpatterns_distribution.pdf",
        width = 4, height = 2)
