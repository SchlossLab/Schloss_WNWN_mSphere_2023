#!/usr/bin/env Rscript

library(tidyverse)
library(ggtext)

add_nl <- function(x) {

  paste0("&Ntilde;~L~ = ", format(as.numeric(x), big.mark = "\\\\,"))

}

pretty_metric <- c(
    mean = "Mean fraction of negative\ndistances for each simulation",
    has_negative = "Fraction of simulations\nwith a negative distance"
    )

read_tsv("data/deseq_negative_values.tsv.gz") %>%
  filter(filter == "filter") %>%
  filter(model == "gp") %>%
  filter(distribution == "random") %>%
  pivot_longer(c(mean, has_negative),
              names_to = "metric", values_to = "value") %>%
  mutate(lci = if_else(metric == "has_negative", 0, lci),
        uci = if_else(metric == "has_negative", 0, uci),
        metric = factor(pretty_metric[metric], levels = pretty_metric)
        ) %>%
  ggplot(aes(x = effect_size, y = value, group = 1)) +
  geom_line() +
  geom_point() +
  geom_linerange(aes(ymin = lci, ymax = uci), alpha = 0.25) +
  labs(x = "Effect Size",
      y = NULL) +
  facet_grid(metric ~ n_seqs,
            switch = "y", labeller = labeller(n_seqs = add_nl)) +
  theme_light() +
  theme(
    strip.placement = "outside",
    strip.text.x = element_markdown(),
    strip.text.y = element_markdown(color = "black"),
    strip.background.y = element_blank()
  )

ggsave("results/figures/deseq_negative_values.pdf", width = 8, height = 5)
