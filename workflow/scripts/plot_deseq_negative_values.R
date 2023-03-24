#!/usr/bin/env Rscript

library(tidyverse)
library(ggtext)

add_nl <- function(x) {

  paste0("N~L~ = ", format(as.numeric(x), big.mark = "\\\\,"))

}

pretty_metric = c(
    mean_a = "Mean fraction of negative\ndistances for each simulation",
    frac_negative = "Fraction of simulations\nwith a negative distance"
    )

read_tsv("data/deseq_negative_values.tsv.gz") %>%
  filter(filter == "filter") %>%
  group_by(fraction, n_seqs) %>%
  summarize(mean_a = mean(sim_a_fraction),
            lq_a = quantile(sim_a_fraction, prob = 0.025),
            uq_a = quantile(sim_a_fraction, prob = 0.975),
            frac_negative = mean(sim_a_negative),
            .groups = "drop") %>%
  pivot_longer(c(mean_a, frac_negative),
              names_to = "metric", values_to = "value") %>%
  mutate(lq_a = if_else(metric == "frac_negative", 0, lq_a),
        uq_a = if_else(metric == "frac_negative", 0, uq_a),
        metric = factor(pretty_metric[metric], levels = pretty_metric)
        ) %>%
  ggplot(aes(x = fraction, y = value, group = 1)) +
  geom_line() +
  geom_point() +
  geom_linerange(aes(ymin = lq_a, ymax = uq_a), alpha = 0.25) +
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
