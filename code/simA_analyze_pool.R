library(tidyverse)

accuracy <- read_tsv("data/simulation_cluster_accuracy.tsv") %>%
  group_by(fraction, n_seqs, filter, transform, distance, cluster) %>%
  summarize(subset = mean(subset), all = mean(all), .groups = "drop")

# no difference between RLE, TMM, and upperquartile - run with upperquartile
accuracy %>%
  filter(transform %in% c("RLE", "TMM", "upperquartile")) %>% 
  group_by(cluster, fraction, n_seqs, filter, transform, distance) %>%
  summarize(diff = max(all) - min(all), .groups = "drop") %>%
  group_by(cluster) %>%
  summarize(min_diff = min(diff, na.rm = TRUE),
            max_diff = max(diff, na.rm = TRUE))


accuracy <- accuracy %>%
  filter(!transform %in% c("RLE", "TMM"))


# difference between bcv and logFC? yes
accuracy %>%
  select(-subset) %>%
  filter(distance == "bcv" | distance == "logFC") %>%
  pivot_wider(names_from = "distance", values_from = "all") %>%
  mutate(diff = bcv - logFC) %>%
  filter(transform == "upperquartile") %>%
  arrange(desc(diff))


# difference between filter or no filter
# difference between 00 and 15
# difference between rarefaction and subsample

accuracy %>%
  filter(transform %in% c("rarefaction00", "subsample00")) %>%
  filter(cluster == "pam") %>%
  select(-subset) %>%
  pivot_wider(names_from = transform, values_from = all) %>%
  mutate(diff = rarefaction00 - subsample00) %>%
  arrange(desc(diff)) %>%
  filter(filter == "nofilter")




accuracy %>%
  filter(transform %in% c("rarefaction00", "rarefaction15")) %>%
  filter(cluster == "pam") %>%
  select(-subset) %>%
  pivot_wider(names_from = transform, values_from = all) %>%
  mutate(diff = rarefaction00 - rarefaction15) %>%
  arrange(desc(diff)) %>%
  filter(filter == "nofilter" & n_seqs ==10000) %>%
  print(n = Inf)


accuracy %>%
  filter(distance == "uunifrac" & fraction == "s1" & cluster == "pam" & filter == "nofilter") %>%
  select(-subset) %>%
  arrange(desc(all)) %>%
  print(n = Inf)

# I feel like their measure of accuracy is garbage. Maybe calculate MCC instead?

pretty_distances <- c(bray = "Bray-Curtis",
                      euclidean = "Euclidean",
                      poisson = "Poisson",
                      logFC = "Mean Squared Diference",
                      uunifrac = "Unweighted UniFrac",
                      wunifrac = "Weighted UniFrac")

pretty_transform <- c(deseq = "DeSeq VS",
                    none = "Raw counts",
                    proportion = "Relative abundance",
                    subsample15 = "Sub-sampled (15%)",
                    upperquartile = "Upper Quartile Log Fold Change",
                    rarefaction00 = "Rarefaction")

read_tsv("data/simulation_cluster_accuracy.tsv") %>%
  filter(filter == "filter" & !transform %in% c("RLE", "TMM")) %>%
  filter(cluster == "pam") %>%
  filter(fraction != "s1") %>%
#  filter(n_seqs %in% c(1000, 2000, 10000)) %>%
  mutate(fraction = as.numeric(fraction)) %>%
  filter((distance == "bray" &
            transform %in% c("proportion", "subsample15", "none", "deseq", "rarefaction00")) |
         (distance == "euclidean" &
            transform %in% c("none", "deseq", "subsample15", "rarefaction00")) |
         (distance == "poisson" &
            transform %in% c("subsample15", "none", "rarefaction00")) |
         (distance == "logFC" &
            transform %in% c("upperquartile", "subsample15", "rarefaction00")) |
         (distance == "uunifrac" &
            transform %in% c("subsample15", "proportion", "rarefaction00")) |
         (distance == "wunifrac" &
            transform %in% c("subsample15", "none", "proportion", "deseq", "rarefaction00"))
         ) %>%
  group_by(fraction, n_seqs, transform, distance) %>%
  summarize(mean = mean(all), sd = sd(all), n = n()) %>%
  mutate(distance = factor(pretty_distances[distance],
                           levels = pretty_distances),
         transform = factor(pretty_transform[transform],
                            levels = pretty_transform)) %>%
  ggplot(aes(x = fraction, y = mean, group = transform,
             color = transform, shape = transform, fill = transform)) +
  geom_line(position = position_dodge(width = 0.075)) +
  geom_linerange(aes(ymin = mean - sd, ymax = mean + sd),
                 alpha = 0.6, position = position_dodge(width = 0.075),
                 show.legend = FALSE) +
  geom_point(position = position_dodge(width = 0.075),
             size = 2) +
  facet_grid(n_seqs ~ distance) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(
    values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#000000")
  ) +
  scale_color_manual(
    values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#000000")
  ) +
  scale_shape_manual(
    values = c(21, 22, 23, 24, 25, 21)
  ) +
  coord_cartesian(ylim = c(0.4, 1.05)) +
  labs(x = "Effect Size",
       y = "Accuracy",
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

ggsave("figure_4.pdf", width = 11, height = 10)



read_tsv("data/simulation_cluster_accuracy.tsv") %>%
  filter(filter == "filter" & !transform %in% c("RLE", "TMM")) %>%
  filter(cluster == "pam" & distance == "bray") %>%
  filter(fraction == 1) %>%
  group_by(n_seqs, transform, distance) %>%
  summarize(mean = mean(all), sd = sd(all)) %>%
  group_by(n_seqs) %>%
  arrange(mean, .by_group = TRUE) %>%
  print(n = 50)

read_tsv("data/simulation_cluster_accuracy.tsv") %>%
  filter(filter == "filter" & !transform %in% c("RLE", "TMM")) %>%
  filter(distance == "bray") %>%
  filter(fraction == "1.15" & n_seqs == 10000) %>%
  group_by(cluster, n_seqs, transform, distance) %>%
  summarize(mean = mean(all), sd = sd(all)) %>%
  group_by(cluster) %>%
  arrange(mean, .by_group = TRUE) %>%
  print(n = 50)

# why is proportion so good? it's only so good with pam - it kind of sucks with
# hclust and kmeans



read_tsv("data/simulation_adonis.tsv") %>%
  filter(filter == "filter" & !transform %in% c("RLE", "TMM")) %>%
  filter(distance == "bray") %>%
  filter(fraction == "s1" & n_seqs == 5000) %>%
  group_by(transform, distance) %>%
  summarize( frac_sig = mean(p_value < 0.05))

read_tsv("data/simulation_adonis.tsv") %>%
  filter(filter == "filter" & !transform %in% c("RLE", "TMM")) %>%
  filter(distance == "bray") %>%
  filter(fraction == "1.15" & n_seqs == 5000) %>%
  group_by(transform, distance) %>%
  summarize(frac_sig = mean(p_value < 0.05))


read_tsv("data/simulation_alpha.tsv") %>%
  filter(filter == "filter" & !transform %in% c("RLE", "TMM")) %>%
  filter(fraction == "s1" & n_seqs == 5000) %>%
  group_by(transform) %>%
  summarize(sobs_sig = mean(sobs_p < 0.05),
            shannon_sig = mean(shannon_p < 0.05))

read_tsv("data/simulation_alpha.tsv") %>%
  filter(filter == "filter" & !transform %in% c("RLE", "TMM")) %>%
  filter(fraction == "1" & n_seqs == 5000) %>%
  group_by(transform) %>%
  summarize(sobs_sig = mean(sobs_p < 0.05),
            shannon_sig = mean(shannon_p < 0.05))
  
read_tsv("data/simulation_alpha.tsv") %>%
  filter(filter == "filter" & !transform %in% c("RLE", "TMM")) %>%
  filter(fraction == "1.15" & n_seqs == 5000) %>%
  group_by(transform) %>%
  summarize(sobs_sig = mean(sobs_p < 0.05),
            shannon_sig = mean(shannon_p < 0.05)) #rarefacgtion00 the only one always significant

