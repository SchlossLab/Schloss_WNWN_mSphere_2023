library(tidyverse)

accuracy <- read_tsv("simulation_clusters.tsv") %>%
  rename(subset = fracCorrectPred, all = fracCorrect) %>%
  group_by(fraction, n_seqs, filter, transform, distance, method) %>%
  summarize(subset = mean(subset), all = mean(all), .groups = "drop")

# no difference between RLE, TMM, and upperquartile - run with upperquartile
accuracy %>%
  filter(transform %in% c("RLE", "TMM", "upperquartile")) %>% 
  group_by(method, fraction, n_seqs, filter, transform, distance) %>%
  summarize(diff = max(all) - min(all), .groups = "drop") %>%
  group_by(method) %>%
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
  filter(method == "pam") %>%
  select(-subset) %>%
  pivot_wider(names_from = transform, values_from = all) %>%
  mutate(diff = rarefaction00 - subsample00) %>%
  arrange(desc(diff)) %>%
  filter(filter == "nofilter")




accuracy %>%
  filter(transform %in% c("rarefaction00", "rarefaction15")) %>%
  filter(method == "pam") %>%
  select(-subset) %>%
  pivot_wider(names_from = transform, values_from = all) %>%
  mutate(diff = rarefaction00 - rarefaction15) %>%
  arrange(desc(diff)) %>%
  filter(filter == "nofilter" & n_seqs ==10000) %>%
  print(n = Inf)


accuracy %>%
  filter(distance == "uunifrac" & fraction == "s1" & method == "pam" & filter == "nofilter") %>%
  select(-subset) %>%
  arrange(desc(all)) %>%
  print(n = Inf)

# I feel like their measure of accuracy is garbage. Maybe calculate MCC instead?


read_tsv("simulation_clusters.tsv") %>%
  rename(subset = fracCorrectPred, all = fracCorrect) %>%
  filter(filter == "filter" & !transform %in% c("RLE", "TMM")) %>%
  filter(method == "pam" & distance == "bray") %>%
  filter(fraction == 1) %>%
  group_by(n_seqs, transform, distance) %>%
  summarize(mean = mean(all),
            lci = quantile(all, 0.025),
            uci = quantile(all, 0.975)) %>%
  group_by(n_seqs) %>%
  arrange(mean, .by_group = TRUE) %>%
  print(n = Inf)

read_tsv("simulation_clusters.tsv") %>%
  rename(subset = fracCorrectPred, all = fracCorrect) %>%
  filter(filter == "filter" & !transform %in% c("RLE", "TMM")) %>%
  filter(distance == "bray") %>%
  filter(fraction == "1.15" & n_seqs == 10000) %>%
  group_by(method, n_seqs, transform, distance) %>%
  summarize(mean = mean(all),
            lci = quantile(all, 0.025),
            uci = quantile(all, 0.975)) %>%
  group_by(method) %>%
  arrange(mean, .by_group = TRUE) %>%
  print(n = Inf)

# why is proportion so good? it's only so good with pam - it kind of sucks with
# hclust and kmeans



read_tsv("simulation_adonis.tsv") %>%
  filter(filter == "nofilter" & !transform %in% c("RLE", "TMM")) %>%
  filter(distance == "bray") %>%
  filter(fraction == "s1" & n_seqs == 50000) %>%
  group_by(transform, distance) %>%
  summarize( frac_sig = mean(p_value < 0.05))

# subsample and rarefaction are near 5% false positives, everything else
# is at 100% false positive


read_tsv("simulation_adonis.tsv") %>%
  filter(filter == "nofilter" & !transform %in% c("RLE", "TMM")) %>%
  filter(distance == "bray") %>%
  filter(fraction == "1" & n_seqs == 50000) %>%
  group_by(transform, distance) %>%
  summarize( frac_sig = mean(p_value < 0.05))

# everything near 5% false positives


read_tsv("simulation_adonis.tsv") %>%
  filter(filter == "nofilter" & !transform %in% c("RLE", "TMM")) %>%
  filter(distance == "bray") %>%
  filter(fraction == "1.15" & n_seqs == 50000) %>%
  group_by(transform, distance) %>%
  summarize(frac_sig = mean(p_value < 0.05), .groups = "drop") %>%
  arrange(frac_sig) %>%
  print(n = Inf)
  
# proportion, rarefactions, and subsamples all at 100% power, none,
# uq, and deseq below 10%


read_tsv("simulation_alpha.tsv") %>%
  filter(filter == "filter" & !transform %in% c("RLE", "TMM")) %>%
  filter(fraction == "s1" & n_seqs == 5000) %>%
  group_by(transform) %>%
  summarize(sobs_sig = mean(sobs_p < 0.05),
            shannon_sig = mean(shannon_p < 0.05))
  
# subsample and rarefaction are near 5% false positives, everything else
# is at 100% false positive


read_tsv("simulation_alpha.tsv") %>%
  filter(filter == "filter" & !transform %in% c("RLE", "TMM")) %>%
  filter(fraction == "1" & n_seqs == 5000) %>%
  group_by(transform) %>%
  summarize(sobs_sig = mean(sobs_p < 0.05),
            shannon_sig = mean(shannon_p < 0.05))

# everything near 5% false positives

  
read_tsv("simulation_alpha.tsv") %>%
  filter(filter == "filter" & !transform %in% c("RLE", "TMM")) %>%
  filter(fraction == "1.15" & n_seqs == 5000) %>%
  group_by(transform) %>%
  summarize(sobs_sig = mean(sobs_p < 0.05),
            shannon_sig = mean(shannon_p < 0.05)) #rarefacgtion00 the only one always significant

# rarefactoins both at high powers, everything else below 10% for richness,
# lag in diversity
