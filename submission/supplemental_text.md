---
bibliography: references.bib
output: pdf_document
csl: asm.csl
geometry: margin=1.0in
header-includes:
 - \usepackage{setspace}
 - \doublespacing
 - \usepackage{helvet} % Helvetica font
 - \renewcommand*\familydefault{\sfdefault} % Use the sans serif version of the font
 - \usepackage[T1]{fontenc}
 - \usepackage[shortcuts]{extdash}
 - \usepackage[document]{ragged2e}
---

## Supplementary text

Quotes from selected papers that cited WNWN and conflated "rarefying"/"rarefy" and "rarefaction". Many of these papers also attributed "rarefaction" to WNWN. Rarefying, rarefy, and rarefaction are bolded to make the words easier to identify. All references have been updated to reflect the numbering in the current study. Citations to WNWN are indicated by [@McMurdie2014] and all numbering is in reference to the list of references in this Supplement. Quotes are ordered chronologically.

* "Additionally, **rarefaction** has recently been shown to introduce errors in analyses, and alternatives to rarefaction have been proposed [@McMurdie2014]" - Quoted from [@Goodrich2014]. Implies that rarefaction and WNWN's rarefying are the same thing.

* "**Rarefaction** is analytically problematic and poses multiple statistical problems: (i) omission of available valid data, (ii) the estimation of overdispersion is more difficult due to data loss, (iii) loss of power (type II error), (iv) dependence on an arbitrary threshold and (v) additional uncertainty due to the randomness in **rarefaction** [@McMurdie2014]." - Quoted from [@Blint2016]. Implies that rarefaction and WNWN's rarefying are the same thing.

* "McMurdie and Holmes [@McMurdie2014] penalized the **rarefying** technique for dropping the lowest fifteenth percentile of sample library sizes in their simulations by counting the dropped samples as "incorrectly clustered." Because the 15th percentile was used to set **rarefaction** depth, this capped clustering accuracy at 85%." - Quoted from [@Weiss2017]. An example of a paper from the QIIME development community using rarefaction and rarefying interchangeably.

* "**Rarefaction** has a limited ability in this regard since the total sum constraint still exists after **rarefaction**. In addition, it suffers from a great power loss due to the discard of a large number of reads [@McMurdie2014]." - Quoted from [@Chen2018]. Implies that rarefaction and WNWN's rarefying are the same thing.

* "Another popular normalization approach is **rarefaction**, which consists on subsampling the same number of reads for each sample so that all samples have the same number of total counts. **Rarefaction** is not recommended because it entails the loss of important information [@McMurdie2014]." - Quoted from [@Calle2019]. Implies that rarefaction and WNWN's rarefying are the same thing.

* "Unfortunately, **rarefaction** is neither justifiable nor necessary, a view framed statistically by McMurdie and Holmes [@McMurdie2014] in the context of comparison of relative abundances." - Quoted from [@Willis2019]. Implies that rarefaction and WNWN's rarefying are the same thing.

* "Some studies perform **rarefaction** to adjust for differences in library size due to unexhaustive metagenomic sampling. Although several pipelines provide this functionality, it has been found inadmissible for metagenomics microbiome studies as it discards many reads leading to decreased sensitivity in differential abundance testing [@McMurdie2014] and biased estimates for alpha diversity [@Willis2019]." - Quoted from [@Eetemadi2020]. Implies that rarefaction and WNWN's rarefying are the same thing.

* "Another widespread practice is to **rarefy** the count data to force the samples to have the same number of total sequence reads [@Hughes2005], at the expense of discarding vast amounts of information [@McMurdie2014]." - Quoted from [@Leite2020]. Hughes and Hellmann [@Hughes2005] describe traditional rarefaction with multiple subsamplings.

* "Note that the terms **rarefying** and **rarefaction** are used interchangeably in microbiome literature [@McMurdie2014]. **Rarefying** was first recommended for microbiome data to deal with rare taxa [@Lozupone2010], which impact some measures of alpha and beta diversities [@McMurdie2014]." - Quoted from [@Lin2020]. As the quote indicates they use rarefying and rarefaction interchangeably although the cited WNWN indicates that they are not the same. In addition, the first review article to advocate the use of rarefaction was Hughes et al. [@Hughes2001] and examples of its use can be found earlier [@Dunbar1999; @Moyer1998];

* "Ways to tackle this include total library size normalization and **rarefaction**, with both remaining debated to date [@McMurdie2014; @Weiss2017]." - Quoted from [@Alteio2021]. Implies that rarefaction and WNWN's rarefying are the same thing.

* "**Rarefaction** is a widely used normalization technique that involves the random subsampling of sequences from the initial sample library to a selected normalized library size. This process is often dismissed as statistically invalid because subsampling effectively discards a portion of the observed sequences, yet it remains prevalent in practice and the suitability of **rarefying**, relative to many other normalization approaches, for diversity analysis has been argued" - Quoted by [@Cameron2021]. Uses rarefaction and rarefying interchangeably. This paper describes a new method as "repeated rarefying", which is the traditional rarefaction approach.

* "Confronted with technical variation as well as the overall increase in raw sequencing data generated per sample over the years, **rarefaction** (or downsampling) was suggested to standardize within and across dataset comparisons ...  However, sequencing depth-based downsizing procedures were soon criticized, not only for being wasteful and discarding information on low-abundant taxa [@McMurdie2014], but also for being unsuited when applied to communities characterized by substantial variation in cell density [@Vandeputte2017]." - Quoted by [@LlornsRico2021].  Implies that rarefaction and WNWN's rarefying are the same thing.

* "Some have suggested that data should be **rarefied** [@Weiss2017; @McKnight2018] to a common sampling depth, typically to the level of the sample with fewest sequences, while others argue that such **rarefaction** is 'inadmissible' and favor approaches that transform or scale sequence counts [@McMurdie2014; @Gloor2017]." - Quoted by [@Neu2021]. Uses rarefied and rarefaction interchangeably and implies that both are the same as WNWN's rarefying.

* "**Rarefaction** has been criticized for wasting data since we effectively remove a portion of the data in the downsampling procedure [@McMurdie2014]." - Quoted by [@Hong2022]. Here and elsewhere in this paper the authors use rarefaction and rarefying interchangeably and implies that both are the same as WNWN's rarefying.

* "**Rarefying** (also referred to as **rarefaction**) is a popular but widely criticized technique for correcting uneven sequencing depths. It involves randomly discarding counts in samples until all samples have the same predefined number of total counts." - Quoted by [@Swift2022].  As the quote indicates they use rarefying and rarefaction interchangeably although the cited WNWN indicates that they are not the same. 

* "One commonly used method is to **rarefy** the data; that is, ASVs or OTUs within a sample are randomly subsampled without replacement to a preselected depth that is the same across all samples. The outcome of this is that all samples will have the same number of ASVs and any samples with fewer ASVs than the **rarefaction** level will be removed from the dataset. The level for **rarefaction** can be decided using a **rarefaction** curve, a method in which each sample is subsampled at multiple levels (e.g. 1,000 reads, 2,000 reads, 3,000 reads...), and the number of unique features or another metric of individual sample diversity of each sample at each level is measured and plotted. When the plot begins to level off after an initial climb up, the corresponding number of sequences indicates an appropriate sampling depth. The appropriate number to **rarefy** must then be balanced with the number of samples that may be dropped from the dataset which do not meet that minimum. An advantage of **rarefaction** is that it may be a more appropriate measure of very low-abundance ("rare") ASVs. This can in turn increase the accuracy of the data, as low biomass samples often have contamination and quality concerns [@Kennedy2014]. There are also disadvantages to this method, the most obvious of which is the discarding of valuable data. Clearly, this is less than ideal as the researcher must pay for the samples and sequences, and in cases where the samples are very valuable or difficult to obtain the loss of data may be destructive to the overall experimental integrity. Additionally, the loss of statistical power by removing sequences from a sample could lead to a loss of differences between two samples [@McMurdie2014]. The statistical consequences extend beyond this, as **rarefying** equalizes sample variance by adding artificial uncertainty [@McMurdie2014]." - Quoted by [@Weinroth2022]. The authors use rarefaction and rarefying interchangeably and implies that both are the same as WNWN's rarefying.

\newpage

### References

\setlength{\parindent}{-0.25in}
\setlength{\leftskip}{0.25in}
\noindent

<div id="refs"></div>
\bibliography{ref}
\setlength{\parindent}{0in}
\setlength{\leftskip}{0in}

\newpage

\includegraphics{figure_S1.png}

**Figure S1. Re-running the R markdown files provided in Protocol S1 of WNWN qualitatively reproduced Figure 4 from WNWN.** 

\newpage

\includegraphics{figure_S2.png}

**Figure S2. Re-running the R markdown files provided in Protocol S1 of WNWN qualitatively reproduced Figure 5 from WNWN.**

\newpage

\includegraphics{figure_S3.png}

**Figure S3. Successful reimplementation and expansion of analysis presented in Figure 4 from WNWN in a Snakemake pipeline.** The reimplemented workflow largely borrowed from the original `simulation-cluster-accuracy-server.Rmd` R markdown file provided in WNWN's Protocol S1. The most notable differences include the use of 100 rather than 5 randomizations and the addition of the median sequencing depth (&Ntilde;~L~) of 50,000. The plotting symbols indicate the median of 100 randomizations and the error bars represent the observed 95% confidence interval. Simulations run at the same effect size are dodged to better reveal overlapping data. The sequencing depths were drawn from the GlobalPatterns dataset and sequences were clustered using PAM.

\newpage

\includegraphics{figure_S4.png}

**Figure S4. Successful reimplementation and expansion of analysis presented in Figure 5 from WNWN in a Snakemake pipeline.** The reimplemented workflow largely borrowed from the original `simulation-cluster-accuracy-server.Rmd` R markdown file provided in WNWN's Protocol S1. The most notable differences include the use of 100 rather than 5 randomizations and the addition of the median sequencing depth (&Ntilde;~L~) of 50,000. The plotting symbols indicate the median of 100 randomizations and the error bars represent the observed 95% confidence interval. Simulations run at the same effect size are dodged to better reveal overlapping data. A light gray line is shown to indicate the best possible accuracy for each library size minimum quantile value. The sequencing depths were drawn from the GlobalPatterns dataset and sequences were clustered using PAM.

\newpage

\includegraphics{figure_S5.png}

**Figure S5. DESeq Variance Stabilization of OTU counts resulted in negative values that were used to calculate Bray-Curtis and Weighted UniFrac distances.** The median number of negative OTU counts that had a negative OTU count following normalization increased with the effect size and decreased as &Ntilde;~L~ increased (first row). The error bars indicate the observe 95% confidence interval. The fraction simulated datasets that had a negative OTU count following normalization also increased with effect size, but increased as &Ntilde;~L~ increased (second row). For each effect size there were 100 replicate datasets. 

\newpage

\includegraphics{figure_S6.png}

**Figure S6. Comparison of the normally distributed sequencing depths from the GlobalPatterns dataset and a log-scaled distribution of sequencing depths.** The log-scaled distribution was generated so that each sample in a simulation could have a unique number of sequences and to simulate the skew right distribution commonly seen in microbiome studies.

\newpage

\includegraphics{figure_S7.png}

**Figure S7. Rarefaction with all samples yielded clustering accuracies that were as good or better than removing the smallest 15% of samples across distance calculation methods.** This figure is analogous to Figure 8 except that the sequencing depths for each of the 80 samples in each simulation were drawn without replacement from a log-scaled distribution rather than from the GlobalPatterns sequencing depths.

\newpage

\includegraphics{figure_S8.png}

**Figure S8. Normalization and distance calculation methods vary in their sensitivity to removal of rare OTUs.** Larger values indicate that the clustering accuracy from filtered datasets were larger than those from non-filtered datasets. The median of 100 randomizations did not meaningfully vary from 0.0, but the observed 95% confidence interval varied considerably. Data are shown for a median sequencing depth (&Ntilde;~L~) of 10,000 sequences when individual sequencing depths were sampled with replacement from the GlobalPatterns dataset or without replacement from the log-scaled distribution.

\newpage

\includegraphics{figure_S9.png}

**Figure S9. Rarefaction was consistently as good or better than all other normalization methods at controlling for Type I error and maximizing power to detect differences in treatment group using alpha-diversity metrics regardless of whether sequencing depth was confounded by treatment group when using sequencing depths drawn from a lognormal distribution.** Statistical comparisons of OTU richness and Shannon diversity were performed using the non-parametric Wilcoxon two-sampled test. Type I errors were assessed as the fraction of 100 simulations that yielded a significant P-value (i.e., less than or equal to 0.05) at an effect size of 1.00. Power was assessed as the fraction of 100 simulations that yielded a significant P-Value at an effect size of 1.15. Data are shown for when the case when individual sequencing depths were sampled without replacement from the log-scaled distribution.
