# Re-evaluating *Waste not, want not...*

This is a reanalaysis of McMurdie and Holmes's 2014 manuscript, ["Waste Not, Want Not: Why Rarefying Microbiome Data Is 
Inadmissible"](https://doi.org/10.1371/journal.pcbi.1003531) published in PLOS Computational Biology. The code, primarily
R markdown documents, is contained in the directory downloaded from the journal's website as
[Protocol S1](https://doi.org/10.1371/journal.pcbi.1003531.s001).

## Running code

Make sure conda and mamba are installed. Installing TinyTex may take some fenagling and on Mac OSX seems to need to be installed in the home directory Library (i.e., `~/Library/TinyTex`) using `install.packages`. See [this issue](https://github.com/rstudio/tinytex/issues/24) for more clues on installing `tinytex`.

```bash
conda config --set channel_priority strict
mamba env create -f workflow/envs/nr-base.yml 
conda activate nr-base
```

Use Snakemake (installed in `nr-base` environment) to build project

```bash
snakemake --use-conda --conda-frontend mamba -c1 write_paper
```
