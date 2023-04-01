#!/bin/bash
#SBATCH --job-name=run_wnwn
#SBATCH --output=%x-%j.out
#SBATCH --mail-user=pschloss@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --partition=standard
#SBATCH --account=pschloss1
#SBATCH --time=300:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8g



### Conda environment:
source ~/miniconda3/etc/profile.d/conda.sh
conda activate nr-base

### Run script
#time snakemake -nr
#time snakemake --use-conda --conda-frontend mamba -nr data/sim_a/1_1000_1.filter.RDS
#time snakemake --profile config/slurm --batch all_targets=9/10
time snakemake --profile config/slurm
 