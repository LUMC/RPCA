#!/bin/bash
#
#SBATCH --job-name="RPCA"
#SBATCH --mem=50g
#SBATCH --ntasks-per-node=10
#SBATCH --time=8:0:0
#SBATCH --partition=all,highmem
#SBATCH --output=RPCA-%j.out

CONF="config/config.yaml"
snakemake -p --use-singularity --singularity-prefix "resources"  --singularity-args "--bind /exports" --use-conda -j 10 all --configfile $CONF
