#!/bin/bash
#
#SBATCH --job-name="RPCA"
#SBATCH --mem=30g
#SBATCH --ntasks-per-node=4
#SBATCH --time=1:0:0
#SBATCH --partition=all,highmem
#SBATCH --output=RPCA-%j.out

CONF="config/config.yaml"
snakemake -p --use-singularity --singularity-prefix "resources"  --singularity-args "--bind /exports" --use-conda -j 4 all --configfile $CONF
