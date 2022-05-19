#!/bin/bash
#
#SBATCH --job-name="RPCA"
#SBATCH --mem=100g
#SBATCH --ntasks-per-node=10
#SBATCH --time=24:0:0
#SBATCH --partition=all,highmem
#SBATCH --output=RPCA-%j.out
#SBATCH --mail-user="d.bayraktar@lumc.nl"
#SBATCH --mail-type=END

CONF="config/config.yaml"

# Create rulegraph
#snakemake -n --dry-run -p --use-singularity --singularity-prefix "resources"  --singularity-args "--bind /exports" --use-conda -j 10 all --configfile $CONF --rulegraph | dot -Tpdf > rulegraph.pdf

# Create Snakemake report (after creating all results)
#snakemake -p --use-singularity --singularity-prefix "resources"  --singularity-args "--bind /exports" --use-conda -j 10 all --report report.html --configfile $CONF

# Run the pipeline
snakemake -p --use-singularity --singularity-prefix "resources" --singularity-args "--bind /exports" --use-conda -j 10 all --configfile $CONF

# Dryrun and state reason for running rules
#snakemake -n --reason --dry-run -p --use-singularity --singularity-prefix "resources"  --singularity-args "--bind /exports" --use-conda -j 1 all --configfile $CONF
