#!/bin/bash

set -e

CONF="config/config.yaml"

# Run the pipeline using slurm
snakemake -p --profile slurm --use-singularity --singularity-prefix "resources" --singularity-args "--bind /exports" all --configfile $CONF
