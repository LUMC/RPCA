# RNA-seq Pipeline Comparison and Analyses

A Snakemake workflow for running TALON, FLAIR, and pipeline-nanopore-ref-isoforms. Performs a comparative analyses of results using tools such as GFFcompare using GFFcompare.

### Flowchart

<img title="RPCA_flowchart" src="RPCA.png" alt="">

### Dependencies

1. snakemake 7.3.1
2. Singularity 3.7.0

### Installation

Clone the repository to desired location.

### How to run

1. Set parameters in the config.yaml
2. run: snakemake -p --use-singularity --singularity-prefix "resources"  --singularity-args "--bind /exports" --use-conda -j 10 all

### License

MIT, see LICENSE