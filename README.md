# RNA-seq Pipeline Comparison and Analyses

A Snakemake workflow for running TALON, FLAIR, and pipeline-nanopore-ref-isoforms. Performs a comparative analyses of results using tools such as GFFcompare using GFFcompare.

# <a name="Flowchart"><a/>Flowchart

<img title="RPCA_flowchart" src="RPCA.png" alt="">

# <a name="Dependencies"><a/>Dependencies

1. Snakemake 7.3.1
2. Singularity 3.7.0

### Installation

Clone the repository to desired location.

# <a name="How to run"><a/>How to run

1. Set parameters in ```config.yaml```
2. run: ```snakemake -p --use-singularity --singularity-prefix "resources"  --singularity-args "--bind *" --use-conda -j ** all```

Note * : You should provide your own directory for the --bind command so that the data is accesible from the singularity containers.
Note ** : Specify number of available threads here.

# <a name="License"><a/>License

MIT, see LICENSE