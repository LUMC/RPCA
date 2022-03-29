# RNA-seq Pipeline Comparison and Analyses

A Snakemake workflow for running TALON, FLAIR, and pipeline-nanopore-ref-isoforms. Performs a comparative analyses of results using tools such as GFFcompare.

# <a name="Flowchart"><a/>Flowchart

<img title="RPCA_flowchart" src="RPCA.png" alt="">

# <a name="Dependencies"><a/>Dependencies

1. Snakemake 7.3.1 <br>
[A full snakemake instalation is recommended](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#full-installation)
2. Singularity 3.7.0

# <a name="Installation"><a/>Installation

Clone the repository to desired location.

# <a name="How to run"><a/>How to run

1. Set parameters in ```config.yaml```
2. run: ```snakemake -p --use-singularity --singularity-prefix "resources"  --singularity-args "--bind *" --use-conda -j ** all```

Note * : You should provide your own directory for the --bind command so that the data is accesible from the singularity containers. <br>
Note ** : Specify number of available threads here.

# <a name="Troubleshooting"><a/>Troubleshooting

### Conda environment fails to build
Try running the workflow with an older version of snakemake such as version 5.3.2


# <a name="License"><a/>License

MIT, see LICENSE