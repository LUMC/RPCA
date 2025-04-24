FROM condaforge/miniforge3:24.11.3-2

RUN mamba install -c conda-forge -c bioconda -y \
        python \
        pysam \
        tqdm \
        minimap2 \
        samtools \
        seqkit \
        jupyter=1.0.0 \
        numpy \
        pandas \
        matplotlib \
        matplotlib-venn \
        upsetplot \
        plotly && \
    pip install df2img && \
    pip cache purge && \
    mamba clean -a -y