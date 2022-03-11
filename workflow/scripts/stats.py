import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

# this will probably give a error, not sure how it works when input is multiple from wildcards
files = snakemake.input.alignstats
# dictionary for dataframe
stats_dict = defaultdict(list)
# populate dictionary with alignment stats
for file in files:
    with open(file) as file:
        for line in file:
            line = line.strip('\n').split('\t')
            stats_dict[line[0]].append(float(line[1]))
# convert to dataframe
stats = pd.DataFrame.from_dict(stats_dict)

# plot read alignments
cols = ['reads mapped:', 'reads unmapped:']
stats[cols].plot(kind='bar', figsize=(10, 5))
plt.yscale('log')
plt.ylabel('Number of reads')
plt.xlabel('sample number')
plt.title('Read alignment per sample')
plt.legend(loc='upper left')
plt.savefig(snakemake.output.readalign, dpi=200)

# plot mismatches
stats['mismatches:'].plot(kind='bar', figsize=(10, 5))
plt.yscale('log')
plt.ylabel('Number of bases')
plt.xlabel('sample number')
plt.title('Mismatches per sample (Adenovirus)')
plt.savefig(snakemake.output.mismatch, dpi=200)

# plot average read length
stats['average length:'].plot(kind='bar', figsize=(10, 5))
plt.ylabel('Number of bases')
plt.xlabel('sample number')
plt.title('Average read length per sample')
plt.savefig(snakemake.output.readlen, dpi=200)


# dictionary for dataframe
subreads_counts = {}
# populate dictionary
with open(snakemake.input.annocount, 'r') as file:
    next(file)
    next(file)
    for line in file:
        line = line.split()
        a = [int(line[6])]
        subreads_counts[line[0]] = sum(a)

cols = ['chr', 'start', 'end', 'strand', 'length', 'T0', 'T1', 'T2',
        'T3', 'T4', 'T5']
# convert dictionary to dataframe
subreads_df = pd.DataFrame.from_dict(subreads_counts,
                                     orient='index', columns=cols)
# plot total count per sample
sample_totals = [subreads_df['T0'].sum(),
                 subreads_df['T1'].sum(),
                 subreads_df['T2'].sum(),
                 subreads_df['T3'].sum(),
                 subreads_df['T4'].sum(),
                 subreads_df['T5'].sum()]
plt.bar(['T0','T1','T2','T3','T4','T5'],sample_totals)
plt.ylabel('Number of transcripts')
plt.xlabel('sample number')
plt.title('Total count per sample')
plt.savefig(snakemake.output.anno_count, dpi=200)

# plot count per annotated gene
labels, data = [*zip(*subreads_counts.items())]
plt.bar(labels,data)
plt.ylabel('Number of raw counts')
plt.xlabel('feature')
plt.title('Raw counts per feature for sample 5')
plt.xticks(rotation=30, ha='right')
plt.savefig(snakemake.output.anno_count_all, dpi=200)

# dictionary for data
sample_dict = defaultdict(list)
# read in data
talon_samples = np.genfromtxt(snakemake.input.talon_count,
                              delimiter='\t', skip_header=True)
flair_samples = np.genfromtxt(snakemake.input.flair_count,
                              delimiter='\t', skip_header=True)
oxford_samples = np.genfromtxt(snakemake.input.oxford_count,
                               delimiter=',', skip_header=True)
subreads_samples = np.genfromtxt(snakemake.input.annocount,
                                 delimiter='\t', skip_header=True)
# populate dict
for i, j in enumerate(range(1, len(oxford_samples[0]))):
    sample_dict[i + 1].append(np.nansum(oxford_samples[:, j]))
for i, j in enumerate(range(1, len(flair_samples[0]))):
    sample_dict[i + 1].append(np.nansum(flair_samples[:, j]))
for i, j in enumerate(range(11, len(talon_samples[0]))):
    sample_dict[i + 1].append(np.nansum(talon_samples[:, j]))
for i, j in enumerate(range(7, len(subreads_samples[0]))):
    sample_dict[i + 1].append(np.nansum(subreads_samples[:, j]))

samples = pd.DataFrame.from_dict(sample_dict,
                                 orient='index',
                                 columns=['oxford', 'flair', 'talon',
                                          'subread'])

# plot total number of counts per sample per pipeline
samples.plot(kind="bar", figsize=(10, 5))
plt.yscale('log')
plt.ylabel('Raw total count')
plt.xlabel('sample number')
plt.title('Total number of counts per sample per pipeline')
plt.savefig(snakemake.output.count_comp, dpi=200)

