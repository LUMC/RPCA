"""
Dogukan Bayraktar
Python Version 3.9.7
Compare exon positions of transcripts matched by GFFcompare
"""
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt


def tracking(gfftracking):
    """
    Read the data from the GFFcompare gffcmp.tracking file and
    format as a dictionary.

    :return: dictionary format: tcons_XXXX : [[transcript_id 1 ],
    [transcript_id 2 ], [transcript_id 3]]
    """
    tcons = {}
    with open(gfftracking) as file:
        for line in file:
            line = line.split()
            transcripts = line[4::]
            temp_list = []
            # '-' means that at least one pipeline did not have a matching transcript
            # only take lines with three transcripts
            if '-' not in transcripts:
                for transcript in transcripts:
                    temp_list.append(transcript.split('|')[1])
                tcons[line[0]] = temp_list
    return tcons


def gtf(file):
    """
    Read transcript data from a GTF file and format as dictionary.

    :param file: GTF file
    :return: Dictionary output: key = transcript_id, value =
    [(exon_start, exon_end), (exon_start, exon_end),
    (exon_start, exon_end) etc...]
    """
    gtf_dict = defaultdict(list)
    with open(file) as GTF:
        for line in GTF:
            line = line.split()
            if line[2] == 'exon':
                transcript_id = line[11][1:-2]
                start = int(line[3])
                end = int(line[4])
                exon = (start, end)
                gtf_dict[transcript_id].append(exon)
    return gtf_dict


def merge(tcons, oxford, flair, talon):
    """
    Place exon positions of transcripts into dictionary by matching
    transcript id

    :param tcons: tcons_XXXX : [[transcript_id 1 ], [transcript_id 2 ]]
    :param oxford: key = transcript_id, value = [(exon_start, exon_end)]
    :param flair: key = transcript_id, value = [(exon_start, exon_end)]
    :param talon: key = transcript_id, value = [(exon_start, exon_end)]
    :return: dictionary output: Key = tcons_XXXX, value = [[(exon_start,
     exon_end)], [(exon_start, exon_end)], [(exon_start, exon_end)]]
    """
    tcons_exon = {}
    for key, value in tcons.items():
        # no need to select terminal exons for single-exon transcripts
        if len(oxford[value[0]]) == 1:
            tcons_exon[key] = [oxford[value[0]], flair[value[1]], talon[value[2]]]
        # Slice terminal exons for multi-exon transcripts
        else:
            ox_start_end = [oxford[value[0]][0], oxford[value[0]][-1]]
            fl_start_end = [flair[value[1]][0], flair[value[1]][-1]]
            ta_start_end = [talon[value[2]][0], talon[value[2]][-1]]
            tcons_exon[key] = [ox_start_end, fl_start_end, ta_start_end]
    return tcons_exon


def calculate_difference(transcript_exons):
    """
    Calculate the differences between exon positions

    :param transcript_exons: Key = tcons_XXXX, value = [[(exon_start,
     exon_end)], [(exon_start, exon_end)], [(exon_start, exon_end)]]
    :return: Differences between pipeline exon positions
    """
    oxford_flair = []
    oxford_talon = []
    flair_talon = []
    for key, value in transcript_exons.items():
        value = np.array(value, dtype=object)
        oxford = value[0]
        flair = value[1]
        talon = value[2]
        a = np.subtract(oxford, flair)
        b = np.subtract(oxford, talon)
        c = np.subtract(flair, talon)
        oxford_flair.append(a)
        oxford_talon.append(b)
        flair_talon.append(c)
    return oxford_flair, oxford_talon, flair_talon


def start_end(compare):
    """
    Create four dictionaries. Start dict contains differences of all
    starting positions of exons. End dict contains differences of all
    ending positions of exons.

    :param compare: Dictionary with calculated differences.
    :return: Dictionaries with start and end positions.
    """
    start_exon_start = defaultdict(list)
    start_exon_end = defaultdict(list)
    end_exon_start = defaultdict(list)
    end_exon_end = defaultdict(list)
    single_exon_start = defaultdict(list)
    single_exon_end = defaultdict(list)
    for positions in compare:
        # single exon transcripts
        if len(positions) == 1:
            single_exon_start['0'].append(abs(positions[0][0]))
            single_exon_end['0'].append(abs(positions[0][1]))
        else:
            # multi exon transcripts
            for number, position in enumerate(positions):
                # start exon
                if number == 0:
                    start_exon_start[number].append(abs(position[0]))
                    start_exon_end[number].append(abs(position[1]))
                    # end exon
                else:
                    end_exon_start[number].append(abs(position[0]))
                    end_exon_end[number].append(abs(position[1]))
    return start_exon_start, start_exon_end, end_exon_start, end_exon_end, single_exon_start, single_exon_end


def plot(start_exon_start, start_exon_end, end_exon_start, end_exon_end, single_exon_start, single_exon_end, output):
    """
        Plots the exon differences into 3 figures with 6 subplots each.
    """
    plt.rcParams["figure.figsize"] = [15, 15]
    plt.rcParams["figure.autolayout"] = True
    plt.rcParams.update({'font.size': 15})

    fig_oxford_flair, axs = plt.subplots(nrows=3, ncols=2)
    title = output.split('/')[-1].split('.')[0]
    fig_oxford_flair.suptitle(f'Comparison of {title} start and end exons from three-way matched transcripts')

    labels_start, data_start = [*zip(*single_exon_start.items())]
    axs[0, 0].boxplot(data_start, showfliers=True)
    axs[0, 0].set(xlabel="Exon number", ylabel="Number of bases")
    axs[0, 0].set_title('Single-exon, start position differences')
    axs[0, 0].set_xticklabels(labels_start, rotation=45)

    labels_start, data_start = [*zip(*single_exon_end.items())]
    axs[0, 1].boxplot(data_start, showfliers=True)
    axs[0, 1].set(xlabel="Exon number", ylabel="Number of bases")
    axs[0, 1].set_title('Single-exon, end position differences')
    axs[0, 1].set_xticklabels(labels_start, rotation=45)

    labels_start, data_start = [*zip(*start_exon_start.items())]
    axs[1, 0].boxplot(data_start, showfliers=True)
    axs[1, 0].set(xlabel="Exon number", ylabel="Number of bases")
    axs[1, 0].set_title('Start exon, start position differences')
    axs[1, 0].set_xticklabels(labels_start, rotation=45)

    labels_start, data_start = [*zip(*start_exon_end.items())]
    axs[1, 1].boxplot(data_start, showfliers=True)
    axs[1, 1].set(xlabel="Exon number", ylabel="Number of bases")
    axs[1, 1].set_title('Start exon, end position differences')
    axs[1, 1].set_xticklabels(labels_start, rotation=45)

    labels_start, data_start = [*zip(*end_exon_start.items())]
    axs[2, 0].boxplot(data_start, showfliers=True)
    axs[2, 0].set(xlabel="Exon number", ylabel="Number of bases")
    axs[2, 0].set_title('End exon, start position differences')
    axs[2, 0].set_xticklabels(labels_start, rotation=45)

    labels_start, data_start = [*zip(*end_exon_end.items())]
    axs[2, 1].boxplot(data_start, showfliers=True)
    axs[2, 1].set(xlabel="Exon number", ylabel="Number of bases")
    axs[2, 1].set_title('End exon, end position differences')
    axs[2, 1].set_xticklabels(labels_start, rotation=45)

    plt.savefig(output, dpi=200)


def main():
    tcons = tracking(snakemake.input.tracking)
    oxford = gtf(snakemake.input.oxford)
    flair = gtf(snakemake.input.flair)
    talon = gtf(snakemake.input.talon)

    transcript_exons = merge(tcons, oxford, flair, talon)
    oxford_flair, oxford_talon, flair_talon = calculate_difference(transcript_exons)

    start_exon_start, start_exon_end, end_exon_start, end_exon_end, single_exon_start, single_exon_end = start_end(oxford_flair)
    plot(start_exon_start, start_exon_end, end_exon_start, end_exon_end, single_exon_start, single_exon_end, snakemake.output.oxford_flair)

    start_exon_start, start_exon_end, end_exon_start, end_exon_end, single_exon_start, single_exon_end = start_end(oxford_talon)
    plot(start_exon_start, start_exon_end, end_exon_start, end_exon_end, single_exon_start, single_exon_end, snakemake.output.oxford_talon)

    start_exon_start, start_exon_end, end_exon_start, end_exon_end, single_exon_start, single_exon_end = start_end(flair_talon)
    plot(start_exon_start, start_exon_end, end_exon_start, end_exon_end, single_exon_start, single_exon_end, snakemake.output.flair_talon)


main()
