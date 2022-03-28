"""
Dogukan Bayraktar
Python Version 3.9.7
Compare exon positions of transcripts matched by GFFcompare
"""
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate


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
        tcons_exon[key] = [oxford[value[0]], flair[value[1]], talon[value[2]]]
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
    Create two dictionaries. Start dict contains differences of all
    starting positions of exons. End dict contains differences of all ending
    positions of exons.

    :param compare: Dictionary with calculated differences.
    :return: Dictionaries with start and end positions.
    """
    exon_start = defaultdict(list)
    exon_end = defaultdict(list)
    for positions in compare:
        for number, position in enumerate(positions):
            exon_start[number].append(abs(position[0]))
            exon_end[number].append(abs(position[1]))
    return exon_start, exon_end


def writefile(start_dict, end_dict, outfile, header):
    """
    Writes a file with statistics on the exons.

    :param start_dict: Starting positions differences
    :param end_dict: Ending position differences
    :param outfile: Output location
    :param header: Name of comparison
    :return: File with statistics
    """
    for start_pos, end_pos in zip(start_dict.items(), end_dict.items()):
        table = [
            [header],
            ['', 'START', 'END'],
            ['EXON:', start_pos[0], end_pos[0]],
            ['Number of exons:', len(start_pos[1]), len(end_pos[1])],
            ['Average', np.average(start_pos[1]), np.average(end_pos[1])],
            ['Median:', np.median(start_pos[1]), np.median(end_pos[1])],
            ['Minimum:', np.min(start_pos[1]), np.min(end_pos[1])],
            ['Maximum:', np.max(start_pos[1]), np.max(end_pos[1])]
        ]
        with open(outfile, 'w') as file:
            file.write(tabulate(table, tablefmt="plain") + '\n\n')


def plot_boxplots(oxford_flair_start, oxford_flair_end,
                  oxford_talon_start, oxford_talon_end
                  , flair_talon_start, flair_talon_end):
    labels, data = [*zip(*oxford_flair_start.items())]
    plt.figure(figsize=(20, 10))
    plt.boxplot(data)
    plt.yscale('log')
    plt.ylabel('Number of bases')
    plt.xlabel('Exon number')
    plt.title('Exon starting position differences between flair and oxford for all transcripts')
    plt.show()

    labels, data = [*zip(*oxford_flair_end.items())]
    plt.figure(figsize=(20, 10))
    plt.boxplot(data)
    plt.yscale('log')
    plt.ylabel('Number of bases')
    plt.xlabel('Exon number')
    plt.title('Exon ending position differences between flair and oxford for all transcripts')
    plt.show()

    labels, data = [*zip(*oxford_talon_start.items())]
    plt.figure(figsize=(20, 10))
    plt.boxplot(data)
    plt.yscale('log')
    plt.ylabel('Number of bases')
    plt.xlabel('Exon number')
    plt.title('Exon starting position differences between oxford and talon for all transcripts')
    plt.show()

    labels, data = [*zip(*oxford_talon_end.items())]
    plt.figure(figsize=(20, 10))
    plt.boxplot(data)
    plt.yscale('log')
    plt.ylabel('Number of bases')
    plt.xlabel('Exon number')
    plt.title('Exon ending position differences between oxford and talon for all transcripts')
    plt.show()

    labels, data = [*zip(*flair_talon_start.items())]
    plt.figure(figsize=(20, 10))
    plt.boxplot(data)
    plt.yscale('log')
    plt.ylabel('Number of bases')
    plt.xlabel('Exon number')
    plt.title('Exon starting position differences between flair and talon for all transcripts')
    plt.show()

    labels, data = [*zip(*flair_talon_end.items())]
    plt.figure(figsize=(20, 10))
    plt.boxplot(data)
    plt.yscale('log')
    plt.ylabel('Number of bases')
    plt.xlabel('Exon number')
    plt.title('Exon ending position differences between flair and talon for all transcripts')
    plt.show()


def main():
    tcons = tracking('/home/dogukan/Downloads/exoncompare_test/gffcmp.tracking')
    oxford = gtf('/home/dogukan/Downloads/exoncompare_test/oxford_counts.gtf')
    flair = gtf('/home/dogukan/Downloads/exoncompare_test/flair_counts.gtf')
    talon = gtf('/home/dogukan/Downloads/exoncompare_test/talon_counts.gtf')

    transcript_exons = merge(tcons, oxford, flair, talon)
    oxford_flair, oxford_talon, flair_talon = calculate_difference(transcript_exons)

    oxford_flair_start, oxford_flair_end = start_end(oxford_flair)
    oxford_talon_start, oxford_talon_end = start_end(oxford_talon)
    flair_talon_start, flair_talon_end = start_end(flair_talon)

    plot_boxplots(oxford_flair_start, oxford_flair_end,
                  oxford_talon_start, oxford_talon_end,
                  flair_talon_start, flair_talon_end)

    writefile(oxford_flair_start, oxford_flair_end, '/home/dogukan/Downloads/exoncompare_test/oxford_flair.txt', 'oxford vs flair')
    writefile(oxford_talon_start, oxford_talon_end, '/home/dogukan/Downloads/exoncompare_test/oxford_talon.txt', 'oxford vs talon')
    writefile(flair_talon_start, flair_talon_end, '/home/dogukan/Downloads/exoncompare_test/flair_talon.txt', 'flair vs talon')


main()
