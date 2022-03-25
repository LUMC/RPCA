"""
Dogukan Bayraktar
Python Version 3.9.7
Compare exon positions of transcripts matched by GFFcompare
"""
from collections import defaultdict
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate


def tracking(gfftracking):
    """
    Function to read data from the GFFcompare gffcmp.tracking file and
    format as dictionary.
    :return: dictionary format: tcons_XXXX : [[transcript_id 1 ],
    [transcript_id 2 ], [transcript_id 3]]
    """
    tcons = {}
    with open(gfftracking) as file:
        for line in file:
            line = line.split()
            transcripts = line[4::]
            temp_list = []
            # only take lines with three transcripts
            if '-' not in transcripts:
                for transcript in transcripts:
                    temp_list.append(transcript.split('|')[1])
                tcons[line[0]] = temp_list
    return tcons


def gtf(file):
    """
    Function to read transcript data from a GTF file.
    :param file: GTF file
    :return: dictonary output: transcript_id: [(exon_start, exon_end),
    (exon_start, exon_end), (exon_start, exon_end) etc...]
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
    Merges output from the two functions above.
    :param tcons: Dictonary with transcript id's as values
    :param oxford: oxford exon positions per transcript
    :param flair: flair exon positions per transcript
    :param talon: talon exon positions per transcript
    :return: dictonary output: tcons_XXXX : [[(exon_start, exon_end)],
    [(exon_start, exon_end)], [(exon_start, exon_end)]]
    """
    tcons_exon = {}
    for key, value in tcons.items():
        tcons_exon[key] = [oxford[value[0]], flair[value[1]], talon[value[2]]]
    return tcons_exon


def calculate_difference(transcript_exons):
    """
    Checks the differences between exon positions
    :param transcript_exons:
    :return: differences between pipeline exon positions
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
    Split list into start and end dictonaries
    :param compare:
    :return:
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
    Function writes a file with statistics on the exons.
    :param start_dict:
    :param end_dict:
    :param outfile:
    :param header:
    :return:
    """
    for item1, item2 in zip(start_dict.items(), end_dict.items()):
        table = [
            [header],
            ['', 'START', 'END'],
            ['EXON:', item1[0], item2[0]],
            ['Number of exons:', len(item1[1]), len(item2[1])],
            ['Average', np.average(item1[1]), np.average(item2[1])],
            ['Median:', np.median(item1[1]), np.median(item2[1])],
            ['Minimum:', np.min(item1[1]), np.min(item2[1])],
            ['Maximum:', np.max(item1[1]), np.max(item2[1])]
        ]
        with open(outfile, 'a') as file:
            file.write(tabulate(table, tablefmt="plain") + '\n\n')


def main():
    tcons = tracking('')
    oxford = gtf('')
    flair = gtf('')
    talon = gtf('')

    transcript_exons = merge(tcons, oxford, flair, talon)
    oxford_flair, oxford_talon, flair_talon = calculate_difference(transcript_exons)

    oxford_flair_start, oxford_flair_end = start_end(oxford_flair)
    oxford_talon_start, oxford_talon_end = start_end(oxford_talon)
    flair_talon_start, flair_talon_end = start_end(flair_talon)

    writefile(oxford_flair_start, oxford_flair_end, 'oxford_flair.txt', 'oxford vs flair')
    writefile(oxford_talon_start, oxford_talon_end, 'oxford_talon.txt', 'oxford vs talon')
    writefile(flair_talon_start, flair_talon_end, 'flair_talon.txt', 'flair vs talon')


main()
