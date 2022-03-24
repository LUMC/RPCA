from collections import defaultdict
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate


def tracking():
    """
    Function to read data from the GFFcompare gffcmp.tracking file.
    dictionary format: tcons_XXXX : [[transcript_id 1 ], [transcript_id 2 ], [transcript_id 3]]
    :return:
    """
    tcons = {}
    with open('/home/dogukan/Documents/afstuderen/notebook/selas/selas_new/selas_new_gffcmp.tracking') as file:
        for line in file:
            line = line.split()
            transcripts = line[4::]
            temp_list = []
            for transcript in transcripts:
                if transcript == '-':
                    temp_list.append('-')
                else:
                    temp_list.append(transcript.split('|')[1])
            tcons[line[0]] = temp_list
    return tcons


def gtf(file):
    """
    Function to read transcript data from a GTF file.
    dictonary output: transcript_id: [(exon_start, exon_end), (exon_start, exon_end), (exon_start, exon_end) etc...]
    :param file:
    :return:
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
    dictonary output: tcons_XXXX : [[(exon_start, exon_end)], [(exon_start, exon_end)], [(exon_start, exon_end)]]
    :param tcons:
    :param oxford:
    :param flair:
    :param talon:
    :return:
    """
    tcons_exon = {}
    for key, value in tcons.items():
        tcons_exon[key] = [oxford[value[0]], flair[value[1]], talon[value[2]]]
    return tcons_exon


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
    pass


main()
