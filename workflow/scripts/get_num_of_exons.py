"""
Dogukan Bayraktar
Python Version 3.10
Gets number of exons per transcript from a GTF file
"""

from collections import defaultdict
from argparse import ArgumentParser

parser = ArgumentParser(description='Gets number of exons per transcript from a GTF file')
parser.add_argument('--GTF', type=str,  help='Full path to GTF file')
parser.add_argument('--output', type=str, help='Name of the output file')
args = parser.parse_args()

if __name__ == "__main__":

    transcripts = defaultdict(int)
    with open(args.GTF, "r") as content:
        for line in content:
            if line.startswith("#"):
                continue
            line = line.split()
            feature = line[2]
            if feature == "exon":
                transcript_id = line[13][1:-2]
                transcripts[transcript_id] += 1

    exons = defaultdict(int)
    for num_exons in transcripts.values():
        exons[num_exons] += 1

    with open(args.output, "w") as outfile:
        outfile.write(f'Transcript_exon_length\tNumber_of_transcripts\n')
        sorted_exons = sorted(exons.items())
        for item in sorted_exons:
            num_exon = item[0]
            count = item[1]
            outfile.write(f'{num_exon}\t{count}\n')
