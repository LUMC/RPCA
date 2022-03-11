import logging

# Setup argparse
# import argparse
# parser = argparse.ArgumentParser()
# parser.add_argument('gtf', help='Path to GTF file', type=str)
# parser.add_argument('abundance', help='Path to file containing transcript_id and counts', type=str)
# parser.add_argument('output', help='Name of combined output GTF', type=str)
# args = parser.parse_args()


def read_abundance(abundace_file):
    # Read the abundance file, create and populate counts dictionary
    counts = {}
    with open(abundace_file) as file:
        for line in file:
            line = line.split()
            # line[0] = transcript id
            # line[1] = counts
            counts[line[0]] = line[1]
    return counts


def combine(gtf_file, counts, output):
    # Setup logging
    logging.basicConfig(filename=output + '_missing_transcripts.log', level=logging.DEBUG,
                        format='%(asctime)s %(levelname)s %(name)s %(message)s')
    logger = logging.getLogger(__name__)

    # Read the GTF file
    with open(gtf_file) as GTF, open(output, 'w') as outfile:
        for line in GTF:
            # Check if line is a transcript
            if line.split()[2] == 'transcript':
                try:
                    # Search for transcript id in counts
                    count = counts[line.split()[11][1:-2]]
                except KeyError as err:
                    # Transcript is missing from the abundance file
                    count = '0'
                    # Add missing transcript id to missing_transcripts.log
                    logger.error(err)
                # Add transcript count to end of line and write to file
                outfile.write(line[:-1] + ' ' + 'TPM' + ' ' + '"' + count + '"' + ';' + '\n')
            else:
                # Line is not a transcript
                outfile.write(line)


def main():
    abundance_files = [snakemake.input.oxford_count,
                       snakemake.input.flair_count,
                       snakemake.input.talon_count]
    transcript_files = [snakemake.input.oxford_transcripts,
                        snakemake.input.flair_transcripts,
                        snakemake.input.talon_transcripts]
    output_files = [snakemake.output.oxford_combined,
                    snakemake.output.flair_combined,
                    snakemake.output.talon_combined]

    for abundance, transcript, output in zip(abundance_files, transcript_files, output_files):
        counts_dict = read_abundance(abundance)
        combine(transcript, counts_dict, output)
        counts_dict.clear()


main()
