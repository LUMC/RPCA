"""
Dogukan Bayraktar
Python Version 3.9.7
Adds transcript counts to GTF files for FLAIR, TALON, and OXFORD.
"""
import logging


def format_oxford(abundance):
    """
    Takes oxford abundance file and produces dictionary with counts per
    transcript.

    :param abundance: File containing transcript counts
    :return: Dictionary: key = transcript id, value = counts
    """
    ox_dict = {}
    with open(abundance, 'r') as file:
        # Skip header
        next(file)
        for line in file:
            line = line.strip("\n").split(',')
            transcript_id = line[0]
            # Can contain an empty string instead of 0
            counts = sum([int(count) for count in line[1:] if count != ''])
            ox_dict[transcript_id] = counts
    return ox_dict


def format_talon(abundance):
    """
    Takes talon abundance file and produces dictionary with counts per
    transcript.

    :param abundance: File containing transcript counts
    :return: Dictionary: key = transcript id, value = counts
    """
    talon_dict = {}
    with open(abundance, 'r') as file:
        # Skip header
        next(file)
        for line in file:
            line = line.strip('\n').split('\t')
            transcript_id = line[3]
            # Always has a round number
            counts = sum(map(int, line[11:]))
            talon_dict[transcript_id] = counts
    return talon_dict


def format_flair(abundance):
    """
    Takes flair abundance file and produces dictionary with counts per
    transcript.

    :param abundance: File containing transcript counts
    :return: Dictionary: key = transcript id, value = counts
    """
    flair_dict = {}
    with open(abundance, 'r') as file:
        # Skip header
        next(file)
        for line in file:
            line = line.strip('\n').split('\t')
            transcript_id = line[0]
            # Numbers are always rounded but saves them as x.0
            counts = sum(map(float, line[1:]))
            flair_dict[transcript_id] = counts
    return flair_dict


def combine(gtf_file, counts, output):
    """
    Loops through the GTF file and checks the counts dictionary
    to add counts to the line and then writes it to a new file.

    Flair combines the transcript id with the gene id for known
    transcripts using '_' as a delimiter. This can cause problems
    when the names in the annotation file contain '_'.

    Flair and Oxford can have transcripts present in the GTF file
    that are missing from the abundance file. This is caused by minimap2
    mapping ambiguity in the case of flair, and the long-read mode from
    stringtie for oxford. If there are any missing transcripts these
    will be written to the missing transcript log.

    :param gtf_file: The GTF file produced by a pipeline.
    :param counts: Dictionary: key = transcript id, value = counts
    :param output: Output directory and file name
    :return: Writes a file to the output param
    """
    # Setup logging
    logging.basicConfig(filename=output + '_missing_transcripts.log',
                        level=logging.DEBUG,
                        format='%(asctime)s %(levelname)s %(name)s %(message)s')
    logger = logging.getLogger(__name__)

    # Read the GTF file
    with open(gtf_file) as GTF, open(output, 'w') as outfile:
        for line in GTF:
            # Check if line is a transcript
            if line.split()[2] == 'transcript':
                try:
                    # Check for transcript id in counts
                    transcript_id = line.split()[11][1:-2]
                    count = counts[transcript_id]
                except KeyError:
                    # Transcript is missing from the abundance file
                    try:
                        # This is for the flair '_' delimiter issue
                        gene_id = line.split()[9][1:-2]
                        new_transcript_id = transcript_id + '_' + gene_id
                        count = counts[new_transcript_id]
                    except KeyError as err:
                        # transcript is actually missing
                        # Add missing transcript id to missing transcripts.log
                        logger.error(err)
                        count = 0
                # Add transcript count to end of line and write to file
                # counts are not in TPM, hijacking the column for use with GFFcompare
                outfile.write(line[:-1] + ' ' + 'TPM' + ' ' + '"' + str(count) + '"' + ';' + '\n')
            else:
                # Line is not a transcript
                outfile.write(line)


def main():
    """
    Takes abundance file and passes it to a format function. The
    resulting dictionary is passed to combine. Combine produces a new
    GTF that has counts for transcripts.
    :return: -
    """
    flair_abundance = snakemake.input.flair_count
    flair_gtf = snakemake.input.flair_transcripts
    flair_dict = format_flair(flair_abundance)
    combine(flair_gtf, flair_dict, snakemake.output.flair_combined)

    ox_abundance = snakemake.input.oxford_count
    ox_gtf = snakemake.input.oxford_transcript
    ox_dict = format_oxford(ox_abundance)
    combine(ox_gtf, ox_dict, snakemake.output.oxford_combined)

    talon_abundance = snakemake.input.talon_count
    talon_gtf = snakemake.input.talon_transcripts
    talon_dict = format_talon(talon_abundance)
    combine(talon_gtf, talon_dict, snakemake.output.talon_combined)


main()
