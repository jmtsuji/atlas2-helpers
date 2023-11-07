#!/usr/bin/env python
# Generate MAG table of percentages of mapped reads from ATLAS2 output
# Copyright Jackson M. Tsuji, Neufeld Research Group, 2020

# TODO - document functions

# Imports
import sys
import os
import logging
import argparse
import pandas as pd

VERSION = '0.2.0'

# GLOBAL variables (default relative paths within the ATLAS dir)
GENOME_COUNTS_FILEPATH = 'genomes/counts/counts_genomes.parquet'
ASSEMBLY_STATS_FILEPATH = 'stats/combined_contig_stats.tsv'
READ_COUNTS_FILEPATH = 'stats/read_counts.tsv'
GTDB_TAXONOMY_FILEPATH = 'genomes/taxonomy/gtdb_taxonomy.tsv'
CHECKM2_COMPLETENESS_FILEPATH = 'genomes/genome_quality.tsv'

# Set up the logger
logging.basicConfig(level=logging.INFO, format='[ %(asctime)s UTC ]: %(levelname)s: %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)


# Load the read_counts.tsv table, make custom read count stat, and get key columns
def load_read_counts(read_counts_filepath):

    # Import read totals and filter to QC only with pe and se reads
    read_counts = pd.read_csv(read_counts_filepath, sep='\t')
    read_counts = read_counts[read_counts['Step'] == 'QC']
    read_counts = read_counts[['Sample', 'Reads_pe', 'Reads_se']]

    # Re-calculate totals in a way compatible with read mapping output (i.e., total = R1 + R2 + se)
    read_counts['Reads'] = read_counts['Reads_pe'] * 2 + read_counts['Reads_se']
    read_counts = read_counts[['Sample', 'Reads']]

    return read_counts


# Load the combined_contig_stats.tsv table and get key columns
def load_assembly_stats(assembly_stats_filepath):

    # Import assembly stats and filter to assembled reads column
    assembly_stats = pd.read_csv(assembly_stats_filepath, sep='\t')
    assembly_stats = assembly_stats.rename(columns={assembly_stats.columns[0]: 'Sample',
                                           'Assembled_Reads': 'Reads'})
    assembly_stats = assembly_stats[['Sample', 'Reads']]

    return assembly_stats


# Load and bind rows of GTDBTk taxonomy tables
def load_gtdb_taxonomy_table(gtdb_classification_filepath):

    taxonomy_table = pd.read_csv(gtdb_classification_filepath, sep='\t')

    # Make all the taxonomy names title case
    taxonomy_column_names = list(taxonomy_table.drop(columns='user_genome').columns)
    taxonomy_column_names_title_case = []

    for taxonomy_column_name in taxonomy_column_names:
        taxonomy_column_names_title_case.append(taxonomy_column_name.title())

    taxonomy_column_name_mapping = dict(zip(taxonomy_column_names, taxonomy_column_names_title_case))

    taxonomy_table = taxonomy_table.rename(columns={'user_genome': 'MAG'})\
        .rename(columns=taxonomy_column_name_mapping)

    return taxonomy_table


# Load CheckM completeness data and get key columns
def load_checkm2_completeness_table(checkm2_completeness_table_filepath):

    # TODO - check if I should include the Contamination_general and Contamination_specific columns
    checkm2_completeness_table = pd.read_csv(checkm2_completeness_table_filepath,
                                             sep='\t')[['Bin Id', 'Completeness', 'Contamination']]
    checkm2_completeness_table = checkm2_completeness_table.rename(columns={'Bin Id': 'MAG'})

    return checkm2_completeness_table


# Normalize the genome count table to a particular table of total reads
# The genome count table must have only the following columns: MAG (id) and sample names
def normalize_genome_count_table(genome_count_table, totals_table):

    read_counts = totals_table.set_index('Sample')['Reads']
    genome_count_table_normalized = genome_count_table.set_index('MAG').div(read_counts, axis=1).multiply(100)

    # See the summed totals of percent mapped reads
    column_sums = genome_count_table_normalized.sum(axis=0).round(1)
    logger.info("Total normalized values:")
    for index, value in column_sums.items():
        logger.info(f'{index}: {value}')

    genome_count_table_normalized = genome_count_table_normalized.reset_index()

    return genome_count_table_normalized


def main(args):
    # Set user variables
    atlas_dir = args.atlas_dir
    output_table_filepath = args.output_file
    genome_counts_filepath = args.genome_counts
    read_counts_filepath = args.reads_total
    assembly_stats_filepath = args.reads_assembled
    gtdb_taxonomy_filepath = args.gtdb_taxonomy_table
    checkm2_completeness_filepath = args.checkm2_table

    # Auto set non-specified key flags if atlas_dir has been set
    # But keep any manual user overrides
    if atlas_dir is not False:
        if genome_counts_filepath is False:
            genome_counts_filepath = os.path.join(atlas_dir, GENOME_COUNTS_FILEPATH)
        if (read_counts_filepath is False) and (assembly_stats_filepath is False):
            # Set read_counts_filepath by default
            read_counts_filepath = os.path.join(atlas_dir, READ_COUNTS_FILEPATH)
        if gtdb_taxonomy_filepath is False:
            gtdb_taxonomy_filepath = os.path.join(atlas_dir, GTDB_TAXONOMY_FILEPATH)
        if checkm2_completeness_filepath is False:
            checkm2_completeness_filepath = os.path.join(atlas_dir, CHECKM2_COMPLETENESS_FILEPATH)

    # Startup messages
    logger.info('Running ' + os.path.basename(sys.argv[0]))
    logger.info('##### Settings #####')
    logger.info('ATLAS directory:   ' + str(atlas_dir))
    logger.info('Output filepath: ' + str(output_table_filepath))
    logger.info('Genome counts filepath: ' + str(genome_counts_filepath))
    logger.info('Read counts filepath: ' + str(read_counts_filepath))
    logger.info('Assembly stats filepath: ' + str(assembly_stats_filepath))
    logger.info('GTDB taxonomy filepath: ' + str(gtdb_taxonomy_filepath))
    logger.info('CheckM2 completeness filepath: ' + str(checkm2_completeness_filepath))
    logger.info('####################')

    # TODO - check the integrity of the input files/folders

    # Load genome counts
    genome_counts = pd.read_parquet(genome_counts_filepath)
    genome_counts = genome_counts.rename(columns={'Sample': 'MAG'})

    # Normalize the counts
    if (read_counts_filepath is not False) and (assembly_stats_filepath is not False):
        logger.error(
            'Cannot provide both -r and -R; please be clear on which normalization method you want. Exiting...')
        sys.exit(1)
    elif read_counts_filepath is not False:
        # Then normalize to total counts
        total_read_counts = load_read_counts(read_counts_filepath)
        logger.info('Normalizing MAG abundances by total read counts')
        genome_counts_normalized = normalize_genome_count_table(genome_counts, total_read_counts)
    elif assembly_stats_filepath is not False:
        # Then normalize to assembled reads
        assembly_stats = load_assembly_stats(assembly_stats_filepath)
        logger.info('Normalizing MAG abundances by assembled read counts')
        genome_counts_normalized = normalize_genome_count_table(genome_counts, assembly_stats)
    else:
        logger.error('At least one of -r and -R MUST be provided. Exiting...')
        sys.exit(1)

    # Bind on taxonomy
    if gtdb_taxonomy_filepath is not False:
        # Load and add taxonomy
        logger.info('Adding GTDB taxonomy info to the table')
        taxonomy_table = load_gtdb_taxonomy_table(gtdb_taxonomy_filepath)
        genome_counts_normalized = pd.merge(genome_counts_normalized, taxonomy_table,
                                            how='left', on='MAG', validate='1:1')

        # Sort by taxonomy
        logger.debug('Sorting table by taxonomy')
        genome_counts_normalized = genome_counts_normalized.sort_values(
            by=['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'MAG'])

    # Add on CheckM2 info
    if checkm2_completeness_filepath is not False:
        # Load and add CheckM completeness info
        logger.info('Adding CheckM completeness info to the table')
        checkm_completeness_table = load_checkm2_completeness_table(checkm2_completeness_filepath)
        genome_counts_normalized = pd.merge(genome_counts_normalized, checkm_completeness_table,
                                            how='left', on='MAG', validate='1:1')

    # Write output
    logger.info('Writing output table')
    pd.DataFrame.to_csv(genome_counts_normalized, output_table_filepath, sep='\t', index=False)

    logger.info(os.path.basename(sys.argv[0]) + ': done.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=f'Creates a TSV-format table of MAG abundances from ATLAS2. '
                    f'Copyright Jackson M. Tsuji, 2023. '
                    f'Version: {VERSION}')

    parser.add_argument('-o', '--output_file', required=True,
                        help='The path to the output TSV-format MAG table')
    parser.add_argument('-a', '--atlas_dir', required=False, default=False,
                        help='The path to the ATLAS directory. If you set this, it will set the inputs file '
                             'flags (-g, -r, -t, -c) automatically. '
                             'However, you can override with custom flags if you desire.')
    parser.add_argument('-g', '--genome_counts', required=False, default=False,
                        help='The path to the genome counts file "raw_counts_genomes.tsv" generated by ATLAS.')
    parser.add_argument('-r', '--reads_total', required=False, default=False,
                        help='The path to the total read counts file "read_counts.tsv" generated by ATLAS. '
                             'Specify this flag if you want to normalize to the total number of unassembled reads. '
                             'Cannot specify both this and -R.')
    parser.add_argument('-R', '--reads_assembled', required=False, default=False,
                        help='The path to the assembly stats file "combined_contig_stats.tsv" generated by ATLAS.'
                             'Specify this flag if you want to normalize to the number of assembled reads. '
                             'Cannot specify both this and -r.')
    parser.add_argument('-t', '--gtdb_taxonomy_table', required=False, default=False,
                        help='Paths to the parsed GTDB classification file "gtdb_taxonomy.tsv" generated by ATLAS.')
    parser.add_argument('-c', '--checkm2_table', required=False, default=False,
                        help='The path to the CheckM2 file "genome_quality.tsv" generated by ATLAS.')

    command_line_args = parser.parse_args()
    main(command_line_args)
