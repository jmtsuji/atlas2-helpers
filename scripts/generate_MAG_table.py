#!/usr/bin/env python
# Generate MAG table of percentages of mapped reads from ATLAS2 output
# Copyright Jackson M. Tsuji, Neufeld Research Group, 2019

# Imports
import sys
import os
import time
import logging
import argparse

import pandas as pd
import re

# GLOBAL variables (defaults)
GENOME_COUNTS_FILEPATH = "genomes/counts/raw_counts_genomes.tsv"
READ_COUNTS_FILEPATH = "stats/read_counts.tsv"
ASSEMBLY_STATS_FILEPATH = "stats/combined_contig_stats.tsv"
GTDBTK_TAXONOMY_FILEPATHS = ('genomes/taxonomy_gtdbtk/gtdbtk.bac120.summary.tsv',
                                   'genomes/taxonomy_gtdbtk/gtdbtk.ar122.summary.tsv')
CHECKM_COMPLETENESS_FILEPATH = 'genomes/checkm/completeness.tsv'
OUTPUT_TABLE_FILEPATH = "../test_table.tsv"

# Set up the logger
logging.basicConfig(level=logging.INFO, format='[ %(asctime)s UTC ]: %(levelname)s: %(module)s: %(message)s')
logging.Formatter.converter = time.gmtime
logger = logging.getLogger(__name__)

def load_read_counts(read_counts_filepath):
    # Import read totals and filter to QC only with pe and se reads
    read_counts = pd.read_csv(read_counts_filepath, sep = '\t', header = 0)
    read_counts = read_counts[read_counts['Step'] == 'QC']
    read_counts = read_counts[['Sample', 'Reads_pe', 'Reads_se']]

    # Re-calculate totals in a way compatible with read mapping output (i.e., total = R1 + R2 + se)
    read_counts['Reads_total'] = read_counts['Reads_pe'] * 2 + read_counts['Reads_se']
    read_counts = read_counts[['Sample', 'Reads_total']]

    return(read_counts)

def load_assembly_stats(assembly_stats_filepath):
    # Import assembly stats and filter to assembled reads column
    assembly_stats = pd.read_csv(assembly_stats_filepath, sep = '\t', header = 0)
    assembly_stats.rename(columns = {assembly_stats.columns[0]: 'Sample'}, inplace = True)
    assembly_stats = assembly_stats[['Sample', 'Assembled_Reads']]

    return(assembly_stats)

# Loads and binds rows of GTDBTk taxonomy tables
def load_gtdbtk_taxonomy_table(gtdbtk_classification_filepaths):
    taxonomy_table_list = []
    
    for gtdbtk_classification_filepath in gtdbtk_classification_filepaths:
        taxonomy_table_single = pd.read_csv(gtdbtk_classification_filepath, sep = '\t', header = 0)
        taxonomy_table_single.rename(columns = {'user_genome': 'MAG_ID'}, inplace = True)
        taxonomy_table_single = taxonomy_table_single[['MAG_ID', 'classification']]
        taxonomy_table_list.append(taxonomy_table_single)
        
    taxonomy_table = pd.concat(taxonomy_table_list)
    return(taxonomy_table)

def parse_gtdbtk_taxonomy_entry(taxonomy_entry, resolve = True):
    # EXAMPLE: 'd__Bacteria;p__Actinobacteriota;c__Acidimicrobiia;o__Microtrichales;f__;g__;s__'
    
    taxonomy_split = str(taxonomy_entry).split(sep = ';')
    if len(taxonomy_split) != 7:
        print('Taxonomy entry is ' + str(len(taxonomy_split)) + ' long, not 7 as expected. Exiting...')
        print('Entry was: ' + str(taxonomy_entry))
        # sys.exit(1)

    # Remove header pieces
    # TODO - confirm they are in the right order (d, p, c, o, f, g, s)
    taxonomy_split = [ re.sub("[dpcofgs]__", "", level) for level in taxonomy_split ]

    # Fill in empty parts, if they exist
    if '' in taxonomy_split and resolve is True:
        # print('Resolving ambiguities')

        # Get the positions of the empty spots
        empty_taxa = []
        for taxonomy_level in taxonomy_split:
            if taxonomy_level == '':
                empty_taxa.append(True)
            else:
                empty_taxa.append(False)

        # Get the numeric index of the first empty taxon
        # See https://stackoverflow.com/a/9542768, accessed Sept. 18, 2019
        first_empty_taxon = empty_taxa.index(True)

        if False in empty_taxa[first_empty_taxon:]:
            print('There seems to be an empty entry in the middle of your taxonomy levels. Cannot resolve. Exiting...')
            print('Entry was: ' + str(taxonomy_entry))
            # sys.exit(1)

        filler_entry = 'Unresolved_' + taxonomy_split[(first_empty_taxon-1)]

        for taxonomy_level_index in range(first_empty_taxon, 7):
            taxonomy_split[taxonomy_level_index] = filler_entry

    return(taxonomy_split)

def load_and_parse_gtdbtk_taxonomy_table(gtdbtk_classification_filepaths, resolve = True):
    # Load the table
    taxonomy_table = load_gtdbtk_taxonomy_table(gtdbtk_classification_filepaths)
    
    # Parse all taxonomy entries to a map object, then convert to DataFrame
    # Pass additional argments to map via lambad function - see https://stackoverflow.com/a/31995624 (accessed Sept. 18, 2019)
    # TODO - is the lambda function the best option?
    # TODO - make multithreaded using pool: https://stackoverflow.com/a/23537302
    taxonomy_entries_parsed = map(lambda entry: parse_gtdbtk_taxonomy_entry(entry, resolve = resolve), 
                                  taxonomy_table['classification'].tolist())
    taxonomy_table_parsed = pd.DataFrame(taxonomy_entries_parsed, 
                                  columns = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'])

    # Add MAG ID column
    taxonomy_table_parsed.insert(loc = 0, column = 'MAG_ID', value = taxonomy_table['MAG_ID'].tolist())
    
    return(taxonomy_table_parsed)


def load_checkm_completeness_table(checkm_completeness_table_filepath):

    # Load CheckM data
    checkm_completeness_table = pd.read_csv(checkm_completeness_table_filepath, sep = '\t', header = 0)
    checkm_completeness_table.rename(columns = {'Bin Id': 'MAG_ID',
                                               'Strain heterogeneity': 'Strain_heterogeneity'},
                                     inplace = True)
    checkm_completeness_table = checkm_completeness_table[['MAG_ID', 'Completeness', 'Contamination', 'Strain_heterogeneity']]

    return(checkm_completeness_table)

def normalize_genome_count_table(genome_count_table, totals_table):
    normalized_totals = []
    
    for row_tuple in totals_table.itertuples(index = False):
        sample_ID = row_tuple[0]
        total_reads = row_tuple[1]

        #print(sample_ID + ': ' + str(total_reads))

        normalized_totals.append(genome_count_table[sample_ID] / total_reads * 100)
        #count_table_normalized[sample_ID] = count_table[sample_ID] / total_reads * 100
        
    # Bind columns and add MAG IDs
    genome_count_table_normalized = pd.concat(normalized_totals, axis = 1)
    genome_count_table_normalized.insert(loc = 0, column = 'MAG_ID',
                                         value = genome_count_table['MAG_ID'].tolist())
    
    # See the summed totals of percent mapped reads
    column_sums = genome_count_table_normalized.drop(columns = 'MAG_ID').apply(sum, axis = 0)
    print("Column sums:")
    print(column_sums)
    
    return(genome_count_table_normalized)

def main(args):
    # Set user variables
    # TODO - set default if some kind of ATLAS dir variable is set to True
    # TODO add variable normalization_method

    genome_counts_filepath = args.genome_counts
    read_counts_filepath = args.read_counts
    assembly_stats_filepath = args.assembly_stats
    gtdbtk_taxonomy_filepaths = args.taxonomy_tables
    checkm_completeness_filepath = args.checkm_table
    output_table_filepath = args.output

    # Set sort_features to True if rename_features is True
    if rename_features is True:
        sort_features = True

    # Startup messages
    logger.info('Running ' + os.path.basename(sys.argv[0]))
    logger.info('Feature table filepath: ' + feature_table_filepath)
    logger.info('Representative sequences filepath: ' + str(rep_seq_filepath))
    logger.info('Taxonomy filepath: ' + str(taxonomy_filepath))
    logger.info('Feature ID colname: ' + str(feature_id_colname))
    logger.info('Sort Feature IDs roughly by relative abundance?: ' + str(sort_features))
    logger.info('Rename Feature IDs sequentially?: ' + str(rename_features))

    # Load genome counts
    genome_counts = pd.read_csv(genome_counts_filepath, sep='\t', header=0)
    genome_counts.rename(columns={'Sample': 'MAG_ID'}, inplace=True)

    # Load supplementary data
    total_read_counts = load_read_counts(read_counts_filepath)
    assembly_stats = load_assembly_stats(assembly_stats_filepath)
    taxonomy_table = load_and_parse_gtdbtk_taxonomy_table(gtdbtk_taxonomy_filepaths, resolve = True)
    checkm_completeness_table = load_checkm_completeness_table(checkm_completeness_filepath)

    # Normalize the counts
    if normalization_method == 'total':
        genome_counts_normalized = normalize_genome_count_table(genome_counts, total_read_counts)
    elif normalization_method == 'assembled':
        genome_counts_normalized = normalize_genome_count_table(genome_counts, assembly_stats)
    else:
        # TODO - make this statement more informative
        logger.error('Selected normalization method is not valid. Exiting...')
        sys.exit(1)

    # Bind on other info
    genome_counts_normalized = pd.merge(genome_counts_normalized, taxonomy_table,
                                        how = 'left', on = 'MAG_ID',
                                        sort = False, validate = 'one_to_one')
    genome_counts_normalized = pd.merge(genome_counts_normalized, checkm_completeness_table,
                                        how = 'left', on = 'MAG_ID',
                                        sort = False, validate = 'one_to_one')

    # Sort by taxonomy
    genome_counts_normalized = genome_counts_normalized.sort_values(
        by = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'MAG_ID'])

    # Write output
    logger.info('Writing output table')
    pd.DataFrame.to_csv(genome_counts_normalized, output_table_filepath, sep = '\t', index = False)

    logger.info('Done.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Creates a TSV-format table of MAG abundances from ATLAS2.'
                    'Copyright Jackson M. Tsuji, Neufeld Research Group, 2019.')

    # TODO - make these parsers correct
    parser.add_argument('-f', '--feature_table', required=True,
                        help='The path to the input TSV feature table file.')
    parser.add_argument('-s', '--rep_seqs', required=False, default=False,
                        help='The path to the input FastA representative sequences file. Sequences will be added as the '
                             '"ReprSequences" column. You can optionally omit this flag and not have sequences added to the table.')
    parser.add_argument('-t', '--taxonomy', required=False, default=False,
                        help='The path to the input taxonomy file. Taxonomy will be added as the "Consensus.Lineage" column. '
                             'You can optionally omit this flag and not have taxonomy added to the table.')
    parser.add_argument('-o', '--output_feature_table', required=True,
                        help='The path to the output TSV feature table.')
    parser.add_argument('-N', '--feature_id_colname', required=False, default='Feature ID',
                        help='The name of the first column of the output ASV table. [Default: "Feature ID"]')
    parser.add_argument('-S', '--sort_features', required=False, action='store_true',
                        help='Optionally sort Feature IDs roughly based on overall abundance.')
    parser.add_argument('-R', '--rename_features', required=False, action='store_true',
                        help='Optionally rename the Feature IDs sequentially, roughly based on overall abundance. '
                             'Automatically sets --sort_features')
    # TODO - add optional flag to parse taxonomy into 7 ranks
    # TODO - add option to auto-detect if a QZA file is provided instead of the unpackaged file. Deal with the converstions. Same for if a BIOM file is provided.

    command_line_args = parser.parse_args()
    main(command_line_args)