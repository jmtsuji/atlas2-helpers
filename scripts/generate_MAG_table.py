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
# TODO - change from tuple to list
GTDBTK_TAXONOMY_FILEPATHS = ('genomes/taxonomy_gtdbtk/gtdbtk.bac120.summary.tsv',
                                   'genomes/taxonomy_gtdbtk/gtdbtk.ar122.summary.tsv')
CAT_TAXONOMY_FILEPATH = 'genomes/taxonomy/taxonomy_names.tsv'
CHECKM_COMPLETENESS_FILEPATH = 'genomes/checkm/completeness.tsv'

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
# TODO - change from tuple to list
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

def load_CAT_taxonomy_table(CAT_taxonomy_table_filepath):
    # Load table
    taxonomy_table = pd.read_csv(CAT_taxonomy_table_filepath, sep = '\t', header = 0)
    taxonomy_table.rename(columns = {'# bin': 'MAG_ID',
                                     'superkingdom': 'Domain',
                                     'phylum': 'Phylum',
                                     'class': 'Class',
                                     'order': 'Order',
                                     'family': 'Family',
                                     'genus': 'Genus',
                                     'species': 'Species'},
                          inplace = True)
    taxonomy_table = taxonomy_table[['MAG_ID', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']]

    # Example
    ## bin	superkingdom	phylum	class	order	family	genus	species
    #MAG001    Bacteria: 1.00    Nitrospirae: 0.95    Nitrospira: 0.91    Nitrospirales: 0.91    Nitrospiraceae: 0.91    Nitrospira: 0.88    not classified

    # TODO:
    ## - Get rid of the ': [score]' in the taxonomy entry
    ## - Resolve 'not classified' as done in the GTDB taxonomy

    return(taxonomy_table)

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
    atlas_dir = args.atlas_dir
    output_table_filepath = args.output_file

    # Auto set other filepaths if atlas_dir has been set
    if atlas_dir is False:
        genome_counts_filepath = args.genome_counts
        read_counts_filepath = args.reads_total
        assembly_stats_filepath = args.reads_assembled
        CAT_taxonomy_filepath = args.cat_taxonomy_table
        gtdbtk_taxonomy_filepaths = args.gtdb_taxonomy_tables
        checkm_completeness_filepath = args.checkm_table
    else:
        genome_counts_filepath = os.path.join(atlas_dir, GENOME_COUNTS_FILEPATH)
        read_counts_filepath = os.path.join(atlas_dir, READ_COUNTS_FILEPATH)
        # TODO - add support for assembly_stats_filepath
        assembly_stats_filepath = False
        # assembly_stats_filepath = os.path.join(atlas_dir, ASSEMBLY_STATS_FILEPATH)
        gtdbtk_taxonomy_filepaths = []
        for entry in GTDBTK_TAXONOMY_FILEPATHS:
            gtdbtk_taxonomy_filepaths.append(os.path.join(atlas_dir, entry))
        # TODO - add support for CAT taxonomy filepath
        CAT_taxonomy_filepath = False
        # CAT_taxonomy_filepath = os.path.join(atlas_dir, CAT_TAXONOMY_FILEPATH)
        checkm_completeness_filepath = os.path.join(atlas_dir, CHECKM_COMPLETENESS_FILEPATH)

    # Startup messages
    logger.info('Running ' + os.path.basename(sys.argv[0]))
    logger.info('ATLAS directory: ' + str(atlas_dir))
    logger.info('Output filepath: ' + str(output_table_filepath))
    logger.info('Genome counts filepath: ' + str(genome_counts_filepath))
    logger.info('Read counts filepath: ' + str(read_counts_filepath))
    logger.info('Assembly stats filepath: ' + str(assembly_stats_filepath))
    # TODO - somehow print out the list
    logger.info('GTDB taxonomy filepaths: ' + ', '.join(list(gtdbtk_taxonomy_filepaths)))
    logger.info('CAT taxonomy filepath: ' + str(CAT_taxonomy_filepath))
    logger.info('CheckM completeness filepath: ' + str(checkm_completeness_filepath))

    # TODO - check the integrity of the input files/folders

    # Load genome counts
    genome_counts = pd.read_csv(genome_counts_filepath, sep='\t', header=0)
    genome_counts.rename(columns={'Sample': 'MAG_ID'}, inplace=True)

    # Normalize the counts
    if (read_counts_filepath is not False) and (assembly_stats_filepath is not False):
        logger.error('Cannot provide both -r and -R; please be clear on which normalization method you want. Exiting...')
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
    if (CAT_taxonomy_filepath is not False) and (gtdbtk_taxonomy_filepaths is not False):
        logger.error('Must specify only -T or -t, or neither. You specified both flags. Exiting...')
        sys.exit(1)
    elif gtdbtk_taxonomy_filepaths is not False:
        # Load and add taxonomy
        logger.info('Adding GTDB taxonomy info to the table')
        # TODO - remove tuple requirement and make the function receive lists instead
        taxonomy_table = load_and_parse_gtdbtk_taxonomy_table(tuple(gtdbtk_taxonomy_filepaths), resolve=True)
        # TODO - expose 'resolve' setting to user
        genome_counts_normalized = pd.merge(genome_counts_normalized, taxonomy_table,
                                            how = 'left', on = 'MAG_ID',
                                            sort = False, validate = 'one_to_one')

        # Sort by taxonomy
        logger.info('Sorting table by taxonomy')
        genome_counts_normalized = genome_counts_normalized.sort_values(
            by=['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'MAG_ID'])
    elif CAT_taxonomy_filepath is not False:
        # Load and add taxonomy
        logger.info('Adding CAT taxonomy info to the table')
        taxonomy_table = load_CAT_taxonomy_table(CAT_taxonomy_filepath)
        # TODO - add 'resolve' setting once available in the function itself
        genome_counts_normalized = pd.merge(genome_counts_normalized, taxonomy_table,
                                            how='left', on='MAG_ID',
                                            sort=False, validate='one_to_one')

        # Sort by taxonomy
        logger.info('Sorting table by taxonomy')
        genome_counts_normalized = genome_counts_normalized.sort_values(
            by=['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'MAG_ID'])

    # Add on CheckM info
    if checkm_completeness_filepath is not False:
        # Load and add CheckM completeness info
        logger.info('Adding CheckM completeness info to the table')
        checkm_completeness_table = load_checkm_completeness_table(checkm_completeness_filepath)
        genome_counts_normalized = pd.merge(genome_counts_normalized, checkm_completeness_table,
                                            how = 'left', on = 'MAG_ID',
                                            sort = False, validate = 'one_to_one')

    # Write output
    logger.info('Writing output table')
    pd.DataFrame.to_csv(genome_counts_normalized, output_table_filepath, sep = '\t', index = False)

    logger.info(os.path.basename(sys.argv[0]) + ': done.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Creates a TSV-format table of MAG abundances from ATLAS2. '
                    'Copyright Jackson M. Tsuji, Neufeld Research Group, 2019.')

    parser.add_argument('-o', '--output_file', required=True,
                        help='The path to the output TSV-format MAG table')
    parser.add_argument('-a', '--atlas_dir', required=False, default = False,
                        help='The path to the ATLAS directory. If you set this, it will set all five inputs file '
                             'flags (-g, -r, -s, -c, -t) automatically.')
    parser.add_argument('-g', '--genome_counts', required=False, default=False,
                        help='The path to the genome counts file "raw_counts_genomes.tsv" generated by ATLAS.')
    parser.add_argument('-r', '--reads_total', required=False, default=False,
                        help='The path to the total read counts file "read_counts.tsv" generated by ATLAS. '
                        'Specify this flag if you want to normalize to the total number of unassembled reads. '
                        'Must specify either this or -R.')
    parser.add_argument('-R', '--reads_assembled', required=False, default=False,
                        help='The path to the assembly stats file "combined_contig_stats.tsv" generated by ATLAS.'
                             'Specify this flag if you want to normalize to the number of assembled reads. '
                             'Must specify either this or -r.')
    parser.add_argument('-t', '--cat_taxonomy_table', required=False, default=False,
                        help='Paths to the classification output file from the CAT classifier, '
                             'e.g., "taxonomy_names.tsv". Cannot specify both this and -T')
    parser.add_argument('-T', '--gtdb_taxonomy_tables', required=False, default=False, action='append',
                        help='Paths to the two classification output files from the GTDBTk classifier, '
                             'e.g., "gtdbtk.ar122.summary.tsv" and "gtdbtk.bac120.summary.tsv". '
                             'You have to specify each one separately and use this flag twice. '
                             'Cannot specify both this and -t.')
    parser.add_argument('-c', '--checkm_table', required=False, default=False,
                        help='The path to the checkm file "completeness.tsv" generated by ATLAS.')

    command_line_args = parser.parse_args()
    main(command_line_args)