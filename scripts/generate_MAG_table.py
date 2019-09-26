#!/usr/bin/env python
# Generate MAG table of percentages of mapped reads from ATLAS2 output
# Copyright Jackson M. Tsuji, Neufeld Research Group, 2019

# TODO - document functions

# Imports
import sys
import os
import time
import logging
import argparse

import pandas as pd
import re

# GLOBAL variables (default relative paths within the ATLAS dir)
GENOME_COUNTS_FILEPATH = "genomes/counts/raw_counts_genomes.tsv"
READ_COUNTS_FILEPATH = "stats/read_counts.tsv"
CAT_TAXONOMY_FILEPATH = 'genomes/taxonomy/taxonomy_names.tsv'
CHECKM_COMPLETENESS_FILEPATH = 'genomes/checkm/completeness.tsv'

# Set up the logger
logging.basicConfig(level=logging.INFO, format='[ %(asctime)s UTC ]: %(levelname)s: %(message)s')
logging.Formatter.converter = time.gmtime
logger = logging.getLogger(__name__)

# Load the read_counts.tsv table, make custom read count stat, and get key columns
def load_read_counts(read_counts_filepath):
    # Import read totals and filter to QC only with pe and se reads
    read_counts = pd.read_csv(read_counts_filepath, sep = '\t', header = 0)
    read_counts = read_counts[read_counts['Step'] == 'QC']
    read_counts = read_counts[['Sample', 'Reads_pe', 'Reads_se']]

    # Re-calculate totals in a way compatible with read mapping output (i.e., total = R1 + R2 + se)
    read_counts['Reads_total'] = read_counts['Reads_pe'] * 2 + read_counts['Reads_se']
    read_counts = read_counts[['Sample', 'Reads_total']]

    return(read_counts)

# Load the combined_contig_stats.tsv table and get key columns
def load_assembly_stats(assembly_stats_filepath):
    # Import assembly stats and filter to assembled reads column
    assembly_stats = pd.read_csv(assembly_stats_filepath, sep = '\t', header = 0)
    assembly_stats.rename(columns = {assembly_stats.columns[0]: 'Sample'}, inplace = True)
    assembly_stats = assembly_stats[['Sample', 'Assembled_Reads']]

    return(assembly_stats)

# Load and bind rows of GTDBTk taxonomy tables
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

# Parse a single GTDBTk taxonomy entry
def parse_gtdbtk_taxonomy_entry(taxonomy_entry, resolve = True):
    # EXAMPLE: 'd__Bacteria;p__Actinobacteriota;c__Acidimicrobiia;o__Microtrichales;f__;g__;s__'
    
    taxonomy_split = str(taxonomy_entry).split(sep = ';')
    if len(taxonomy_split) != 7:
        logger.error('Taxonomy entry is ' + str(len(taxonomy_split)) + ' long, not 7 as expected. Exiting...')
        loger.error('Entry was: ' + ', '.join(taxonomy_entry))
        sys.exit(1)

    # Remove header pieces
    # TODO - confirm they are in the right order (d, p, c, o, f, g, s)
    taxonomy_split = [ re.sub("[dpcofgs]__", "", level) for level in taxonomy_split ]

    # Fill in empty parts, if they exist
    if '' in taxonomy_split and resolve is True:
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
            logger.error('There seems to be an empty entry in the middle of your taxonomy levels. Cannot resolve. Exiting...')
            logger.error('Entry was: ' + ', '.join(taxonomy_entry))
            sys.exit(1)

        filler_entry = 'Unresolved_' + taxonomy_split[(first_empty_taxon-1)]

        for taxonomy_level_index in range(first_empty_taxon, 7):
            taxonomy_split[taxonomy_level_index] = filler_entry

    return(taxonomy_split)

# Load and parse any number of GTDBTk taxonomy tables (recommended: 2 tables, bacteria and archaea)
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

# Parse a single taxonomy entry
def parse_CAT_taxonomy_cell(CAT_taxonomy_cell):
    # e.g., 'Nitrospira: 0.88' -> 'Nitrospira'
    # e.g., 'not classified' -> 'not classified' (no change)

    if (CAT_taxonomy_cell == 'not classified') or (CAT_taxonomy_cell == 'NA'):
        logger.debug('Skipping empty taxonomy entry')
    elif CAT_taxonomy_cell.find(':') == -1:
        logger.warning('The taxonomy cell "' + CAT_taxonomy_cell + '" seems to be missing the expected colon. ' +
                           'Cannot parse this properly, but will continue on nevertheless, leaving as-is')
    else:
        # Split by the colon
        cell_split = CAT_taxonomy_cell.split(sep=':')
        if len(cell_split) != 2:
            logger.warning('The taxonomy cell "' + CAT_taxonomy_cell + '" seems to not have the expected layout. ' +
                           'Will continue on and leave the cell as-is.')
        else:
            CAT_taxonomy_cell = cell_split[0]

    return(CAT_taxonomy_cell)

# Parse one row of length-7 CAT taxonomy entries
# Input a list of taxonomy entries (string list) and output a list of the corrected entries
def parse_CAT_taxonomy_row(CAT_taxonomy_row, resolve = True):
    # E.g., Bacteria: 1.00    Nitrospirae: 0.95    Nitrospira: 0.91    Nitrospirales: 0.91    Nitrospiraceae: 0.91    not classified    not classified
    # Goal: Bacteria    Nitrospirae    Nitrospira    Nitrospirales    Nitrospiraceae    Unresolved_Nitrospiraceae    Unresolved_Nitrospiraceae

    # Correct each cell
    CAT_taxonomy_row = list(map(parse_CAT_taxonomy_cell, CAT_taxonomy_row))

    # Fill in empty parts, if they exist
    if ('not classified' in CAT_taxonomy_row) or ('NA' in CAT_taxonomy_row) and resolve is True:

        # Get the positions of the empty spots
        empty_taxa = []
        for taxonomy_level in CAT_taxonomy_row:
            if (taxonomy_level == 'not classified') or (taxonomy_level == 'NA'):
                empty_taxa.append(True)
            else:
                empty_taxa.append(False)

        # Get the numeric index of the first empty taxon
        # See https://stackoverflow.com/a/9542768, accessed Sept. 18, 2019
        first_empty_taxon = empty_taxa.index(True)

        # Only look at entries 1-6 for the initial check; sometimes CAT will leave a species level classification in at position 7
        # But CAT does some odd things with taxonomy sometimes, so there might be a difficult to parse entry as well.
        if False in empty_taxa[first_empty_taxon:len(empty_taxa)-1]:
            # TODO - is there a more elegant way to handle these entries?
            logger.warning('There seems to be an empty entry in the middle of one of your taxonomy levels. Cannot resolve...')
            logger.warning('Entry will remain as: "' + ', '.join(CAT_taxonomy_row) + '"')
        else:
            filler_entry = 'Unresolved_' + CAT_taxonomy_row[(first_empty_taxon - 1)]
            # Fill in entries up to 6 with the filler entry
            for taxonomy_level_index in range(first_empty_taxon, 6):
                CAT_taxonomy_row[taxonomy_level_index] = filler_entry
            # Also fill in entry 7 (remember, zero index!) with the filler if there is no special entry there
            if empty_taxa[6] is True:
                CAT_taxonomy_row[6] = filler_entry

    return(CAT_taxonomy_row)

# Load and parse the CAT taxonomy table
def load_and_parse_CAT_taxonomy_table(CAT_taxonomy_table_filepath, resolve = True):
    # Load table
    # Need 'keep_default_na = False' for the support functions to work.
    taxonomy_table = pd.read_csv(CAT_taxonomy_table_filepath, sep = '\t', header = 0, keep_default_na = False)
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

    # Now parse the table
    parsed_rows = []
    for index, taxonomy_table_row in taxonomy_table.iterrows():
        parsed_row = [taxonomy_table_row[0]] + parse_CAT_taxonomy_row(taxonomy_table_row[1:8], resolve = resolve)
        parsed_rows.append(parsed_row)

    taxonomy_table = pd.DataFrame.from_records(parsed_rows, columns =
                                               ['MAG_ID', 'Domain', 'Phylum', 'Class',
                                                'Order', 'Family', 'Genus', 'Species'])

    return(taxonomy_table)

# Load CheckM completeness data and get key columns
def load_checkm_completeness_table(checkm_completeness_table_filepath):

    # Load CheckM data
    checkm_completeness_table = pd.read_csv(checkm_completeness_table_filepath, sep = '\t', header = 0)
    checkm_completeness_table.rename(columns = {'Bin Id': 'MAG_ID',
                                               'Strain heterogeneity': 'Strain_heterogeneity'},
                                     inplace = True)
    checkm_completeness_table = checkm_completeness_table[['MAG_ID', 'Completeness', 'Contamination', 'Strain_heterogeneity']]

    return(checkm_completeness_table)

# Normalize the genome count table to a particular table of total reads
def normalize_genome_count_table(genome_count_table, totals_table):
    normalized_totals = []
    
    for row_tuple in totals_table.itertuples(index = False):
        sample_ID = row_tuple[0]
        total_reads = row_tuple[1]

        normalized_totals.append(genome_count_table[sample_ID] / total_reads * 100)
        #count_table_normalized[sample_ID] = count_table[sample_ID] / total_reads * 100
        
    # Bind columns and add MAG IDs
    genome_count_table_normalized = pd.concat(normalized_totals, axis = 1)
    genome_count_table_normalized.insert(loc = 0, column = 'MAG_ID',
                                         value = genome_count_table['MAG_ID'].tolist())
    
    # See the summed totals of percent mapped reads
    column_sums = genome_count_table_normalized.drop(columns = 'MAG_ID').apply(sum, axis = 0)
    column_sums = pd.DataFrame(column_sums, columns = ['normalized sum'])
    logger.info("Total normalized values:")
    for rowname, row_data in column_sums.iterrows():
        logger.info(rowname + ': ' + str(round(row_data[0], 1)))
    
    return(genome_count_table_normalized)

def main(args):
    # Set user variables
    atlas_dir = args.atlas_dir
    output_table_filepath = args.output_file
    genome_counts_filepath = args.genome_counts
    read_counts_filepath = args.reads_total
    assembly_stats_filepath = args.reads_assembled
    CAT_taxonomy_filepath = args.cat_taxonomy_table
    gtdbtk_taxonomy_filepaths = args.gtdb_taxonomy_tables
    checkm_completeness_filepath = args.checkm_table

    # Auto set non-specified key flags if atlas_dir has been set
    # But keep any manual user overrides
    if atlas_dir is not False:
        if genome_counts_filepath is False:
            genome_counts_filepath = os.path.join(atlas_dir, GENOME_COUNTS_FILEPATH)
        if (read_counts_filepath is False) and (assembly_stats_filepath is False):
            # Set read_counts_filepath by default
            read_counts_filepath = os.path.join(atlas_dir, READ_COUNTS_FILEPATH)
        if (gtdbtk_taxonomy_filepaths is False) and (CAT_taxonomy_filepath is False):
            CAT_taxonomy_filepath = os.path.join(atlas_dir, CAT_TAXONOMY_FILEPATH)
        if checkm_completeness_filepath is False:
            checkm_completeness_filepath = os.path.join(atlas_dir, CHECKM_COMPLETENESS_FILEPATH)

    # Startup messages
    logger.info('Running ' + os.path.basename(sys.argv[0]))
    logger.info('##### Settings #####')
    logger.info('ATLAS directory:   ' + str(atlas_dir))
    logger.info('Output filepath: ' + str(output_table_filepath))
    logger.info('Genome counts filepath: ' + str(genome_counts_filepath))
    logger.info('Read counts filepath: ' + str(read_counts_filepath))
    logger.info('Assembly stats filepath: ' + str(assembly_stats_filepath))
    if gtdbtk_taxonomy_filepaths is not False:
        logger.info('GTDB taxonomy filepaths: ' + ', '.join(list(gtdbtk_taxonomy_filepaths)))
    else:
        logger.info('GTDB taxonomy filepaths: ' + str(gtdbtk_taxonomy_filepaths))
    logger.info('CAT taxonomy filepath: ' + str(CAT_taxonomy_filepath))
    logger.info('CheckM completeness filepath: ' + str(checkm_completeness_filepath))
    logger.info('####################')

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
        taxonomy_table = load_and_parse_CAT_taxonomy_table(CAT_taxonomy_filepath, resolve=True)
        # TODO - expose 'resolve' setting to the user
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
    parser.add_argument('-t', '--cat_taxonomy_table', required=False, default=False,
                        help='Path to the classification output file from the CAT classifier, '
                             'e.g., "taxonomy_names.tsv". Cannot specify both this and -T')
    parser.add_argument('-T', '--gtdb_taxonomy_tables', required=False, default=False, nargs='+',
                        help='Paths to the two classification output files from the GTDBTk classifier, '
                             'separated by spaces e.g., "gtdbtk.ar122.summary.tsv gtdbtk.bac120.summary.tsv". '
                             'Cannot specify both this and -t.')
    parser.add_argument('-c', '--checkm_table', required=False, default=False,
                        help='The path to the checkm file "completeness.tsv" generated by ATLAS.')

    command_line_args = parser.parse_args()
    main(command_line_args)