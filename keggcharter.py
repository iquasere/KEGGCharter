#!/usr/bin/env python

from argparse import ArgumentParser, ArgumentTypeError
import numpy as np
import os
import pandas as pd
from pathlib import Path
from subprocess import run
import sys
from io import StringIO
from time import time, gmtime, strftime
from Bio.KEGG.REST import kegg_link, kegg_list, kegg_get
from Bio.KEGG.KGML import KGML_parser
from matplotlib import pyplot as plt
from tqdm import tqdm
from glob import glob
import json
import re

from keggpathway_map import KEGGPathwayMap, expand_by_list_column

__version__ = "0.7.0"


def get_arguments():
    parser = ArgumentParser(
        description="""KEGGCharter - A tool for representing genomic potential and transcriptomic expression into 
        KEGG pathways""", epilog="Input file must be specified.")
    parser.add_argument("-o", "--output", help="Output directory", default='KEGGCharter_results')
    parser.add_argument(
        "-rd", "--resources-directory", default=sys.path[0], help="Directory for storing KGML and CSV files.")
    parser.add_argument(
        "-mm", "--metabolic-maps", help="IDs of metabolic maps to output",
        default=','.join(keggcharter_prokaryotic_maps()))
    parser.add_argument("-qcol", "--quantification-columns", help="Names of columns with quantification")
    parser.add_argument(
        "-tls", "--taxa-list", help="List of taxa to represent in genomic potential charts (comma separated)")  # TODO - must be tested
    parser.add_argument(
        "-not", "--number-of-taxa", type=int, default=10,
        help="Number of taxa to represent in genomic potential charts (comma separated)")
    parser.add_argument("-keggc", "--kegg-column", help="Column with KEGG IDs.")
    parser.add_argument("-koc", "--ko-column", help="Column with KOs.")
    parser.add_argument("-ecc", "--ec-column", help="Column with EC numbers.")
    # TODO - test this argument without UniProt shenanigans
    parser.add_argument(
        "-tc", "--taxa-column", default='Taxonomic lineage (GENUS)',
        help="Column with the taxa designations to represent with KEGGCharter."
             " NOTE - for valid taxonomies, check: https://www.genome.jp/kegg/catalog/org_list.html")
    parser.add_argument(
        "-iq", "--input-quantification", action="store_true",
        help="If input table has no quantification, will create a mock quantification column")
    parser.add_argument(
        "-it", "--input-taxonomy", default=None,
        help="If no taxonomy column exists and there is only one taxon in question.")
    parser.add_argument(
        "--step", default=40, type=int, help="Number of IDs to submit per request through the KEGG API [40]")
    parser.add_argument(
        "--map-all", default=False, action="store_true",
        help="Ignore KEGG taxonomic information. All functions for all KOs will be represented,"
             " even if they aren't attributed by KEGG to the specific species.")
    parser.add_argument(
        "--include-missing-genomes", default=False, action="store_true",
        help="Map the functions for KOs identified for organisms not present in KEGG Genomes.")
    parser.add_argument(
        "--resume", action="store_true", default=False,
        help="If data inputed has already been analyzed by KEGGCharter.")
    parser.add_argument('-v', '--version', action='version', version='KEGGCharter ' + __version__)

    required_named = parser.add_argument_group('required named arguments')
    required_named.add_argument("-f", "--file", required=True, help="TSV or EXCEL table with information to chart")

    special_functions = parser.add_argument_group('Special functions')
    special_functions.add_argument(
        "--show-available-maps", action="store_true", default=False,
        help="Outputs KEGG maps IDs and descriptions to the console (so you may pick the ones you want!)")

    args = parser.parse_args()

    args.output = args.output.rstrip('/')
    for directory in [args.output] + [f'{args.resources_directory}/{folder}' for folder in ['', 'kc_kgmls', 'kc_csvs']]:
        if not os.path.isdir(directory):
            Path(directory).mkdir(parents=True, exist_ok=True)
            print(f'Created {directory}')
    if not hasattr(args, 'quantification_columns'):
        input_quantification = str2bool(
            'No quantification columns specified! Do you want to use mock quantification? '
            'See https://github.com/iquasere/KEGGCharter#mock-imputation-of-quantification-and-taxonomy for more '
            'details. [y/N]')
        if input_quantification:
            args.input_quantification = True
        else:
            exit('No quantification columns specified!')
    return args


def str2bool(v):
    if v.lower() == 'auto':
        return 'auto'
    elif v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise ArgumentTypeError('Boolean value expected.')


def timed_message(message):
    print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ': ' + message)


def run_command(bash_command, print_command=True, stdout=None, stderr=None):
    if print_command:
        print(bash_command)
    run(bash_command.split(), stdout=stdout, stderr=stderr)


def download_organism(directory):
    if not os.path.isfile(f"{directory}/organism"):
        run_command(f'wget http://rest.kegg.jp/list/organism -P {directory}')


def read_input_file(file):
    timed_message('Reading input data.')
    if file.endswith('.xlsx'):
        return pd.read_excel(file)
    return pd.read_csv(file, sep='\t', low_memory=False)


def further_information(data, output, kegg_column=None, ko_column=None, ec_column=None, step=150):
    """
    Adds KOs and EC numbers to the input data
    """
    data = get_cross_references(data, kegg_column=kegg_column, ko_column=ko_column, ec_column=ec_column, step=step)
    main_column = kegg_column if kegg_column is not None else ko_column if ko_column is not None else ec_column
    data = condense_data(data, main_column)
    data.to_csv(output, sep='\t', index=False)
    timed_message(f'Results saved to {output}')
    return data, main_column


# Conversion functions
def id2id(input_ids, in_col, out_col, in_type, out_type, step=150):
    """
    Converts KEGG IDs, KOs or EC numbers to KOs and EC numbers through the KEGG API
    :param input_ids: (list) - IDs to convert
    :param in_col: (str) - name of column with IDs to convert (used for merging DFs)
    :param out_col: (str) - name of column to return
    :param in_type: (str) - type of input IDs
    :param out_type: (str) - database to convert IDs to
    :param step: (int) - will convert "step" KEGG IDs at a time
    :return: (list) - (list,list) - KEGG ID genes converted and ko IDs
    """
    result = pd.DataFrame(columns=[in_col, out_col])
    if len(input_ids) == 0:
        return result
    type_to_desc = {'kegg': 'KEGG IDs', 'ko': 'KOs', 'ec': 'EC numbers', 'enzyme': 'EC numbers'}
    desc = f'Converting {len(input_ids)} {type_to_desc[in_type]} to {type_to_desc[out_type]} through the KEGG API'
    for i in tqdm(range(0, len(input_ids), step), desc=desc, ascii=' >='):
        j = min(i + step, len(input_ids))
        try:
            result = pd.concat([result, pd.read_csv(
                kegg_link(out_type, input_ids[i:j]), sep='\t', names=[in_col, out_col])])
        except:
            try:        # try a second time. The API is robust, and it's unusual to fail a third time
                result = pd.concat([result, pd.read_csv(
                    kegg_link(out_type, input_ids[i:j]), sep='\t', names=[in_col, out_col])])
            except Exception as e:
                print(f'IDs conversion broke at index {i}; Error: {e}')
    if in_type == 'ko':
        result[in_col] = result[in_col].apply(lambda x: x.split('ko:')[-1])
    elif in_type == 'ec':
        result[in_col] = result[in_col].apply(lambda x: x.split('ec:')[-1])
    if out_type == 'ko':
        result[out_col] = result[out_col].apply(lambda x: x.split('ko:')[-1])
    elif out_type == 'enzyme':
        result[out_col] = result[out_col].apply(lambda x: x.split('ec:')[-1])
    result = result.groupby(in_col)[out_col].agg(lambda x: ','.join(map(str, x))).reset_index()
    return result


def ids_xref(data, in_col, out_col, in_type='kegg', step=150):
    data[f'{in_col}_split'] = data[in_col].apply(lambda x: x.split(',') if type(x) != float else x)    # split by comma
    data = expand_by_list_column(data, column=f'{in_col}_split')
    ids = ','.join(data[f'{in_col}_split'].dropna().unique()).split(',')       # KEGGCharter only accepts "," as separator because KEGG uses the same separator
    if in_type == 'kegg':
        new_ids = id2id(ids, f'{in_col}_split', out_col, in_type='kegg', out_type='ko', step=step)
    elif in_type == 'ko':
        new_ids = id2id(ids, f'{in_col}_split', out_col, in_type='ko', out_type='enzyme', step=step)
    elif in_type == 'ec':
        new_ids = id2id(ids, f'{in_col}_split', out_col, in_type='ec', out_type='ko', step=step)
    else:
        raise ValueError('ids_type must be one of: kegg, ko, ec')
    data = pd.merge(data, new_ids, on=f'{in_col}_split', how='outer')
    data = data.drop(columns=[f'{in_col}_split'])
    other_cols = [col for col in data.columns if col not in [in_col, out_col]]
    result = pd.concat([
        data[other_cols].drop_duplicates(),
        data.groupby(in_col)[out_col].agg(lambda x: ','.join(map(str, set([val for val in x if type(val) != float])))
                                          ).reset_index()], axis=1)[data.columns].replace('', np.nan)
    return result


def get_cross_references(data, kegg_column=None, ko_column=None, ec_column=None, step=150):
    # KEGG ID to KO -> if KO column is not set, KEGGCharter will get them through the KEGG API
    ko_cols = []
    ec_cols = []
    if kegg_column:
        data = ids_xref(data, in_col=kegg_column, out_col='KO (kegg-column)', in_type='kegg', step=step)
        data = ids_xref(data, in_col='KO (kegg-column)', out_col='EC (kegg-column)', in_type='ko', step=step)
        ko_cols.append('KO (kegg-column)'); ec_cols.append('EC (kegg-column)')
    if ko_column:
        data = ids_xref(data, in_col=ko_column, out_col='EC (ko-column)', in_type='ko', step=step)
        data = ids_xref(data, in_col='EC (ko-column)', out_col='KO (ko-column)', in_type='ec', step=step)
        ko_cols.append(ko_column); ko_cols.append('KO (ko-column)'); ec_cols.append('EC (ko-column)')
    if ec_column:
        data = ids_xref(data, in_col=ec_column, out_col='KO (ec-column)', in_type='ec', step=step)
        data = ids_xref(data, in_col='KO (ec-column)', out_col='EC (ec-column)', in_type='ko', step=step)
        ko_cols.append('KO (ec-column)'); ec_cols.append(ec_column); ec_cols.append('EC (ec-column)')
    # join all unique KOs in a column
    data['KO (KEGGCharter)'] = data[ko_cols].apply(
        lambda x: ','.join(set([elem for elem in x if elem is not np.nan])), axis=1)
    data['KO (KEGGCharter)'] = data['KO (KEGGCharter)'].apply(lambda x: ','.join(sorted(set(x.split(',')))))
    # join all unique ECs in a column
    data['EC number (KEGGCharter)'] = data[ec_cols].apply(
        lambda x: ','.join(set([elem for elem in x if elem is not np.nan])), axis=1)
    data['EC number (KEGGCharter)'] = data['EC number (KEGGCharter)'].apply(
        lambda x: ','.join(sorted(set(x.split(',')))))
    if not (kegg_column or ko_column or ec_column):
        exit('Need to specify a column with either KEGG IDs, KOs or EC numbers!')
    return data


def condense_data(data, main_column):
    onlykos = data[data['KO (KEGGCharter)'].notnull() & (data['EC number (KEGGCharter)'].isnull())][
        [main_column, 'KO (KEGGCharter)']]
    onlykos = onlykos.groupby(main_column).agg({'KO (KEGGCharter)': lambda x: ','.join(set(x))}).reset_index()
    onlykos['EC number (KEGGCharter)'] = [np.nan] * len(onlykos)
    wecs = data[data['EC number (KEGGCharter)'].notnull()][[main_column, 'KO (KEGGCharter)', 'EC number (KEGGCharter)']]
    wecs = wecs.groupby(main_column).agg(
        {'KO (KEGGCharter)': lambda x: ','.join(set([elem for elem in x if elem is not np.nan])),
         'EC number (KEGGCharter)': lambda x: ','.join(set(x))}).reset_index()
    del data['KO (KEGGCharter)']
    del data['EC number (KEGGCharter)']
    return pd.merge(data, pd.concat([onlykos, wecs]), on=main_column, how='left').drop_duplicates()


def prepare_data_for_charting(data, mt_cols=None, ko_column='KO (KEGGCharter)', ko_from_uniprot=False):
    """
    This function expands the dataframe by the KO column, so that each row has only one KO.
    """
    # TODO - this function is outdated, refresh it
    nokos = data[data[ko_column].isnull()]
    wkos = data[data[ko_column].notnull()].reset_index(drop=True)
    if ko_from_uniprot:
        # If coming from UniProt ID mapping, KOs will be in the form "KXXXXX;"
        wkos[ko_column] = wkos[ko_column].apply(lambda x: x.rstrip(';'))
    # Expand KOs column if some elements are in the form KO1,KO2,...
    wkos[ko_column] = wkos[ko_column].apply(lambda x: x.split(','))
    # Divide the quantification by the number of KOs in the column
    if mt_cols is not None:
        for col in mt_cols:
            wkos[col] = wkos[col] / wkos[ko_column].apply(lambda x: len(x))
    wkos = expand_by_list_column(wkos, column=ko_column)
    data = pd.concat([wkos, nokos])
    timed_message('Data expanded by KO column.')
    return data


# Get metabolic maps from KEGG Pathway
def keggcharter_prokaryotic_maps(file=f'{sys.path[0]}/KEGGCharter_prokaryotic_maps.txt'):
    return open(file).read().split('\n')


def kegg_metabolic_maps():
    """
    Creates a dic with all specific kegg pathway maps and their description
    :return: pandas.DataFrame with Map ID as index and maps names as
    sole column
    """
    maps = pd.read_csv(StringIO(kegg_list("pathway").read()), sep='\t', names=['Map ID', 'Description'])
    maps['Map ID'] = [ide.split('map')[-1] for ide in maps['Map ID']]
    return maps


def write_kgml(mmap, output, organism='ko', max_tries=3):
    """
    This function is confusing, in that it both writes the KGML, and parses it.
    Still, it works, and for now that's enough.
    """
    tries = 0
    while tries < max_tries:
        try:
            data = kegg_get(f"{organism}{mmap}", "kgml").read()
            with open(output, 'w') as f:
                if len(data) > 1:
                    f.write(data)
                    return KGML_parser.read(data)
                return None
        except Exception as e:
            print(f'Error: {e}.{" Trying again..." if tries < max_tries else ""}')
            tries += 1
    return None


def glob_re(pattern, strings):
    return filter(re.compile(pattern).match, strings)


def write_kgmls(mmaps, out_dir, max_tries=3, org='ko'):
    maps_done = [
        filename.split(org)[-1].rstrip('.xml') for filename in glob_re(fr'{org}\d+\.xml', os.listdir(out_dir))]
    mmap_to_orthologs = {}
    for name in maps_done:          # maps already done will have their orthologs and genes put in
        parsed = KGML_parser.read(open(f'{out_dir}/{org}{name}.xml'))
        mmap_to_orthologs[name] = [orth.id for orth in parsed.orthologs] + [gene.id for gene in parsed.genes]
    mmaps = [mmap for mmap in mmaps if mmap not in maps_done]
    i = 1
    if len(mmaps) == 0:
        return mmap_to_orthologs
    for mmap in tqdm(mmaps, desc=f'Getting [{len(mmaps)}] KGMLs for taxon [{org}]', ascii=' >='):
        tries = 0
        done = False
        while tries < max_tries and not done:
            orthologs = [orth.id for orth in write_kgml(mmap, f'{out_dir}/{org}{mmap}.xml', organism=org).orthologs]
            genes = [gene.id for gene in write_kgml(mmap, f'{out_dir}/{org}{mmap}.xml', organism=org).genes]
            mmap_to_orthologs[mmap] = orthologs + genes
            done = True
        i += 1
    return mmap_to_orthologs


def set_text_boxes_kgml(kgml_filename, out_filename, desc=''):
    handler = KGML_parser.read(open(kgml_filename))
    # Set text in boxes to EC numbers
    with open(out_filename, 'w') as f:
        for ortholog_rec in tqdm(handler.orthologs, desc=desc, ascii=' >='):
            lines = []
            kos = ortholog_rec.name.split()
            lines += kegg_link("enzyme", kos).read().split('\n')
            ecs = [line.split('\t')[1] for line in lines if len(line) > 0]
            if len(ecs) > 0:
                f.write(f'{",".join(ecs)}\n')
            else:
                f.write(f'{",".join(kos)}\n')


def taxon2prefix(taxon_name, organism_df):
    """
    :param taxon_name: str - e.g. Pongo abelii (Sumatran orangutan)
    :param organism_df: pandas.DataFrame - organism file, index = taxa names, column names = KEGG prefix, ...
    :returns str - KEGG prefix of taxon name
    """
    if taxon_name == 'ko':
        return 'ko'
    if taxon_name in organism_df.index:                     # taxon is found as it is
        if len(organism_df.loc[taxon_name]) > 1:
            return organism_df.loc[taxon_name, 'prefix'][0]
        return organism_df.loc[taxon_name, 'prefix']
    if taxon_name.split(' (')[0] in organism_df.index:      # Homo sapiens (human) -> Homo sapiens
        if len(organism_df.loc[taxon_name.split(' (')[0]]) > 1:
            return organism_df.loc[taxon_name.split(' (')[0], 'prefix'][0]
        return organism_df.loc[taxon_name.split(' (')[0], 'prefix']
    possible_prefixes = organism_df[organism_df.index.str.contains(taxon_name)].prefix.tolist()
    if len(possible_prefixes) > 0:
        return possible_prefixes[0]                         # select the first strain
    return None         # not found in taxon to KEGG prefix conversion


def get_taxon_maps(kegg_prefix, max_tries=3):
    if kegg_prefix is None:
        return []
    tries = 0
    while tries < max_tries:
        try:
            df = pd.read_csv(StringIO(kegg_list("pathway", kegg_prefix).read()), sep='\t', header=None)
            return df[[0]].apply(lambda x: x.split(kegg_prefix)[1]).tolist()
        except:
            tries += 1
    return []


def parse_organism(file):
    return pd.read_csv(file, sep='\t', usecols=[1, 2], header=None, index_col=1, names=['prefix', 'name'])


def download_resources(
        data, resources_directory, taxa_column, metabolic_maps, map_all=False, map_non_kegg_genomes=True):
    """
    Download all resources for a given dataframe
    :param data: pandas.DataFrame - dataframe with taxa names in taxa_column
    :param resources_directory: str - directory where to save the resources
    :param taxa_column: str - column name in dataframe with taxa names
    :param metabolic_maps: list - metabolic maps to download
    :param map_all: bool - if True, attribute all maps and all functions to all taxa, only limit by the identifications
    :param map_non_kegg_genomes: bool - if True, map non-KEGG genomes to KEGG orthologs
    :return: taxon_to_mmap_to_orthologs - dic with taxon name as key and dic with metabolic maps as values
    """
    timed_message('Downloading resources')
    download_organism(resources_directory)
    taxa = ['ko'] + data[taxa_column].unique().tolist()
    if np.nan in taxa:
        taxa.remove(np.nan)
    taxa_df = parse_organism(f'{resources_directory}/organism')
    taxon_to_mmap_to_orthologs = {}  # {'Keratinibaculum paraultunense' : {'00190': ['1', '2']}}
    if map_all:     # attribute all maps and all functions to all taxa, only limit by the data
        mmap_to_orthologs = write_kgmls(metabolic_maps, f'{resources_directory}/kc_kgmls', org='ko')
        taxon_to_mmap_to_orthologs = {taxon: mmap_to_orthologs for taxon in taxa}
    else:
        kegg_prefixes = [(taxon, taxon2prefix(taxon, taxa_df)) for taxon in tqdm(
            taxa, desc='Obtaining KEGG prefixes from inputted taxa', ascii=' >=')]
        kegg_prefixes = [pair for pair in kegg_prefixes if pair[1] is not None]
        for taxon, kegg_prefix in tqdm(
                kegg_prefixes, desc=f'Getting information for {len(kegg_prefixes) - 1} taxa', ascii=' >='):
            if kegg_prefix is not None:
                taxon_mmaps = get_taxon_maps(kegg_prefix)
                taxon_mmaps = [mmap for mmap in taxon_mmaps if mmap in metabolic_maps]  # select only inputted maps
                taxon_to_mmap_to_orthologs[taxon] = write_kgmls(
                    taxon_mmaps, f'{resources_directory}/kc_kgmls', org=kegg_prefix)
            else:
                if map_non_kegg_genomes:
                    taxon_to_mmap_to_orthologs[taxon] = write_kgmls(
                        metabolic_maps, f'{resources_directory}/kc_kgmls', org='ko')
                else:
                    taxon_to_mmap_to_orthologs[taxon] = {}
    with open(f'{resources_directory}/taxon_to_mmap_to_orthologs.json', 'w') as f:
        json.dump(taxon_to_mmap_to_orthologs, f)
    timed_message('Generated taxon_to_mmap_to_orthologs.json')
    return taxon_to_mmap_to_orthologs


def get_mmaps2taxa(taxon_to_mmap_to_orthologs):
    """
    :param taxon_to_mmap_to_orthologs: dict - {'Keratinibaculum paraultunense' : {'00190': ['1', '2']}}
    :returns dict - {'00190': ['Keratinibaculum paraultunense', 'Keratinibaculum paraultunense']}
    """
    mmaps2taxa = {}
    for org, mmaps2orthologs in taxon_to_mmap_to_orthologs.items():
        for mmap in mmaps2orthologs.keys():
            if mmap in mmaps2taxa.keys():
                mmaps2taxa[mmap].append(org)
            else:
                mmaps2taxa[mmap] = [org]
    return mmaps2taxa


def chart_map(
        kgml_filename, ec_list, data, taxon_to_mmap_to_orthologs, mmaps2taxa, output=None, ko_column=None,
        taxa_column=None, quantification_columns=None, number_of_taxa=10,
        grey_taxa='Other taxa'):
    mmap = KGML_parser.read(open(kgml_filename))
    kegg_pathway_map = KEGGPathwayMap(pathway=mmap, ec_list=ec_list)
    kegg_pathway_map.genomic_potential_taxa(
        data, quantification_columns, ko_column, taxon_to_mmap_to_orthologs, mmaps2taxa=mmaps2taxa,
        taxa_column=taxa_column, output_basename=f'{output}/potential', number_of_taxa=number_of_taxa,
        grey_taxa=grey_taxa)
    mmap = KGML_parser.read(open(kgml_filename))        # need to re-read the file because it's modified by the function
    kegg_pathway_map = KEGGPathwayMap(pathway=mmap, ec_list=ec_list)
    kegg_pathway_map.differential_expression_sample(
        data, quantification_columns, ko_column, mmaps2taxa, taxa_column=taxa_column,
        output_basename=f'{output}/differential', log=False)
    plt.close()


def get_pathway_and_ec_list(rd, mmap):
    download = True
    pathway = None
    if os.path.isfile(f'{rd}/kc_kgmls/ko{mmap}.xml') and os.path.isfile(f'{rd}/kc_csvs/ko{mmap}.csv'):
        pathway = KGML_parser.read(open(f'{rd}/kc_kgmls/ko{mmap}.xml'))
        with open(f'{rd}/kc_csvs/ko{mmap}.csv') as f:
            ec_list = f.read().split('\n')
        if len(pathway.orthologs) == len(ec_list) - 1:  # -1 because of newline at the end
            download = False
        else:
            print(f'Lengths of orthologs in KGML and labels for corresponding boxes do not match for map [ko{mmap}]!')
    else:
        print(f'Some resources were not found for map [ko{mmap}]! Going to download them')
    if download:
        try:
            write_kgml(mmap, f'{rd}/kc_kgmls/ko{mmap}.xml')
            print(f'Got KGML for map [ko{mmap}]')
            set_text_boxes_kgml(
                f'{rd}/kc_kgmls/ko{mmap}.xml', f'{rd}/kc_csvs/ko{mmap}.csv',
                desc=f"Getting boxes' labels for map [ko{mmap}]")
            pathway = KGML_parser.read(open(f'{rd}/kc_kgmls/ko{mmap}.xml'))
            with open(f'{rd}/kc_csvs/ko{mmap}.csv') as f:
                ec_list = f.read().split('\n')
        except Exception as e:
            print(f'Could not download resources for [ko{mmap}]! Error: {e}')
            return None, None
    return pathway, ec_list


def read_input():
    args = get_arguments()
    if args.show_available_maps:
        sys.exit(kegg_metabolic_maps().to_string())
    data = read_input_file(args.file)
    if args.input_quantification:
        data['Quantification (KEGGCharter)'] = [1] * len(data)
        args.quantification_columns = 'Quantification (KEGGCharter)'
    if args.input_taxonomy:
        data['Taxon (KEGGCharter)'] = [args.input_taxonomy] * len(data)
        args.taxa_column = 'Taxon (KEGGCharter)'
        args.taxa_list = args.input_taxonomy
    args.metabolic_maps = args.metabolic_maps.split(',')
    args.quantification_columns = args.quantification_columns.split(',')
    timed_message('Arguments valid.')
    return args, data


def main():
    args, data = read_input()
    if not args.resume:
        data, main_column = further_information(
            data,
            f'{args.output}/KEGGCharter_results.tsv',
            kegg_column=args.kegg_column,
            ko_column=args.ko_column,
            ec_column=args.ec_column,
            step=args.step)
    ko_column = args.ko_column if args.ko_column else 'KO (KEGGCharter)'

    if args.resume:
        data = pd.read_csv(f'{args.output}/data_for_charting.tsv', sep='\t', low_memory=False)
        if not args.input_taxonomy:
            with open(f'{args.output}/taxon_to_mmap_to_orthologs.json') as h:
                taxon_to_mmap_to_orthologs = json.load(h)
        else:
            taxon_to_mmap_to_orthologs = None
    else:
        data = prepare_data_for_charting(data, ko_column=ko_column, mt_cols=args.quantification_columns)
        data.to_csv(f'{args.output}/data_for_charting.tsv', sep='\t', index=False)
        if not args.input_taxonomy:
            taxon_to_mmap_to_orthologs = download_resources(
                data, args.resources_directory, args.taxa_column, args.metabolic_maps, map_all=args.map_all,
                map_non_kegg_genomes=args.include_missing_genomes)
            h = open(f"{args.output}/taxon_to_mmap_to_orthologs.json", "w")
            json.dump(taxon_to_mmap_to_orthologs, h)
        else:
            taxon_to_mmap_to_orthologs = None

    mmaps2taxa = get_mmaps2taxa(taxon_to_mmap_to_orthologs) if not args.input_taxonomy else None        # '00190': ['Keratinibaculum paraultunense']
    timed_message(f'Creating KEGG Pathway representations for {len(args.metabolic_maps)} metabolic pathways.')
    for i in range(len(args.metabolic_maps)):
        pathway, ec_list = get_pathway_and_ec_list(args.resources_directory, args.metabolic_maps[i])
        if pathway is not None and ec_list is not None:
            timed_message(f'[{i + 1}/{len(args.metabolic_maps)}] {pathway.title}')
            chart_map(
                f'{args.resources_directory}/kc_kgmls/ko{args.metabolic_maps[i]}.xml',
                ec_list,
                data,
                taxon_to_mmap_to_orthologs,
                mmaps2taxa,
                output=args.output,
                ko_column=ko_column,
                taxa_column=args.taxa_column,
                quantification_columns=args.quantification_columns,
                number_of_taxa=args.number_of_taxa,
                grey_taxa=('Other taxa' if args.input_taxonomy is None else args.input_taxonomy))
        else:
            print(f'Analysis of map {args.metabolic_maps[i]} failed! Map might have been deleted, '
                  f'for more info raise an issue at https://github.com/iquasere/KEGGCharter/issues')
            i += 1

    # TODO - implement multiprocessing for map generation?
    '''
    with multiprocessing.Pool() as p:
        p.starmap(chart_map, [(handler, data, args.output, ko_column, ec_column, args.taxa_column, dic_colors,
                              args.genomic_columns, args.transcriptomic_columns,
                              f'{args.output}/failed_maps.txt') for handler in pathways])
    '''


if __name__ == '__main__':
    start_time = time()
    main()
    print(f'KEGGCharter analysis finished in {strftime("%Hh%Mm%Ss", gmtime(time() - start_time))}')
