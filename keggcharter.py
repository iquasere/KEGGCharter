#!/usr/bin/env python

from argparse import ArgumentParser
import numpy as np
import os
import pandas as pd
from pathlib import Path
from subprocess import run
import sys
from io import StringIO
from time import time, gmtime, strftime, sleep
from Bio.KEGG.REST import kegg_link, kegg_list, kegg_get
from Bio.KEGG.KGML import KGML_parser
from matplotlib import pyplot as plt
from tqdm import tqdm
from glob import glob
import json

from keggpathway_map import KEGGPathwayMap, expand_by_list_column

__version__ = "0.3.2"


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
    parser.add_argument("-gcol", "--genomic-columns", help="Names of columns with genomic identification")
    parser.add_argument(
        "-tcol", "--transcriptomic-columns", help="Names of columns with transcriptomics quantification")
    parser.add_argument(
        "-tls", "--taxa-list", help="List of taxa to represent in genomic potential charts (comma separated)")  # TODO - must be tested
    parser.add_argument(
        "-not", "--number-of-taxa", type=int, default=10,
        help="Number of taxa to represent in genomic potential charts (comma separated)")
    parser.add_argument("-keggc", "--kegg-column", help="Column with KEGG IDs.")
    parser.add_argument("-koc", "--ko-column", help="Column with KOs.")
    parser.add_argument("-ecc", "--ec-column", help="Column with EC numbers.")
    parser.add_argument(
        "-iq", "--input-quantification", action="store_true",
        help="If input table has no quantification, will create a mock quantification column")
    parser.add_argument(
        "-it", "--input-taxonomy", help="If no taxonomy column exists and there is only one taxon in question.")
    # TODO - test this argument without UniProt shenanigans
    parser.add_argument(
        "-tc", "--taxa-column", default='Taxonomic lineage (GENUS)',
        help="Column with the taxa designations to represent with KEGGCharter")
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

    for directory in [args.output]:
        if not os.path.isdir(directory):
            Path(directory).mkdir(parents=True, exist_ok=True)
            print('Created ' + directory)

    return args


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
    if file.endswith('.xlsx'):
        return pd.read_excel(file)
    return pd.read_csv(file, sep='\t')


def further_information(data, kegg_column=None, ko_column=None, ec_column=None):
    data = get_cross_references(data, kegg_column=kegg_column, ko_column=ko_column, ec_column=ec_column)
    main_column = kegg_column if kegg_column is not None else ko_column if ko_column is not None else ec_column
    data = condense_data(data, main_column)
    return data, main_column


# Conversion functions
def id2id(input_ids, column, output_column, output_ids_type, step=150, desc=''):
    """
    Converts KEGG_ID genes to Ortholog KO ID from KEGG
    :param input_ids: (list) - IDs to convert
    :param column: (str) - name of column with IDs to convert (used for merging DFs)
    :param output_column: (str) - name of column to return
    :param output_ids_type: (str) - database to convert IDs to
    :param step: (int) - will convert "step" KEGG IDs at a time
    :param desc: (str) - string to output to tqdm progressbar
    :return: (list) - (list,list) - KEGG ID genes converted and ko IDs
    """
    result = pd.DataFrame(columns=[column, output_column])
    for i in tqdm(range(0, len(input_ids), step), desc=desc):
        j = min(i + step, len(input_ids))
        try:
            result = pd.concat([result, pd.read_csv(
                kegg_link(output_ids_type, input_ids[i:j]), sep='\t', names=[column, output_column])])
        except:
            print(f'IDs conversion broke at index: {i}')
    if output_ids_type == 'ko':
        result[output_column] = result[output_column].apply(lambda x: x.strip('ko:'))
        result[column] = result[column].apply(lambda x: x.strip('ec:'))
    elif output_ids_type == 'enzyme':
        result[column] = result[column].apply(lambda x: x.strip('ko:'))
        result[output_column] = result[output_column].apply(lambda x: x.strip('ec:'))
    return result


def ids_interconversion(data, column, ids_type='kegg'):
    ids = list(set(data[data[column].notnull()][column]))
    base_desc = 'Converting %i %s to %s through the KEGG API'
    if ids_type == 'kegg':
        # sometimes Kegg IDs come as mfc:BRM9_0145;mfi:DSM1535_1468; (e.g. from UniProt IDs mapping).
        # Should be no problem since both IDs likely will always represent same functions
        trimmed_ids = [ide.split(';')[0] for ide in ids]
        relational = pd.DataFrame([ids, trimmed_ids]).transpose()
        new_ids = id2id(trimmed_ids, column, 'KO (KEGGCharter)', 'ko', desc=base_desc % (len(ids), 'KEGG IDs', 'KOs'))
        new_ids = pd.merge(new_ids, relational, left_on=column, right_on=1)  # mcj:MCON_3003;   mcj:MCON_3003
        del new_ids[column]  # mcj:MCON_3003    K07486  mcj:MCON_3003;  mcj:MCON_3003
        del new_ids[1]
        new_ids.columns = ['KO (KEGGCharter)', column]
    elif ids_type == 'ko':
        new_ids = id2id(
            ids, column, 'EC number (KEGGCharter)', 'enzyme', desc=base_desc % (len(ids), 'KOs', 'EC numbers'))
    else:  # ids_type == 'ec':
        ids = [f'ec:{ide}' for ide in ids]
        new_ids = id2id(ids, column, 'KO (KEGGCharter)', 'ko', desc=base_desc % (len(ids), 'EC numbers', 'KOs'))
    return pd.merge(data, new_ids, on=column, how='left')


def get_cross_references(data, kegg_column=None, ko_column=None, ec_column=None):
    # KEGG ID to KO -> if KO column is not set, KEGGCharter will get them through the KEGG API
    if kegg_column:
        data = ids_interconversion(data, column=kegg_column, ids_type='kegg')
        data = ids_interconversion(data, column='KO (KEGGCharter)', ids_type='ko')
    if ko_column:
        data = ids_interconversion(data, column=ko_column, ids_type='ko')
        data = ids_interconversion(data, column='EC number (KEGGCharter)', ids_type='ec')
    if ec_column:
        data = ids_interconversion(data, column=ec_column, ids_type='ec')
        data = ids_interconversion(data, column='KO (KEGGCharter)', ids_type='ko')
    if not (kegg_column or ko_column or ec_column):
        exit('Need to specify a column with either KEGG IDs, KOs or EC numbers!')
    return data


def condense_data(data, main_column):
    onlykos = data[data['KO (KEGGCharter)'].notnull() & (data['EC number (KEGGCharter)'].isnull())]
    onlykos = onlykos.groupby(main_column).agg({'KO (KEGGCharter)': lambda x: ','.join(set(x))}).reset_index()
    onlykos['EC number (KEGGCharter)'] = [np.nan] * len(onlykos)
    wecs = data[data['EC number (KEGGCharter)'].notnull()]
    wecs = wecs.groupby(main_column).agg({'KO (KEGGCharter)': lambda x: ','.join(set(x)),
                                          'EC number (KEGGCharter)': lambda x: ','.join(set(x))}).reset_index()
    del data['KO (KEGGCharter)']
    del data['EC number (KEGGCharter)']
    return pd.merge(data, pd.concat([onlykos, wecs]), on=main_column, how='left').drop_duplicates()


def prepare_data_for_charting(data, mt_cols=None, ko_column='KO (KEGGCharter)', ko_from_uniprot=False):
    nokos = data[data[ko_column].isnull()]
    wkos = data[data[ko_column].notnull()].reset_index(drop=True)
    if ko_from_uniprot:
        # If coming from UniProt ID mapping, KOs will be in the form "KXXXXX;"
        wkos[ko_column] = wkos[ko_column].apply(lambda x: x.rstrip(';'))

    # Expand KOs column if some elements are in the form KO1,KO2,...
    wkos[ko_column] = wkos[ko_column].apply(lambda x: x.split(','))
    if mt_cols is not None:
        for col in mt_cols:
            wkos[col] = wkos[col] / wkos[ko_column].apply(lambda x: len(x))
    wkos = expand_by_list_column(wkos, column=ko_column)
    data = pd.concat([wkos, nokos])
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


def write_kgml(mmap, output, organism='ko'):
    data = kegg_get(f"{organism}{mmap}", "kgml").read()
    with open(output, 'w') as f:
        if len(data) > 1:
            f.write(data)
            return KGML_parser.read(data)
    return None


def write_kgmls(mmaps, out_dir, max_tries=3, org='ko'):
    print(f'[{len(mmaps)}] maps inputted for org [{org}]')
    maps_done = [filename.split('/')[-1].rstrip('.xml').lstrip(org) for filename in glob(f'{out_dir}/{org}*.xml')]
    mmap_to_orthologs = {
        name: [orth.id for orth in KGML_parser.read(open(f'{out_dir}/{org}{name}.xml')).orthologs]
        for name in maps_done}  # maps already done will have their orthologs put in
    print(f'[{len(maps_done)}] KGMLs already obtained for org [{org}]')
    mmaps = [map for map in mmaps if map not in maps_done]
    i = 1
    for mmap in mmaps:
        print(f'[{i}/{len(mmaps)}] Obtaining KGML for map [{org}{mmap}]')
        tries = 0
        done = False
        while tries < max_tries and not done:
            # try:
            orthologs = [orth.id for orth in write_kgml(mmap, f'{out_dir}/{org}{mmap}.xml', organism=org).orthologs]
            mmap_to_orthologs[mmap] = orthologs
            done = True
            print(f'[{org}{mmap}]: success')
            '''
            except:
                if os.path.isfile(f'{out_dir}/{org}{mmap}.xml'):
                    os.remove(f'{out_dir}/{org}{mmap}.xml')
                tries += 1
                print(f'[{org}{mmap}]: failure')
            '''
        i += 1
    return mmap_to_orthologs


def set_text_boxes_kgml(kgml_filename, desc=''):
    handler = KGML_parser.read(open(kgml_filename))
    # Set text in boxes to EC numbers
    with open(kgml_filename.replace('xml', 'csv'), 'w') as f:
        for ortholog_rec in tqdm(handler.orthologs, desc=desc):
            lines = list()
            kos = ortholog_rec.name.split()
            lines += kegg_link("enzyme", kos).read().split('\n')
            ecs = [line.split('\t')[1] for line in lines if len(line) > 0]
            if len(ecs) > 0:
                f.write(f'{",".join(ecs)}\n')
            else:
                f.write(f'{",".join(kos)}\n')


# TODO - deprecated?
def set_text_boxes_kgmls(mmaps, out_dir, max_tries=3, org='ko'):
    maps_done = [filename.split('/')[-1].rstrip('.csv') for filename in glob(f'{out_dir}/{org}*.csv')]
    print(f"[{len(maps_done)}] maps already have boxes' text set")
    mmaps = [map for map in mmaps if map not in maps_done]
    i = 1
    for mmap in mmaps:
        tries = 0
        done = False
        while tries < max_tries and not done:
            try:
                set_text_boxes_kgml(
                    f'{out_dir}/{org}{mmap}.xml', desc=f"[{i}/{len(mmaps)}] Getting boxes' labels for map [{mmap}]")
                done = True
            except:
                print(f'Failed for map. Attempt: {tries + 1}')
                if os.path.isfile(f'{out_dir}/map{mmap}.csv'):
                    os.remove(f'{out_dir}/map{mmap}.csv')
                tries += 1
                sleep(10)
        i += 1


def taxon2prefix(taxon_name, organism_df):
    """
    :param taxon_name: str - e.g. Pongo abelii (Sumatran orangutan)
    :param organism_df: pandas.DataFrame - organism file, index = taxa names, column names = KEGG prefix, ...
    :returns str - KEGG prefix of taxon name
    """
    if taxon_name == 'ko':
        return 'ko'
    if taxon_name in organism_df.index:  # exactly the same
        return organism_df.loc[taxon_name][0]
    if taxon_name.split(' (')[0] in organism_df.index:  # Homo sapiens (human) -> Homo sapiens
        return organism_df.loc[taxon_name.split(' (')[0]][0]
    df = organism_df.reset_index()
    possible_prefixes = df[df.name.str.contains(taxon_name)].prefix.tolist()
    if len(possible_prefixes) > 0:
        return df[df.name.str.contains(taxon_name)].prefix.tolist()[0]  # select the first strain
    print(f'[{taxon_name}] was not found in taxon to KEGG prefix conversion!')
    return None


def get_taxon_maps(kegg_prefix):
    if kegg_prefix is None:
        return []
    df = pd.read_csv(StringIO(kegg_list("pathway", kegg_prefix).read()), sep='\t', header=None)
    return df[0].apply(lambda x: x.split(kegg_prefix)[1]).tolist()


def parse_organism(file):
    return pd.read_csv(file, sep='\t', usecols=[1, 2], header=None, index_col=1, names=['prefix', 'name'])


def download_resources(data, resources_directory, taxa_column, metabolic_maps):
    download_organism(resources_directory)
    taxa = ['ko'] + list(set(data[data[taxa_column].notnull()][taxa_column]))
    taxa_df = parse_organism(f'{resources_directory}/organism')
    i = 1
    taxon_to_mmap_to_orthologs = dict()  # {'Keratinibaculum paraultunense' : {'00190': ['1', '2']}}
    for taxon in taxa:
        print(f'[{i}/{len(taxa)}] Getting information for taxon [{taxon}]')
        kegg_prefix = taxon2prefix(taxon, taxa_df)
        if kegg_prefix is not None:
            taxon_mmaps = get_taxon_maps(kegg_prefix)
            taxon_mmaps = [mmap for mmap in taxon_mmaps if mmap in metabolic_maps]  # select only inputted maps
            taxon_to_mmap_to_orthologs[taxon] = write_kgmls(
                taxon_mmaps, resources_directory, org=kegg_prefix)
        else:
            taxon_to_mmap_to_orthologs[taxon] = dict()
        i += 1
    with open(f'{resources_directory}/taxon_to_mmap_to_orthologs.json', 'w') as f:
        json.dump(taxon_to_mmap_to_orthologs, f)
    return taxon_to_mmap_to_orthologs


def get_mmaps2taxa(taxon_to_mmap_to_orthologs):
    mmaps2taxa = dict()
    for org, mmaps2orthologs in taxon_to_mmap_to_orthologs.items():
        for mmap in mmaps2orthologs.keys():
            if mmap in mmaps2taxa.keys():
                mmaps2taxa[mmap].append(org)
            else:
                mmaps2taxa[mmap] = [org]
    return mmaps2taxa


def chart_map(
        kgml_filename, ec_list, data, taxon_to_mmap_to_orthologs, mmaps2taxa, output=None, ko_column=None,
        taxa_column=None, genomic_columns=None, transcriptomic_columns=None, mmap2taxa=None,
        number_of_taxa=10):
    if genomic_columns:  # when not set is None
        mmap = KGML_parser.read(open(kgml_filename))
        kegg_pathway_map = KEGGPathwayMap(pathway=mmap, ec_list=ec_list)
        kegg_pathway_map.genomic_potential_taxa(
            data, genomic_columns, ko_column, taxon_to_mmap_to_orthologs, mmaps2taxa, taxa_column=taxa_column,
            output_basename=f'{output}/potential', number_of_taxa=number_of_taxa)
    if transcriptomic_columns:  # when not set is None
        mmap = KGML_parser.read(open(kgml_filename))
        kegg_pathway_map = KEGGPathwayMap(pathway=mmap, ec_list=ec_list)
        kegg_pathway_map.differential_expression_sample(
            data, transcriptomic_columns, ko_column, taxon_to_mmap_to_orthologs, mmaps2taxa,
            output_basename=f'{output}/differential',
            log=False, taxa_column=taxa_column)
    plt.close()


def get_pathway_and_ec_list(rd, mmap):
    download = True
    if os.path.isfile(f'{rd}/ko{mmap}.xml') and os.path.isfile(f'{rd}/ko{mmap}.csv'):
        pathway = KGML_parser.read(open(f'{rd}/ko{mmap}.xml'))
        with open(f'{rd}/ko{mmap}.csv') as f:
            ec_list = f.read().split('\n')
        if len(pathway.orthologs) == len(ec_list) - 1:  # -1 because of newline at the end
            download = False
        else:
            print(f'Lengths of orthologs in KGML and labels for corresponding boxes do not match for map [ko{mmap}]!')
    else:
        print(f'Some resources were not found for map [ko{mmap}]! Going to download them')
    if download:
        try:
            write_kgml(mmap, f'{rd}/ko{mmap}.xml')
            print(f'Got KGML for map [ko{mmap}]')
            set_text_boxes_kgml(f'{rd}/ko{mmap}.xml', desc=f"Getting boxes' labels for map [ko{mmap}]")
            pathway = KGML_parser.read(open(f'{rd}/ko{mmap}.xml'))
            with open(f'{rd}/ko{mmap}.csv') as f:
                ec_list = f.read().split('\n')
        except:     # may be 404 not found, but also connection timed out, this way everything works
            print(f'Could not download resources for [ko{mmap}]!')
            return None, None
    return pathway, ec_list


def read_input():
    args = get_arguments()
    timed_message('Arguments valid.')
    if args.show_available_maps:
        sys.exit(kegg_metabolic_maps().to_string())
    data = read_input_file(args.file)
    timed_message('Data successfully read.')
    return args, data


def main():
    args, data = read_input()

    if not args.resume:
        data, main_column = further_information(
            data, kegg_column=args.kegg_column, ko_column=args.ko_column, ec_column=args.ec_column)
        data.to_csv(f'{args.output}/KEGGCharter_results.tsv', sep='\t', index=False)
        timed_message(f'Results saved to {args.output}/KEGGCharter_results.tsv')

    ko_column = 'KO (KEGGCharter)'  # TODO - set ko_column to user defined value
    data = prepare_data_for_charting(data, ko_column=ko_column, mt_cols=args.transcriptomic_columns)

    if args.input_quantification:
        data['Quantification (KEGGCharter)'] = [1] * len(data)
        args.genomic_columns = 'Quantification (KEGGCharter)'

    if args.input_taxonomy:
        data['Taxon (KEGGCharter)'] = [args.input_taxonomy] * len(data)
        args.taxa_column = 'Taxon (KEGGCharter)'
        args.taxa_list = args.input_taxonomy

    metabolic_maps = args.metabolic_maps.split(',')
    args.genomic_columns = args.genomic_columns.split(',')
    if args.transcriptomic_columns:
        args.transcriptomic_columns = args.transcriptomic_columns.split(',')

    taxon_to_mmap_to_orthologs = download_resources(data, args.resources_directory, args.taxa_column, metabolic_maps)
    mmaps2taxa = get_mmaps2taxa(taxon_to_mmap_to_orthologs)    # '00190': ['Keratinibaculum paraultunense']
    timed_message(f'Creating KEGG Pathway representations for {len(metabolic_maps)} metabolic pathways.')
    for i in range(len(metabolic_maps)):
        pathway, ec_list = get_pathway_and_ec_list(args.resources_directory, metabolic_maps[i])
        if pathway is not None and ec_list is not None and metabolic_maps[i] in mmaps2taxa:
            timed_message(f'[{i + 1}/{len(metabolic_maps)}] {pathway.title}')
            #pathway_handler = KEGGPathwayMap(pathway, ec_list)
            chart_map(
                f'{args.resources_directory}/ko{metabolic_maps[i]}.xml', ec_list, data, taxon_to_mmap_to_orthologs,
                mmaps2taxa, output=args.output,
                ko_column=ko_column, taxa_column=args.taxa_column,
                genomic_columns=args.genomic_columns, transcriptomic_columns=args.transcriptomic_columns,
                number_of_taxa=args.number_of_taxa)
        else:
            print(f'Analysis of map {metabolic_maps[i]} failed!')
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
