#!/usr/bin/env python

import PIL
import argparse
import multiprocessing
import numpy as np
import os
import pandas as pd
import pathlib
import re
import subprocess
import sys
from io import StringIO
from itertools import product
from time import gmtime, strftime
import time
import matplotlib.pyplot as plt
import traceback
from Bio.KEGG.REST import kegg_link, kegg_list, kegg_get
from matplotlib import colors, cm
from progressbar import ProgressBar
import glob
from Bio.KEGG.KGML import KGML_parser

from kegg_pathway_map import KEGGPathwayMap

__version__ = "0.1.2"


def get_arguments():
    parser = argparse.ArgumentParser(
        description="""KEGGCharter - A tool for representing genomic potential and transcriptomic expression into 
        KEGG pathways""",
        epilog="Input file must be specified.")
    parser.add_argument("-o", "--output", type=str, help="Output directory",
                        default='KEGGCharter_results')
    parser.add_argument("--tsv", action="store_true", default=False,
                        help="Results will be outputed in TSV format (and not EXCEL).")
    #parser.add_argument("-t", "--threads", type=str, default=str(multiprocessing.cpu_count() - 2),
    #                    help="Number of threads for reCOGnizer to use. Default is number of CPUs available minus 2.")
    parser.add_argument("-mm", "--metabolic-maps", type=str,
                        help="IDs of metabolic maps to output",
                        default=','.join(KEGGCharter_prokaryotic_maps()))
    parser.add_argument("-gcol", "--genomic-columns", type=str,
                        help="Names of columns with genomic identification")
    parser.add_argument("-tcol", "--transcriptomic-columns", type=str,
                        help="Names of columns with transcriptomics quantification")
    parser.add_argument("-tls", "--taxa-list", type=str,
                        help="List of taxa to represent in genomic potential charts (comma separated)")  # TODO - must be tested
    parser.add_argument("-not", "--number-of-taxa", type=str,
                        help="Number of taxa to represent in genomic potential charts (comma separated)",
                        default='10')
    parser.add_argument("-keggc", "--kegg-column", type=str, help="Column with KEGG IDs.")
    parser.add_argument("-koc", "--ko-column", type=str, help="Column with KOs.")
    parser.add_argument("-ecc", "--ec-column", type=str, help="Column with EC numbers.")
    parser.add_argument("--resume", action="store_true", default=False,
                        help="Data inputed has already been analyzed by KEGGCharter.")
    parser.add_argument('-v', '--version', action='version', version='KEGGCharter ' + __version__)

    required_named = parser.add_argument_group('required named arguments')
    required_named.add_argument("-f", "--file", type=str, required=True,
                                help="TSV or EXCEL table with information to chart")
    parser.add_argument("-tc", "--taxa-column", type=str, required=True,
                        default='Taxonomic lineage(GENUS)',
                        help="Column with the taxa designations to represent with KEGGChart")  # TODO - test this argument without UniProt shenanigans

    special_functions = parser.add_argument_group('Special functions')
    special_functions.add_argument("--show-available-maps",
                                   action="store_true", default=False,
                                   help="""Outputs KEGG maps IDs and descriptions to the
                        console (so you may pick the ones you want!)""")

    args = parser.parse_args()

    args.output = args.output.rstrip('/')

    for directory in [args.output]:
        if not os.path.isdir(directory):
            pathlib.Path(directory).mkdir(parents=True, exist_ok=True)
            print('Created ' + directory)

    return args


def timed_message(message):
    print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ': ' + message)


def read_input(file):
    if file.endswith('.xlsx'):
        return pd.read_excel(file)
    return pd.read_csv(file, sep='\t')


def expand_by_list_column(df, column='Pathway'):
    if len(df) == 0:
        return pd.DataFrame()
    lens = [len(item) for item in df[column]]
    dictionary = dict()
    for col in df.columns:
        dictionary[col] = np.repeat(df[col].values, lens)
    dictionary[column] = np.concatenate(df[column].values)
    return pd.DataFrame(dictionary)

    # Functions that deal with taxonomies


def most_abundant_taxa(data, columns, taxa_column, number_of_taxa=10):
    """
    Calculates top genus from samples
    :param columns: list of mg columns to consider for quantification of genus abundance
    :param number_of_taxa: number of top genus to return
    :return: list of top genus
    """
    data = data.groupby(taxa_column)[columns].sum().sum(axis=1).sort_values(ascending=False)
    if number_of_taxa > len(data.index.tolist()):
        number_of_taxa = len(data.index.tolist())
    return data.index.tolist()[:number_of_taxa]


def taxa_colors(hex_values=None, ncolor=1):
    """
    Creates list of hex colors to be used, using matplotlib or using custom colors
    :param hex_values: list of hex colors
    :param ncolor: int indicating the ammount of hex colors that should be created
    :return: returns list with hex color codes
    """
    if not hex_values:  # if no colors are given creates a list of hex colors with ncolor from matplotlib discrete colormaps
        color_scheme = (cm.get_cmap('Pastel2', 8) if ncolor <= 8
                        else cm.get_cmap("Set3", 12) if ncolor <= 12
        else cm.get_cmap("rainbow", ncolor))  # if ncolor > 12 a continuous colormap is used instead
        return [colors.to_hex(color_scheme(i)) for i in range(ncolor)]

    for hex_value in hex_values:
        if not re.search(r'^#(?:[0-9a-fA-F]{3}){1,2}$', hex_value):
            sys.exit(Exception("Colors aren't valid hex codes"))
    return hex_values  # has validated hex values and returns the original list


# Conversion functions
def keggid2ko(kegg_ids, kegg_column, step=150):
    """
    Converts KEGG_ID genes to Ortholog KO ID from KEGG
    :param KEGG_ID: (list) - KEGG ID genes
    :param step: (int) - will convert "step" KEGG IDs at a time
    :return: (list) - (list,list) - KEGG ID genes converted and ko IDs
    """
    result = list()
    pbar = ProgressBar()
    for i in pbar(range(0, len(kegg_ids) - step, step)):
        try:
            result += kegg_link("ko", kegg_ids[i:i + step]).read().split("\n")[:-1]
        except:
            print('KEGG ID to KO broke at index ' + str(i))
            result = [[part[0] + ';', part[1].strip('ko:')] for part in
                      [relation.split('\t') for relation in result]]
            return pd.DataFrame(result, columns=[kegg_column, 'KO (KEGG Charter)'])
    result += kegg_link("ko", kegg_ids[len(kegg_ids) - step:]).read().split("\n")[:-1]
    result = [[part[0] + ';', part[1].strip('ko:')] for part in
              [relation.split('\t') for relation in result]]
    return pd.DataFrame(result, columns=[kegg_column, 'KO (KEGG Charter)'])


def ko2ec(kos, ko_column_name, step=150):
    """
    Converts KOs to EC numbers
    :param kos: list of kegg orthologs
    :return: dic associating ortholog kegg id with list
    of assotiated EC numbers
    """
    result = list()
    pbar = ProgressBar()
    for i in pbar(range(0, len(kos), step)):
        try:
            result += kegg_link("enzyme", kos[i:i + step]).read().split("\n")[:-1]
        except:
            print('KO to EC number broke at index ' + str(i))
            result = [relation.split('\t') for relation in result]
            return list(map(list, zip(*result)))
    result += kegg_link("enzyme", kos[len(kos) - step:]).read().split("\n")[:-1]
    result = [[part[0].strip('ko:'), part[1].upper()] for part in
              [relation.split('\t') for relation in result]]
    return pd.DataFrame(result, columns=[ko_column_name, 'EC number (KEGG Charter)'])


def ec2ko(ecnumbers, ec_column_name, step=150):
    """
    Converts KOs to EC numbers
    :param kos: list of kegg orthologs
    :return: dic associating ortholog kegg id with list
    of assotiated EC numbers
    """
    result = list()
    pbar = ProgressBar()
    for i in pbar(range(0, len(ecnumbers), step)):
        try:
            result += kegg_link("ko", ecnumbers[i:i + step]).read().split("\n")[:-1]
        except:
            print('KO to EC number broke at index ' + str(i))
            result = [relation.split('\t') for relation in result]
            return list(map(list, zip(*result)))
    result += kegg_link("enzyme", ecnumbers[len(ecnumbers) - step:]).read().split("\n")[:-1]
    result = [[part[0], part[1].strip('ko:')] for part in
              [relation.split('\t') for relation in result[:-1]]]
    print(pd.DataFrame(result, columns=[ec_column_name, 'KO (KEGG Charter)']))
    return pd.DataFrame(result, columns=[ec_column_name, 'KO (KEGG Charter)'])


def get_cross_references(data, kegg_column=None, ko_column=None, ec_column=None):
    # KEGG ID to KO -> if KO column is not set, KEGGCharter will get them through the KEGG API
    if not ko_column:
        if kegg_column:
            kegg_ids = data[data[kegg_column].notnull()][kegg_column]
            kegg_ids = [ide.split(';')[0] for ide in
                        kegg_ids]  # TODO - sometimes Kegg IDs come as mfc:BRM9_0145;mfi:DSM1535_1468; (e.g. from UniProt IDs mapping). Should be no problem since both IDs should always be same function
            timed_message('Converting {} KEGG IDs to KOs through the KEGG API.'.format(len(kegg_ids)))
            kos = keggid2ko(kegg_ids, kegg_column)
            data = pd.merge(data, kos, on=kegg_column, how='left')
        elif ec_column:
            ec_numbers = data[data[ec_column].notnull()][ec_column]
            ec_numbers = ['ec:{}'.format(ide) for ide in ec_numbers]
            timed_message('Converting {} EC numbers to KOs through the KEGG API.'.format(len(ec_numbers)))
            kos = ec2ko(ec_numbers, ec_column)
            ec_column = ec_column
            kos[ec_column] = [ide[3:] for ide in kos[ec_column]]
            data = pd.merge(data, kos, on=ec_column, how='left')
        else:
            exit('Need to specify a column with either KEGG IDs, KOs or EC numbers!')
        ko_column = 'KO (KEGG Charter)'
    else:
        ko_column = ko_column
    data[ko_column] = data[ko_column].apply(lambda x: x.rstrip(';') if type(
        x) != float else x)  # If coming from UniProt ID mapping, KOs will be in the form KXXXXX;

    # Expand KOs column if some elements are in the form KO1,KO2,...
    data[ko_column] = [ko.split(',') if type(ko) != float else ko for ko in data[ko_column]]
    wkos = data[data[ko_column].notnull()]
    nokos = data[data[ko_column].isnull()]
    wkos = expand_by_list_column(wkos, column=ko_column)
    data = pd.concat([wkos, nokos])

    # KO to EC number
    kos = data[data[ko_column].notnull()][ko_column].tolist()
    timed_message('Retrieving EC numbers from {} KOs.'.format(len(kos)))
    ecs = ko2ec(kos, ko_column)
    data = pd.merge(data, ecs, on=ko_column, how='left')

    notnull = data[data['EC number (KEGG Charter)'].notnull()]
    notnull = notnull.groupby('Entry').agg(
        {'KO (KEGG Charter)': lambda x: ','.join(set(x)), 'EC number (KEGG Charter)': lambda x: ','.join(set(x))}
    ).reset_index()
    results = data[data.columns[:-2]].drop_duplicates()
    results = pd.merge(results, notnull, on='Entry', how='left')
    return results


# Get metabolic maps from KEGG Pathway
def KEGGCharter_prokaryotic_maps(file=sys.path[0] + '/KEGGCharter_prokaryotic_maps.txt'):
    return open(file).read().split('\n')


def KEGG_metabolic_maps():
    """
    Creates a dic with all specific kegg pathway maps and their description
    :return: pandas.DataFrame with Map ID as index and maps names as
    sole column
    """
    maps = pd.read_csv(StringIO(kegg_list("pathway").read()), sep='\t', names=['Map ID', 'Description'])
    maps['Map ID'] = [ide.split('map')[-1] for ide in maps['Map ID']]
    return maps


def add_blank_space(image_pil, width, height, image_mode='RGB'):
    """
    Resizes an image with white background, keeping image size ratio
    :param image_pil: PIL.Image - image to be resized
    :param width: int - width of final image
    :param height: int - heigth of final image
    :param image_mode: str - image mode of image (RGBA, RGB, ...)
    """
    ratio_w = width / image_pil.width
    ratio_h = height / image_pil.height
    if ratio_w < ratio_h:
        # It must be fixed by width
        resize_width = width
        resize_height = round(ratio_w * image_pil.height)
    else:
        # Fixed by height
        resize_width = round(ratio_h * image_pil.width)
        resize_height = height
    image_resize = image_pil.resize((resize_width, resize_height),
                                    PIL.Image.ANTIALIAS)
    background = PIL.Image.new('RGB', (width, height), (255, 255, 255, 255))
    offset = (round((width - resize_width) / 2),
              round((height - resize_height) / 2))
    background.paste(image_resize, offset)
    return background.convert(image_mode)


def resize_image(image_pil, ratio=None, width=None, height=None):
    """
    Resizes an image with alteration to image size ratio
    :param ratio: int - ratio of resize - how bigger or smaller will the output be?
    :param image_pil: PIL.Image - image to be resized
    :param width: int - width of final image
    :param height: int - heigth of final image
    """
    if ratio:
        return image_pil.resize((image_pil.width * ratio,
                                 image_pil.height * ratio), PIL.Image.ANTIALIAS)
    elif width and height:
        return image_pil.resize((width, height), PIL.Image.ANTIALIAS)
    else:
        return None


def pdf2png(pdf_filename):
    """
    Converts a pdf file to a png file, RGB format - name changes, .pdf to .png
    :param pdf_filename: str - filename of PDF file
    """
    bash_command = 'pdftoppm {} {} -png'.format(pdf_filename, pdf_filename.split('.pdf')[0])
    subprocess.run(bash_command.split())
    os.rename(pdf_filename.replace('.pdf', '-1.png'), pdf_filename.replace('.pdf', '.png'))


def add_legend(kegg_map_file, legend_file, output):
    """
    Merges the two files - KEGG metabolic map and respective legend - into
    one file file
    :param kegg_map_file: str - filename of PDF kegg metabolic map
    :param legend_file: str - filename of PNG legend
    """
    pdf2png(kegg_map_file)
    imgs = [PIL.Image.open(file) for file in
            [kegg_map_file.replace('.pdf', '.png'), legend_file]]
    imgs[0] = imgs[0].convert(
        'RGB')  # KEGG Maps are converted to RGB by pdftoppm, dunno if converting to RGBA adds any transparency
    imgs[1] = resize_image(imgs[1], ratio=5)
    imgs[1] = add_blank_space(imgs[1], imgs[1].width, imgs[0].height)
    imgs_comb = np.hstack([np.asarray(i) for i in imgs])

    # save that beautiful picture
    imgs_comb = PIL.Image.fromarray(imgs_comb)
    imgs_comb.save(output)
    for file in [kegg_map_file, kegg_map_file.replace('.pdf', '.png'), legend_file]:
        os.remove(file)


def create_potential_legend(colors, labels, filename):
    """
    Draws the color to taxa labels of genomic potential representations
    :param colors: list - list of colors of the different taxa
    :param labels: list - list of taxa corresponding to the colors
    :param filename: string - filename to output
    """
    f = lambda m, c: plt.plot([], [], marker=m, color=c, ls="none")[0]
    handles = [f("s", color) for color in colors]
    legend = plt.legend(handles, labels, loc=3, framealpha=1, frameon=True)
    fig = legend.figure
    fig.canvas.draw()
    # The bbox manipulation removes the axis
    bbox = legend.get_window_extent()
    bbox = bbox.from_extents(*(bbox.extents + np.array([-2, -2, 2, 2])))
    bbox = bbox.transformed(fig.dpi_scale_trans.inverted())

    fig.savefig(filename, dpi="figure", bbox_inches=bbox)


def genomic_potential_taxa(kegg_pathway_map, data, samples, dic_colors, ko_column,
                           taxa=None, taxa_column='Taxonomic lineage (GENUS)',
                           output_basename=None, maxshared=10):
    """
    Represents the genomic potential of the dataset for a certain taxa level,
    by coloring each taxon with a unique color
    :param data: pandas.DataFrame with data already processed by KEGGPathway
    :param samples: list of str column names of the dataset correspoding to
    expression values
    :param genera: list of genus to represent
    :param number_of_taxa: int representing the number of diferent taxa to
    be represented in the maps, in case the taxa are not specified
    :param level_of_taxa: str - taxonomic level to represent - SPECIES,
    SUPERKINGDOM, ...
    :param output_basename: str - basename for map outputs
    :param maxshared: int - maximum number of different taxa to represent
    in a single map box
    """
    # for every taxon, check all boxes it is in, and save that info to box2taxon
    box2taxon = dict()
    for taxon in dic_colors.keys():
        df = data[data[taxa_column] == taxon][samples + [ko_column]]
        df = df[df.any(axis=1)]
        for ortholog in df[ko_column]:
            if ortholog in kegg_pathway_map.ko_boxes.keys():
                for box in kegg_pathway_map.ko_boxes[ortholog]:
                    if box in box2taxon.keys():
                        if taxon not in box2taxon[box]:
                            box2taxon[box].append(taxon)
                    else:
                        box2taxon[box] = [taxon]

    # for every box with KOs identified from the most abundant taxa, sub-boxes are created with colours of the corresponding taxa
    kegg_pathway_map.pathway_box_list(box2taxon, dic_colors)

    # boxes with KOs identified but not from the most abundant taxa are still identified
    df = data[samples + [ko_column]]
    df = df[df.any(axis=1)]
    grey_boxes = list()
    for ortholog in df[ko_column]:
        if ortholog in kegg_pathway_map.ko_boxes.keys():
            for box in kegg_pathway_map.ko_boxes[ortholog]:
                if box not in box2taxon.keys() and box not in grey_boxes:
                    grey_boxes.append(box)
    kegg_pathway_map.grey_boxes(grey_boxes)  # TODO - boxes should have EC number as name

    name = kegg_pathway_map.pathway.name.split(':')[-1]
    name_pdf = '{}_{}.pdf'.format(output_basename, name)
    kegg_pathway_map.to_pdf(name_pdf)

    # a new taxon: "Other taxa"
    if len(grey_boxes) > 0:
        dic_colors["Other taxa"] = "#7c7272"

    # TODO - legend should be ajusted for the maps - otherwise, makes no sense to have one legend for each map - they all become the same, except for "Other taxa"
    create_potential_legend(dic_colors.values(), dic_colors.keys(),
                            name_pdf.replace('.pdf', '_legend.png'))

    add_legend(name_pdf, name_pdf.replace('.pdf', '_legend.png'),
               name_pdf.replace(name + '.pdf', kegg_pathway_map.pathway.title.replace('/', '|') + '.png'))


def differential_colorbar(dataframe, filename):
    FIGSIZE = (2, 3)
    mpb = plt.pcolormesh(dataframe, cmap='coolwarm')
    fig, ax = plt.subplots(figsize=FIGSIZE)
    plt.colorbar(mpb, ax=ax)
    ax.remove()
    plt.savefig(filename, bbox_inches='tight')


def differential_expression_sample(kegg_pathway_map, data, samples, ko_column,
                                   output_basename=None, log=True):
    """
    Represents in small heatmaps the expression levels of each sample on the
    dataset present in the given pathway map. The values can be transford to
    a log10 scale
    :param data: pandas.DataFrame with data already processed by KEGGPathway
    :param samples: list - column names of the dataset corresponding to
    expression values
    :param output_folder: string - name of folder to store pdfs
    :param log: bol - convert the expression values to logarithmic scale?
    """
    df = data.groupby(ko_column)[samples + [ko_column]].sum()
    df = df[df.any(axis=1)]
    df['Boxes'] = [kegg_pathway_map.ko_boxes[ko] if ko in
                                                    kegg_pathway_map.ko_boxes.keys() else np.nan for ko in df.index]
    df = df[df['Boxes'].notnull()]
    df = expand_by_list_column(df, column='Boxes')
    if len(df) == 0:
        return 1
    df = df.groupby('Boxes')[samples].sum()

    kegg_pathway_map.pathway_boxes_differential(df, log)

    name = kegg_pathway_map.pathway.name.split(':')[-1]
    name_pdf = '{}_{}.pdf'.format(output_basename, name)
    kegg_pathway_map.to_pdf(name_pdf)

    differential_colorbar(df, name_pdf.replace(".pdf", '_legend.png'))

    add_legend(name_pdf, name_pdf.replace('.pdf', '_legend.png'),
               name_pdf.replace(name + '.pdf', kegg_pathway_map.pathway.title.replace('/', '|') + '.png'))

    return 0


def write_results(data, output_dir, output_type='excel'):
    if output_type == 'excel':
        data.to_excel('{}/KEGGCharter_results.xlsx'.format(output_dir), index=False)
    else:
        data.to_csv('{}/KEGGCharter_results.tsv'.format(output_dir), sep='\t', index=False)


def write_kgml(mmap, out_dir):
    data = kegg_get("ko" + mmap, "kgml").read()
    with open('{}/map{}.xml'.format(out_dir, mmap), 'w') as f:
        if len(data) > 1:
            f.write(data)


def write_kgmls(mmaps, out_dir):
    print('[{}] maps inputted'.format(len(mmaps)))
    maps_done = [filename.split('/map')[-1].rstrip('.xml') for filename in glob.glob('{}/*.xml'.format(out_dir))]
    print('[{}] maps already obtained'.format(len(maps_done)))
    mmaps = [map for map in mmaps if map not in maps_done]

    failed = list()
    i = 0
    for mmap in mmaps:
        print('[{}/{}] Obtaining KGML for [{}]'.format(i, len(mmaps), mmap))
        try:
            write_kgml(mmap, out_dir)
            print('[{}]: success'.format(mmap))
        except:
            failed.append(mmap)
            print('[{}]: failure'.format(mmap))
        i += 1

    with open('{}/failed_maps.txt'.format(out_dir), 'w') as f:
        f.write('\n'.join(failed))


def chart_map(mmap, data, output=None, ko_column=None, ec_column=None, taxa_column=None, dic_colors=None,
              genomic_columns=None, transcriptomic_columns=None):
    timed_message('Handling pathway: {}'.format(mmap.title))
    if genomic_columns:  # when not set is None
        kegg_pathway_map = KEGGPathwayMap(pathway=mmap)
        genomic_potential_taxa(
            kegg_pathway_map, data, genomic_columns, dic_colors, ko_column, taxa_column=taxa_column,
            output_basename=output + '/potential')
    if transcriptomic_columns:  # when not set is None
        kegg_pathway_map = KEGGPathwayMap(pathway=mmap)
        differential_expression_sample(
            kegg_pathway_map, data, transcriptomic_columns, ko_column, output_basename=output + '/differential',
            log=False)


def main():
    # Read input
    args = get_arguments()
    timed_message('Arguments valid.')

    if args.show_available_maps:
        sys.exit(KEGG_metabolic_maps().to_string())
    data = read_input(args.file)
    timed_message('Data successfully read.')

    if not args.resume:
        results = get_cross_references(
            data, kegg_column=args.kegg_column, ko_column=args.ko_column, ec_column=args.ec_column)
        write_results(results, args.output, output_type=('tsv' if args.tsv else 'excel'))
        timed_message('Results saved to {}/KEGGCharter_results.{}'.format(args.output, 'tsv' if args.tsv else 'xlsx'))
    ko_column = args.ko_column if args.ko_column else 'KO (KEGG Charter)'
    ec_column = args.ec_column if args.ec_column else 'EC number (KEGG Charter)'

    # Begin dat chart magic
    metabolic_maps = args.metabolic_maps.split(',')
    timed_message('Creating KEGG Pathway representations for {} metabolic pathways.'.format(len(
        metabolic_maps)))
    write_kgmls(metabolic_maps, sys.path[0])

    # Set colours for taxa if MG data is present
    if hasattr(args, 'genomic_columns'):
        args.genomic_columns = args.genomic_columns.split(',')
        if args.taxa_list is None:
            taxa = most_abundant_taxa(data, args.genomic_columns, args.taxa_column,
                                      number_of_taxa=int(args.number_of_taxa))
        else:
            taxa = args.taxa_list.split(',')

        colors = taxa_colors(ncolor=len(taxa))
        dic_colors = {taxa[i]: colors[i] for i in range(len(taxa))}
    else:
        dic_colors = None

    if args.transcriptomic_columns:
        args.transcriptomic_columns = args.transcriptomic_columns.split(',')

    i = 1
    failed_maps = list()
    differential_no_kos = list()

    # For each metabolic map, will chart genomic potential and differential expression if MG and MT data are available, respectively

    pathways = [KGML_parser.read(open('{}/map{}.xml'.format(sys.path[0], mmap))) for mmap in metabolic_maps]
    for pathway in pathways:
        #try:
        chart_map(pathway, data, args.output, ko_column, ec_column,
                         args.taxa_column, dic_colors, args.genomic_columns,
                         args.transcriptomic_columns)
        #except:
        #    failed_maps.append(mmap)

    '''
    
    with multiprocessing.Pool() as p:
        p.starmap(chart_map, [(handler, data, args.output, ko_column, ec_column, args.taxa_column, dic_colors,
                              args.genomic_columns, args.transcriptomic_columns,
                              '{}/failed_maps.txt'.format(args.output)) for handler in pathways])
    '''
    failed_maps_filename = args.output + '/failed_maps.txt'
    no_kos_filename = args.output + '/no_kos.txt'
    print('{} maps could not be loaded. You can see which ones at {}'.format(
        len(failed_maps), failed_maps_filename))
    open(failed_maps_filename, 'w').write('\n'.join(failed_maps))
    print('{} maps had no KOs attributed to them. You can see which ones at {}'.format(
        len(differential_no_kos), no_kos_filename))
    open(no_kos_filename, 'w').write('\n'.join(differential_no_kos))


if __name__ == '__main__':
    start_time = time.time()
    main()
    print('KEGGCharter analysis finished in {}'.format(strftime("%Hh%Mm%Ss", gmtime(time.time() - start_time))))
