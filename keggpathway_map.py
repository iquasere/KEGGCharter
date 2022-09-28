#!/usr/bin/env python

from Bio.KEGG.KGML import KGML_pathway
from Bio.Graphics.KGML_vis import KGMLCanvas
from PIL import Image
import numpy as np
import os
from subprocess import run
from matplotlib import pyplot as plt, colors, cm
import pandas as pd
from re import search
import sys
import time


def set_bgcolor(pathway_element, color):
    """
    Sets graphic element background color
    :param pathway_element: kegg pathway xml element object
    :param color: color to be used in rgb
    """
    pathway_element.graphics[0].bgcolor = color


def set_fgcolor(pathway_element, color):
    """
    Sets graphic element contour color
    :param pathway_element: kegg pathway xml element object
    :param color:  color to be used in rgb
    """
    pathway_element.graphics[0].fgcolor = color


def conv_value_rgb(value, colormap, norm):
    """
    Normalizes values in a vector and transforms it into corresponding
    hex color using given colormap
    :param value: numpy vector from dataframe apply function with expression
    values
    :param colormap: matplotlib colormap object
    :param norm: matplotlib Normalize or LogNorm object
    :return: returns hex colors in vector format
    """
    return value.apply(norm).apply(colormap)


def conv_rgb_hex(rgb):
    """
    converts rgb into hex color code in vector
    :param rgb: rgb value in vector resulting from pandas dataframe apply
    :return: vector with converted hex color code
    """
    return rgb.apply(colors.to_hex)


def create_tile_box(record):
    """
    Create box graphical element in pathway to draw the box countour and
    give the correct name
    :param record: graphical element to be created
    """
    newrecord = KGML_pathway.Graphics(record)
    newrecord.name = record.graphics[0].name
    newrecord.type = "rectangle"
    newrecord.width = record.graphics[0].width
    newrecord.height = record.graphics[0].height
    newrecord.y = record.graphics[0].y
    newrecord.x = record.graphics[0].x
    newrecord.bgcolor = "#FFFFFF00"
    newrecord.fgcolor = "#000000"
    record.graphics.append(newrecord)
    record.graphics[0].bgcolor = "#FFFFFF00"
    record.graphics[0].fgcolor = "#FFFFFF00"
    record.graphics[0].name = ""


def create_box_heatmap(rec_old, nrboxes, i, paired=True):
    """
    Helper function for creating heatmap, draws one expression value in its
    correct position on the bigger parent box
    :param rec_old: graphical element object to be used as reference
    :param nrboxes: int nr of boxes to be drawed
    :param i: int internal number of movements of the box given by the for loop
    :param paired: boolean is number of samples a pair?
    :return: graphical element object
    """
    if rec_old.graphics[0].width is None:  # it happens, guess it is because it has no KOs
        return 1
    movement_steps = rec_old.graphics[0].width / (nrboxes * (2 if paired else 1))
    newrecord = KGML_pathway.Graphics(rec_old)
    newrecord.name = ""
    newrecord.type = "rectangle"
    adjustment_factor = 1.3 if nrboxes > 2 else 1.1 if nrboxes > 1 else 1  # sub-boxes width, adjusted by a factor that experimentally fixed well in the representations
    newrecord.width = movement_steps * adjustment_factor * (2 if paired else 1)
    newrecord.height = rec_old.graphics[0].height
    newrecord.y = rec_old.graphics[0].y
    newrecord.x = (i * movement_steps) + rec_old.graphics[0].x
    newrecord.fgcolor = "#FFFFFF00"
    return newrecord


def pdf2png(pdf_filename):
    """
    Converts a pdf file to a png file, RGB format - name changes, .pdf to .png
    :param pdf_filename: str - filename of PDF file
    """
    bash_command = f'pdftoppm {pdf_filename} {pdf_filename.split(".pdf")[0]} -png'
    run(bash_command.split())
    os.rename(pdf_filename.replace('.pdf', '-1.png'), pdf_filename.replace('.pdf', '.png'))


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
                                 image_pil.height * ratio), Image.ANTIALIAS)
    elif width and height:
        return image_pil.resize((width, height), Image.ANTIALIAS)
    else:
        return None


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
                                    Image.ANTIALIAS)
    background = Image.new('RGB', (width, height), (255, 255, 255, 255))
    offset = (round((width - resize_width) / 2),
              round((height - resize_height) / 2))
    background.paste(image_resize, offset)
    return background.convert(image_mode)


def expand_by_list_column(df, column='Pathway'):
    if len(df) == 0:
        return pd.DataFrame()
    lens = [len(item) for item in df[column]]
    dictionary = dict()
    for col in df.columns:
        dictionary[col] = np.repeat(df[col].values, lens)
    dictionary[column] = np.concatenate(df[column].values)
    return pd.DataFrame(dictionary)


def taxa_colors(hex_values=None, ncolor=1):
    """
    Creates list of hex colors to be used, using matplotlib or using custom colors
    :param hex_values: list of hex colors
    :param ncolor: int indicating the ammount of hex colors that should be created
    :return: returns list with hex color codes
    """
    if not hex_values:  # if no colors are given creates a list of discrete hex colors
        color_scheme = (
            cm.get_cmap('Pastel2', 8) if ncolor <= 8
            else cm.get_cmap("Set3", 12) if ncolor <= 12
            else cm.get_cmap("rainbow", ncolor))  # if ncolor > 12 a continuous colormap is used instead
        return [colors.to_hex(color_scheme(i)) for i in range(ncolor)]

    for hex_value in hex_values:
        if not search(r'^#(?:[0-9a-fA-F]{3}){1,2}$', hex_value):
            sys.exit(Exception("Colors aren't valid hex codes"))
    return hex_values  # has validated hex values and returns the original list


class KEGGPathwayMap:
    """
    This class retrieves and manipulates KEGG metabolic maps from KEGG Pathway
    """

    def __init__(self, pathway, ec_list):
        """
        Initialize object
        :param pathway: (Bio.KEGG.KGML.KGML_pathway.Pathway)
        :param ec_list: (list) - list from ECs CSV
        """
        self.pathway = pathway
        self.set_pathway(ec_list)

    def __getattr__(self, item):
        m = getattr(self.pathway, item, None)
        if m is None:
            raise AttributeError
        return m

    def set_pathway(self, ec_list):
        """
        Set pathway with KEGG Pathway ID
        """
        self.ko_boxes = {}
        for i in range(len(self.orthologs)):
            set_bgcolor(self.orthologs[i], "#ffffff")  # set all boxes to white
            # self.set_fgcolor(self.pathway.orthologs[i], "#ffffff")             # This might be helpful in the future, if an additional layer of marking is needed
            orthologs_in_box = [ide[3:] for ide in self.orthologs[
                i].name.split()]  # 'ko:K16157 ko:K16158 ko:K16159' -> ['K16157', 'K16158', 'K16159']
            for ortholog in orthologs_in_box:
                if ortholog not in self.ko_boxes.keys():
                    self.ko_boxes[ortholog] = []
                self.ko_boxes[ortholog].append(i)  # {'K16157':[0,13,432], 'K16158':[4,13,545]}

            # Set name as EC number
            ecs = ec_list[i].split(',')
            if len(ecs) > 0:
                self.orthologs[i].graphics[0].name = max(set(ecs), key=ecs.count).split(':')[1]
            else:
                self.orthologs[i].graphics[0].name = orthologs_in_box[0].split(':')[1]

    ############################################################################
    ####                          Operations                                ####
    ############################################################################

    def to_pdf(self, filename, imagemap=True, orthologs=True, compounds=True, maps=True, reactions=True):
        """
        Prints current pathway to PDF file
        :param filename: (str) - PDF filename
        :param imagemap: (bol) - Print imagemap
        :param orthologs: (bol) - Print orthologs
        :param compounds: (bol) - Print compounds
        :param maps: (bol) - Print maps
        :param reactions: (bol) - Print reactions ???
        :return: creates PDF file with current pathway
        """
        # TODO - check reactions parameter
        try:
            KGMLCanvas(self,
                       import_imagemap=imagemap,
                       label_orthologs=orthologs,
                       label_compounds=compounds,
                       label_maps=maps,
                       label_reaction_entries=reactions).draw(filename)
        except PermissionError:
            exit(f'Was forbidden to write to {filename}. Check your permissions for that folder/filename')

    def pathway_box_list(self, taxa_in_box, dic_colors, maxshared=10):
        """
        Represents items in the pathway map
        :param taxa_in_box: dict - {box : list of taxa in box}
        :param dic_colors: list - of colors to assign to taxa
        :param maxshared: int - maximum number of taxa sharing one box
        """
        for box in taxa_in_box.keys():
            nrboxes = len(taxa_in_box[box])
            if nrboxes > maxshared:
                nrboxes = maxshared

            paired = True if nrboxes % 2 == 0 else False
            for i in range(nrboxes):
                newrecord = create_box_heatmap(
                    self.orthologs[box], nrboxes, i * 2 - (nrboxes - 1) if paired else i - int(nrboxes / 2),
                    # if nrboxes = 8, i * 2 - (nrboxes - 1) = -7,-5,-3,-1,1,3,5,7
                    # if nrboxes = 9, i - int(nrboxes / 2) = -4,-3,-2,-1,0,1,2,3,4
                    paired=paired)
                if newrecord != 1:
                    newrecord.bgcolor = dic_colors[taxa_in_box[box][i]]
                    self.orthologs[box].graphics.append(newrecord)
            # TODO - must check more deeply why sometimes width is None
            if self.orthologs[box].graphics[0].width is not None:
                create_tile_box(self.orthologs[box])

    def pathway_boxes_differential(self, dataframe, log=True, colormap="coolwarm"):
        """
        Represents expression values present in a dataframe in the
        pathway map
        :param dataframe: pandas DataFrame with each column representing a sample
        and index corresponding to int list index of the ortholog element in the
        pathway
        :param log: bol providing the option for a log transformation of data
        :param colormap: str representing a costum matplotlib colormap to be used
        """

        if log:
            norm = cm.colors.LogNorm(vmin=dataframe.min().min(), vmax=dataframe.max().max())
        else:
            norm = cm.colors.Normalize(vmin=dataframe.min().min(), vmax=dataframe.max().max())

        colormap = cm.get_cmap(colormap)
        dataframe = dataframe.apply(conv_value_rgb, args=(colormap, norm))  # TODO - Doesn't work if using log
        dataframe = dataframe.apply(conv_rgb_hex)

        dataframe = dataframe[dataframe.columns.tolist()]

        nrboxes = len(dataframe.columns.tolist())  # number of samples

        for box in dataframe.index.tolist():
            colors = dataframe.loc[box].tolist()
            paired = True if nrboxes % 2 == 0 else False
            for i in range(nrboxes):
                newrecord = create_box_heatmap(
                    self.orthologs[box], nrboxes, i * 2 - (nrboxes - 1) if paired else i - int(nrboxes / 2),
                    # if nrboxes = 8, i * 2 - (nrboxes - 1) = -7,-5,-3,-1,1,3,5,7; if nrboxes = 9, i - int(nrboxes / 2) = -4,-3,-2,-1,0,1,2,3,4
                    paired=paired)
                if newrecord != 1:  # TODO - assess why sometimes get 1
                    newrecord.bgcolor = colors[i]
                    self.orthologs[box].graphics.append(newrecord)
            if self.orthologs[box].graphics[
                0].width is not None:  # TODO - should check more deeply why sometimes width is None
                create_tile_box(self.orthologs[box])

    def grey_boxes(self, box_list):
        for i in box_list:
            set_bgcolor(self.orthologs[i], "#7c7272")
            set_fgcolor(self.orthologs[i], "#ffffff")

    def create_potential_legend(self, map_colors, labels, filename):
        """
        Draws the color to taxa labels of genomic potential representations
        :param map_colors: list - list of colors of the different taxa
        :param labels: list - list of taxa corresponding to the colors
        :param filename: string - filename to output
        """
        handles = [plt.plot([], [], marker="s", color=color, ls="none")[0] for color in map_colors]
        legend = plt.legend(handles, labels, loc=3, framealpha=1, frameon=True)
        fig = legend.figure
        fig.canvas.draw()
        # The bbox manipulation removes the axis
        bbox = legend.get_window_extent()
        bbox = bbox.from_extents(*(bbox.extents + np.array([-2, -2, 2, 2])))
        bbox = bbox.transformed(fig.dpi_scale_trans.inverted())

        fig.savefig(filename, dpi="figure", bbox_inches=bbox)

    def most_abundant_taxa(self, data, columns, taxa_column, number_of_taxa=10):
        """
        Calculates top genus from samples
        :param data: data
        :param columns: list of mg columns to consider for quantification of genus abundance
        :param taxa_column: da column with da taxa
        :param number_of_taxa: number of top genus to return
        :return: list of top genus
        """
        data = data.groupby(taxa_column)[columns].sum().sum(axis=1).sort_values(ascending=False)
        if number_of_taxa > len(data.index.tolist()):
            number_of_taxa = len(data.index.tolist())
        return data.index.tolist()[:number_of_taxa]

    def genomic_potential_taxa(
            self, data, samples, ko_column, taxon_to_mmap_to_orthologs, mmaps2taxa,
            taxa_column='Taxonomic lineage (GENUS)', output_basename=None, number_of_taxa=10, grey_taxa='Other taxa'):
        """
        Represents the genomic potential of the dataset for a certain taxa level,
        by coloring each taxon with a unique color
        :param data: pandas.DataFrame with data already processed by KEGGPathway
        :param samples: list of str column names of the dataset correspoding to
        expression values
        :param taxon_to_mmap_to_orthologs: dict - {'Keratinibaculum paraultunense' : {'00190': ['1', '2']}}
        :param mmaps2taxa: dict - of taxa to color
        :param ko_column: str - column with KOs
        :param taxa_column: str - column with taxonomic classification
        :param output_basename: str - basename for map outputs
        :param number_of_taxa: int - number of most abundant taxa to represent in each map
        :param grey_taxa: str - name of taxa to represent in grey
        """
        box2taxon = {}
        if mmaps2taxa is not None:
            # for every taxon, check all boxes it is in, and save that info to box2taxon
            data = data[data[taxa_column].isin(mmaps2taxa[self.name.split('ko')[1]]) &
                        data[ko_column].isin(self.ko_boxes.keys())]
            taxa = self.most_abundant_taxa(data, samples, taxa_column, number_of_taxa=number_of_taxa)
            taxonomy_colors = taxa_colors(ncolor=len(taxa))
            dic_colors = {taxa[i]: taxonomy_colors[i] for i in range(len(taxa))}
            for taxon in dic_colors.keys():
                df = data[data[taxa_column] == taxon][samples + [ko_column]]
                df = df[df.any(axis=1)]
                for ortholog in df[ko_column]:
                    if ortholog in self.ko_boxes.keys():
                        for box in self.ko_boxes[ortholog]:
                            if box in box2taxon.keys() and (
                                    box in taxon_to_mmap_to_orthologs[taxon][self.name.split('ko')[1]]):
                                if taxon not in box2taxon[box]:
                                    box2taxon[box].append(taxon)
                            else:
                                box2taxon[box] = [taxon]
            # for every box with KOs identified from the most abundant taxa, sub-boxes are created with colours of the
            # corresponding taxa
            self.pathway_box_list(box2taxon, dic_colors)
        else:
            dic_colors = {}
        # boxes with KOs identified but not from the most abundant taxa are still identified
        df = data[(data[taxa_column].isin(mmaps2taxa[self.name.split('ko')[1]]) &
                   data[ko_column].isin(self.ko_boxes.keys())) & ~data[taxa_column].isin(taxa)]
        df = df[df.any(axis=1)]
        for ortholog in df[ko_column]:
            if ortholog in self.ko_boxes.keys():
                dic_colors[grey_taxa] = "#7c7272"
                for box in self.ko_boxes[ortholog]:
                    if box not in box2taxon.keys():
                        box2taxon[box].append(grey_taxa)
                    else:
                        box2taxon[box] = [grey_taxa]

        name = self.name.split(':')[-1]
        name_pdf = f'{output_basename}_{name}.pdf'
        self.to_pdf(name_pdf)

        # TODO - legend should be ajusted for the maps - otherwise, makes no sense to have one legend for each map -
        #  they all become the same, except for "Other taxa"
        self.create_potential_legend(dic_colors.values(), dic_colors.keys(), name_pdf.replace('.pdf', '_legend.png'))

        self.add_legend(
            name_pdf, name_pdf.replace('.pdf', '_legend.png'),
            name_pdf.replace(name + '.pdf', self.title.replace('/', '|') + '.png'))

    def differential_colorbar(self, dataframe, filename):
        fig_size = (2, 3)
        mpb = plt.pcolormesh(dataframe, cmap='coolwarm')
        fig, ax = plt.subplots(figsize=fig_size)
        plt.colorbar(mpb, ax=ax)
        ax.remove()
        plt.savefig(filename, bbox_inches='tight')

    def differential_expression_sample(
            self, data, samples, ko_column, mmaps2taxa, taxa_column='Taxonomic lineage (GENUS)',
            output_basename=None, log=True):
        """
        Represents in small heatmaps the expression levels of each sample on the
        dataset present in the given pathway map. The values can be transford to
        a log10 scale
        :param data: pandas.DataFrame with data already processed by KEGGPathway
        :param samples: list - column names of the dataset corresponding to expression values
        :param ko_column: str - column with KOs to represent
        :param mmaps2taxa: dict - of taxa to color
        :param taxa_column: str - column with taxonomic classification
        :param output_basename: string - basename of outputs
        :param log: bol - convert the expression values to logarithmic scale?
        """
        if mmaps2taxa is not None:
            data = data[data[taxa_column].isin(mmaps2taxa[self.name.split('ko')[1]])]
        df = data.groupby(ko_column)[samples + [ko_column]].sum()
        df = df[df.any(axis=1)]
        df['Boxes'] = [self.ko_boxes[ko] if ko in self.ko_boxes.keys() else np.nan for ko in df.index]
        df = df[df['Boxes'].notnull()]
        df = expand_by_list_column(df, column='Boxes')
        if len(df) == 0:
            return 1
        df = df.groupby('Boxes')[samples].sum()

        self.pathway_boxes_differential(df, log)

        name = self.name.split(':')[-1]
        name_pdf = f'{output_basename}_{name}.pdf'
        self.to_pdf(name_pdf)

        self.differential_colorbar(df, name_pdf.replace(".pdf", '_legend.png'))

        self.add_legend(
            name_pdf, name_pdf.replace('.pdf', '_legend.png'),
            name_pdf.replace(name + '.pdf', self.title.replace('/', '|') + '.png'))

        return 0

    def add_legend(self, kegg_map_file, legend_file, output):
        """
        Merges the two files - KEGG metabolic map and respective legend - into one file
        :param kegg_map_file: str - filename of PDF kegg metabolic map
        :param legend_file: str - filename of PNG legend
        """
        pdf2png(kegg_map_file)
        imgs = [Image.open(file) for file in [kegg_map_file.replace('.pdf', '.png'), legend_file]]
        imgs[0] = imgs[0].convert(
            'RGB')  # KEGG Maps are converted to RGB by pdftoppm, dunno if converting to RGBA adds any transparency
        imgs[1] = resize_image(imgs[1], ratio=5)
        imgs[1] = add_blank_space(imgs[1], imgs[1].width, imgs[0].height)
        imgs_comb = np.hstack([np.asarray(i) for i in imgs])

        # save that beautiful picture
        imgs_comb = Image.fromarray(imgs_comb)
        imgs_comb.save(output)
        for file in [kegg_map_file, kegg_map_file.replace('.pdf', '.png'), legend_file]:
            os.remove(file)
