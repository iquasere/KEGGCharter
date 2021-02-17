#!/usr/bin/env python

from Bio.KEGG.KGML import KGML_parser, KGML_pathway
from Bio.Graphics.KGML_vis import KGMLCanvas
from matplotlib import colors, cm
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


class KEGGPathwayMap:
    """
    This class retrieves and manipulates KEGG metabolic maps from KEGG Pathway
    """

    def __init__(self, pathway, ec_list):
        """
        Initialize object
        :param data: pd.DataFrame - data from MOSCA analysis
        :param kgml_file: (str) - KGML filename
        """
        self.pathway = pathway
        self.set_pathway(ec_list)

    def set_pathway(self, ec_list):
        """
        Set pathway with Kegg Pathway ID
        """
        self.ko_boxes = dict()
        for i in range(len(self.pathway.orthologs)):
            set_bgcolor(self.pathway.orthologs[i], "#ffffff")  # set all boxes to white
            # self.set_fgcolor(self.pathway.orthologs[i], "#ffffff")             # This might be helpful in the future, if an additional layer of marking is needed
            orthologs_in_box = [ide[3:] for ide in self.pathway.orthologs[
                i].name.split()]  # 'ko:K16157 ko:K16158 ko:K16159' -> ['K16157', 'K16158', 'K16159']
            for ortholog in orthologs_in_box:
                if ortholog not in self.ko_boxes.keys():
                    self.ko_boxes[ortholog] = list()
                self.ko_boxes[ortholog].append(i)  # {'K16157':[0,13,432], 'K16158':[4,13,545]}

            # Set name as EC number
            ecs = ec_list[i].split(',')
            if len(ecs) > 0:
                self.pathway.orthologs[i].graphics[0].name = max(set(ecs), key=ecs.count).upper()
            else:
                self.pathway.orthologs[i].graphics[0].name = orthologs_in_box[0]

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
        drawn = False
        while not drawn:
            try:
                KGMLCanvas(self.pathway,
                           import_imagemap=imagemap,
                           label_orthologs=orthologs,
                           label_compounds=compounds,
                           label_maps=maps,
                           label_reaction_entries=reactions).draw(filename)
                drawn = True
            except:
                print('Draw failed! Waiting 10 secs...')
                time.sleep(10)

    def pathway_box_list(self, taxa_in_box, dic_colors, maxshared=10):
        """
        Represents items in the pathway map
        :param taxa_in_box: dict - {box : list of taxa in box}
        :param taxa: list - of taxa to be represented and given a specific color
        :param maxshared: int - maximum number of taxa sharing one box
        :param color: list of costum colors to be used to color the elements
        """
        for box in taxa_in_box.keys():
            nrboxes = len(taxa_in_box[box])
            if nrboxes > maxshared:
                nrboxes = maxshared

            paired = True if nrboxes % 2 == 0 else False
            for i in range(nrboxes):
                newrecord = create_box_heatmap(
                    self.pathway.orthologs[box], nrboxes, i * 2 - (nrboxes - 1) if paired else i - int(nrboxes / 2),    # if nrboxes = 8, i * 2 - (nrboxes - 1) = -7,-5,-3,-1,1,3,5,7 # if nrboxes = 9, i - int(nrboxes / 2) = -4,-3,-2,-1,0,1,2,3,4
                    paired=paired)
                if newrecord != 1:
                    newrecord.bgcolor = dic_colors[taxa_in_box[box][i]]
                    self.pathway.orthologs[box].graphics.append(newrecord)
            if self.pathway.orthologs[box].graphics[0].width is not None:  # TODO - should check more deeply why sometimes width is None
                create_tile_box(self.pathway.orthologs[box])

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
                    self.pathway.orthologs[box], nrboxes, i * 2 - (nrboxes - 1) if paired else i - int(nrboxes / 2),    # if nrboxes = 8, i * 2 - (nrboxes - 1) = -7,-5,-3,-1,1,3,5,7; if nrboxes = 9, i - int(nrboxes / 2) = -4,-3,-2,-1,0,1,2,3,4
                    paired=paired)
                if newrecord != 1:  # TODO - assess why sometimes get 1
                    newrecord.bgcolor = colors[i]
                    self.pathway.orthologs[box].graphics.append(newrecord)
            if self.pathway.orthologs[box].graphics[0].width is not None:  # TODO - should check more deeply why sometimes width is None
                create_tile_box(self.pathway.orthologs[box])

    def grey_boxes(self, box_list):
        for i in box_list:
            set_bgcolor(self.pathway.orthologs[i], "#7c7272")
            set_fgcolor(self.pathway.orthologs[i], "#7c7272")
