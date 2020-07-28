#!/usr/bin/env python

from Bio.KEGG.REST import kegg_link, kegg_list
from kegg_pathway_map import KEGGPathwayMap
from progressbar import ProgressBar
from time import gmtime, strftime
from io import StringIO
from matplotlib import colors, cm
import argparse, pandas as pd, numpy as np, os, pathlib, PIL, re, sys, subprocess
import matplotlib.pyplot as plt

__version__ = "0.0.2"

class KEGGCharter:
    
    def __init__(self, **kwargs):
        self.__dict__ = kwargs
        #self.maps = self.KEGG_metabolic_maps()
        
    def get_arguments(self):    
        parser = argparse.ArgumentParser(description="KEGGCharter - A tool for representing genomic potential and transcriptomic expression into KEGG pathways",
                                         epilog="Input file must be specified.")
        parser.add_argument("-f", "--file", type = str,
                            help="TSV or EXCEL table with information to chart")
        parser.add_argument("-o", "--output", type = str, help = "Output directory",
                            default = 'KEGGCharter_results')
        parser.add_argument("--tsv", action = "store_true", default = False,
                            help="Results will be outputed in TSV format (and not EXCEL).")
        parser.add_argument("-mm", "--metabolic-maps", type = str, 
                            help = "IDs of metabolic maps to output",
                            default = ','.join(self.KEGGCharter_prokaryotic_maps()))
        parser.add_argument("-mgc", "--metagenomic-columns", type = str, 
                            help = "Names of columns with metagenomic quantification")
        parser.add_argument("-mtc", "--metatranscriptomic-columns", type = str, 
                            help = "Names of columns with metatranscriptomics quantification")
        parser.add_argument("-tc", "--taxa-column", type = str, default = 'Taxonomic lineage(GENUS)',
                            help = "Column with the taxa designations to represent with KEGGChart")
        parser.add_argument("-tls", "--taxa-list", type = str, 
                            help = "List of taxa to represent in genomic potential charts (comma separated)")
        parser.add_argument("-not", "--number-of-taxa", type = str, 
                            help = "Number of taxa to represent in genomic potential charts (comma separated)",
                            default = '10')
        parser.add_argument("-koc", "--kos-column", type = str, 
                            help = """"If input file has a column "KO (KEGG Charter)",
                            setting this option will make KEGG Charter use those KOs instead
                            (THIS ARGUMENT OVERRIDES KEGG IDS COLUMNS, USING 
                            KOS DIRECTLY INSTEAD!)""")
        parser.add_argument('-v', '--version', action='version', version='KEGGCharter ' + __version__)
        
        uniprotArgs = parser.add_argument_group('UniProt arguments')
        uniprotArgs.add_argument("-utc", "--uniprot-taxonomic-columns", 
                            action = "store_true", default = False,
                            help="""If columns have UniProt names, KEGGCharter will
                            search for UniProt designations (e.g. Taxonomic lineage(GENUS))""")
        uniprotArgs.add_argument("-tl", "--taxonomic-level", type = str,
                            choices=['SPECIES', 'GENUS', 'FAMILY', 'ORDER', 'CLASS', 'PHYLUM', 'SUPERKINGDOM'],
                            help="The taxonomic level to represent")
        
        special_functions = parser.add_argument_group('Special functions')
        special_functions.add_argument("--show-available-maps", 
                            action = "store_true", default = False,
                            help="""Outputs KEGG maps IDs and descriptions to the
                            console (so you may pick the ones you want!)""")
    
        args = parser.parse_args()
        
        args.output = args.output.rstrip('/')
    
        for directory in [args.output]:
            if not os.path.isdir(directory):
                pathlib.Path(directory).mkdir(parents=True, exist_ok=True)
                print('Created ' + directory)
    
        return args
 
    # Helper functions
    '''
    Input:
        bashCommand: str - command to run
        print_command: boolean - if it should print the command
    Output:
        command will be run... hopefully
    '''
    def timed_message(self, message):
        print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ': ' + message)
        
    def read_input(self, file):
        input_format = 'excel' if file.endswith('.xlsx') else 'tsv'
        if input_format == 'excel':
            return pd.read_excel(file)
        return pd.read_csv(file, sep = '\t')
    
    '''
                       rnaseq_a0.17reads_normalized  rnaseq_a1reads_normalized  ...  rnaseq_c3reads_normalized     Boxes
    KO (KEGG Charter)                                                           ...
    K02108                                17.396381                 118.667524  ...                 318.212039   [2, 64]
    K02109                                37.862711                 158.893803  ...                 365.895189   [1, 64]
    >>>
        rnaseq_a0.17reads_normalized  rnaseq_a1reads_normalized  rnaseq_a3reads_normalized  ...  rnaseq_c1reads_normalized  rnaseq_c3reads_normalized  Boxes
    0                      17.396381                 118.667524                 316.978069  ...                 129.690980                 318.212039      2
    1                      17.396381                 118.667524                 316.978069  ...                 129.690980                 318.212039     64
    2                      37.862711                 158.893803                 387.740757  ...                 164.878456                 365.895189      1
    3                      37.862711                 158.893803                 387.740757  ...                 164.878456                 365.895189     64
    '''
    def expand_by_list_column(self, df, column = 'Pathway'):
        if len(df) == 0:
            return pd.DataFrame()
        lens = [len(item) for item in df[column]]
        dictionary = dict()
        for column in df.columns:
            dictionary[column] = np.repeat(df[column].values,lens)
        dictionary[column] = np.concatenate(df[column].values)
        return pd.DataFrame(dictionary) 
        
    # Functions that deal with taxonomies
    def most_abundant_taxa(self, data, samples, number_of_taxa = 10, 
                           level_of_taxa = 'GENUS'):
        '''
        Calculates top genus from samples
        :param samples: list of mg samples to consider for quantification of genus abundance
        :param n: number of top genus to return
        :return: list of top genus
        '''
        data = data.groupby("Taxonomic lineage ({})".format(level_of_taxa))[samples].sum(
            ).sum(axis = 1).sort_values(ascending=False)
        if number_of_taxa > len(data.index.tolist()):
            number_of_taxa = len(data.index.tolist())
        return data.index.tolist()[:number_of_taxa]
    
    def taxa_colors(self, hex_values = None, ncolor = 1):
        '''
        Creates list of hex colors to be used, using matplotlib or using custom colors
        :param hex_values: list of hex colors
        :param ncolor: int indicating the ammount of hex colors that should be created
        :return: returns list with hex color codes
        '''
        if not hex_values:                                                      # if no colors are given creates a list of hex colors with ncolor from matplotlib discrete colormaps
            color_scheme = (cm.get_cmap('Pastel2', 8) if ncolor <= 8
                            else cm.get_cmap("Set3", 12) if ncolor <= 12
                            else cm.get_cmap("rainbow", ncolor))                # if ncolor > 12 a continuous colormap is used instead
            return [colors.to_hex(color_scheme(i)) for i in range(ncolor)]
        else:
            for hex_value in hex_values:
                if not re.search(r'^#(?:[0-9a-fA-F]{3}){1,2}$', hex_value):
                    sys.exit(Exception("Colors aren't valid hex codes"))
            return hex_values                                                    # has validated hex values and returns the original list
                
    
    def get_organisms(self, file):                                              # TODO - is still not in use, but could be to automatically determine most abundant organisms
        '''
        Data obtained from www.genome.jp/kegg/catalog/org_list.html was organized
        to retrieve the identifiers of the organisms available at KEGG
        :param file: string - filename containing the information
        '''
        return pd.read_csv(file, sep = '\t', index_col = 0, header = None)
            
    # Conversion functions
    def keggid2ko(self, kegg_ids, step = 150):
        '''
        Converts KEGG_ID genes to Ortholog KO ID from KEGG
        :param KEGG_ID: (list) - KEGG ID genes
        :param step: (int) - will convert "step" KEGG IDs at a time
        :return: (list) - (list,list) - KEGG ID genes converted and ko IDs
        '''
        result = list()
        pbar = ProgressBar()
        for i in pbar(range(0, len(kegg_ids) - step, step)):
            try:
                result += kegg_link("ko", kegg_ids[i:i+step]).read().split("\n")[:-1]
            except:
                print('KEGG ID to KO broke at index ' + str(i))
                result = [[part[0] + ';', part[1].strip('ko:')] for part in
                   [relation.split('\t') for relation in result]]
                return pd.DataFrame(result, columns = ['Cross-reference (KEGG)', 'KO (KEGG Charter)'])
        result += kegg_link("ko", kegg_ids[len(kegg_ids) - step:]).read().split("\n")[:-1]
        result = [[part[0] + ';', part[1].strip('ko:')] for part in
                   [relation.split('\t') for relation in result]]
        return pd.DataFrame(result, columns = ['Cross-reference (KEGG)', 'KO (KEGG Charter)'])
    
    def ko2ec(self, kos, step = 150):
        '''
        Converts KOs to EC numbers
        :param kos: list of kegg orthologs
        :return: dic associating ortholog kegg id with list
        of assotiated EC numbers
        '''
        result = list()
        pbar = ProgressBar()
        for i in pbar(range(0, len(kos), step)):
            try:
                result += kegg_link("enzyme", kos[i:i+step]).read().split("\n")[:-1]
            except:
                print('KO to EC number broke at index ' + str(i))
                result = [relation.split('\t') for relation in result]
                return list(map(list, zip(*result)))
        result += kegg_link("enzyme", kos[len(kos) - step:]).read().split("\n")[:-1]
        result = [[part[0].strip('ko:'),part[1].upper()] for part in
                   [relation.split('\t') for relation in result]]
        return pd.DataFrame(result, columns = ['KO (KEGG Charter)', 'EC number (KEGG Charter)'])
    
    # Get metabolic maps from KEGG Pathway
    def KEGGCharter_prokaryotic_maps(self, file = sys.path[0] + '/KEGGCharter_prokaryotic_maps.txt'):
        return open(file).read().split('\n')
    
    def KEGG_metabolic_maps(self):
        '''
        Creates a dic with all specific kegg pathway maps and their description
        :return: pandas.DataFrame with Map ID as index and maps names as
        sole column
        '''
        maps = pd.read_csv(StringIO(kegg_list("pathway").read()), sep='\t',
                           names = ['Map ID', 'Description'])
        maps['Map ID'] = [ide.split('map')[-1] for ide in maps['Map ID']]
        return maps
    
    # Methods related with image manipulation to get those maps with those legends
    def add_blank_space(self, image_pil, width, height, image_mode = 'RGB'):
        '''
        Resizes an image with white background, keeping image size ratio
        :param image_pil: PIL.Image - image to be resized
        :param width: int - width of final image
        :param height: int - heigth of final image
        :param image_mode: str - image mode of image (RGBA, RGB, ...)
        '''
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
    
    def resize_image(self, image_pil, ratio = None, width = None, height = None):
        '''
        Resizes an image with alteration to image size ratio
        :param ratio: int - ratio of resize - how bigger or smaller will the output be?
        :param image_pil: PIL.Image - image to be resized
        :param width: int - width of final image
        :param height: int - heigth of final image
        '''
        if ratio:
            return image_pil.resize((image_pil.width * ratio, 
                                     image_pil.height * ratio), PIL.Image.ANTIALIAS)
        elif width and height:
            return image_pil.resize((width, height), PIL.Image.ANTIALIAS)
        else:
            return None
        
    def pdf2png(self, pdf_filename):
        '''
        Converts a pdf file to a png file, RGB format - name changes, .pdf to .png
        :param pdf_filename: str - filename of PDF file
        '''
        bashCommand = 'pdftoppm {} {} -png'.format(pdf_filename, pdf_filename.split('.pdf')[0])
        subprocess.run(bashCommand.split())
        os.rename(pdf_filename.replace('.pdf', '-1.png'), pdf_filename.replace('.pdf', '.png'))
    
    def add_legend(self, kegg_map_file, legend_file, output):
        '''
        Merges the two files - KEGG metabolic map and respective legend - into
        one file file
        :param kegg_map_file: str - filename of PDF kegg metabolic map
        :param legend_file: str - filename of PNG legend
        '''
        self.pdf2png(kegg_map_file)
        imgs = [PIL.Image.open(file) for file in 
                [kegg_map_file.replace('.pdf', '.png'), legend_file]]
        imgs[0] = imgs[0].convert('RGB')                                       # KEGG Maps are converted to RGB by pdftoppm, dunno if converting to RGBA adds any transparency
        imgs[1] = self.resize_image(imgs[1], ratio = 5)
        imgs[1] = self.add_blank_space(imgs[1], imgs[1].width, imgs[0].height)
        imgs_comb = np.hstack([np.asarray(i) for i in imgs])
        
        # save that beautiful picture
        imgs_comb = PIL.Image.fromarray(imgs_comb)
        imgs_comb.save(output)
        
    # Main functions for genomic potential charting
    def create_potential_legend(self, colors, labels, filename, resize_factor = 10):
        '''
        Draws the color to taxa labels of genomic potential representations
        :param colors: list - list of colors of the different taxa
        :param labels: list - list of taxa corresponding to the colors
        :param filename: string - filename to output
        :param size: int - how big you want your legend?
        '''
        f = lambda m, c: plt.plot([], [], marker = m, color = c, ls = "none")[0]
        handles = [f("s", color) for color in colors]
        legend = plt.legend(handles, labels, loc=3, framealpha=1, frameon = True)
        fig = legend.figure
        fig.canvas.draw()
        # The bbox manipulation removes the axis
        bbox = legend.get_window_extent()
        bbox = bbox.from_extents(*(bbox.extents + np.array([-2,-2,2,2])))
        bbox = bbox.transformed(fig.dpi_scale_trans.inverted())
        
        fig.savefig(filename, dpi="figure", bbox_inches=bbox)



    def genomic_potential_taxa(self, kegg_pathway_map, data, samples, dic_colors, 
                               taxa = None, taxa_column = 'Taxonomic lineage (GENUS)', 
                               output_basename = None, maxshared = 10):
        '''
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
        '''        
        # for every taxon, check all boxes it is in, and save that info to box2taxon
        box2taxon = dict()
        for taxon in dic_colors.keys():
            df = data[data[taxa_column] == taxon][samples + ['KO (KEGG Charter)']]
            df = df[df.any(axis=1)]
            for ortholog in df['KO (KEGG Charter)']:
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
        df = data[samples + ['KO (KEGG Charter)']]
        df = df[df.any(axis=1)]
        grey_boxes = list()
        for ortholog in df['KO (KEGG Charter)']:
            if ortholog in kegg_pathway_map.ko_boxes.keys():
                for box in kegg_pathway_map.ko_boxes[ortholog]:
                    if box not in box2taxon.keys() and box not in grey_boxes:
                        grey_boxes.append(box)
        kegg_pathway_map.grey_boxes(grey_boxes)                                 # TODO - boxes should have EC number as name
        
        name = kegg_pathway_map.pathway.name.split(':')[-1]
        name_pdf = '{}_{}.pdf'.format(output_basename, name)
        kegg_pathway_map.to_pdf(name_pdf)
        
        # a new taxon: "Other taxa"
        if len(grey_boxes) > 0:
            dic_colors["Other taxa"] = "#7c7272"

        # TODO - legend should be ajusted for the maps - otherwise, makes no sense to have one legend for each map - they all become the same, except for "Other taxa"
        self.create_potential_legend(dic_colors.values(), dic_colors.keys(), 
                                     name_pdf.replace('.pdf','_legend.png'))
        
        self.add_legend(name_pdf, name_pdf.replace('.pdf','_legend.png'), 
            name_pdf.replace(name + '.pdf', kegg_pathway_map.pathway.title.replace('/','|') + '.png'))
    
    # Main functions for differential expression charting
    def differential_colorbar(self, dataframe, filename):
        FIGSIZE = (2,3)
        mpb = plt.pcolormesh(dataframe,cmap='coolwarm')
        fig,ax = plt.subplots(figsize=FIGSIZE)
        plt.colorbar(mpb,ax=ax)
        ax.remove()
        plt.savefig(filename,bbox_inches='tight')
        
    def differential_expression_sample(self, kegg_pathway_map, data, samples, 
                                       output_basename = None, log = True):
        '''
        Represents in small heatmaps the expression levels of each sample on the
        dataset present in the given pathway map. The values can be transford to
        a log10 scale
        :param data: pandas.DataFrame with data already processed by KEGGPathway
        :param samples: list - column names of the dataset corresponding to
        expression values
        :param output_folder: string - name of folder to store pdfs
        :param log: bol - convert the expression values to logarithmic scale?
        '''
        df = data.groupby('KO (KEGG Charter)')[samples + ['KO (KEGG Charter)']].sum()
        df = df[df.any(axis=1)]
        df['Boxes'] = [kegg_pathway_map.ko_boxes[ko] if ko in 
            kegg_pathway_map.ko_boxes.keys() else np.nan for ko in df.index]
        df = df[df['Boxes'].notnull()]
        df = self.expand_by_list_column(df, column = 'Boxes')
        if len(df) == 0:
            return 1
        df = df.groupby('Boxes')[samples].sum()

        kegg_pathway_map.pathway_boxes_differential(df, log)
        
        name = kegg_pathway_map.pathway.name.split(':')[-1]
        name_pdf = '{}_{}.pdf'.format(output_basename, name)
        kegg_pathway_map.to_pdf(name_pdf)
            
        self.differential_colorbar(df, name_pdf.replace(".pdf",'_legend.png'))
        
        self.add_legend(name_pdf, name_pdf.replace('.pdf','_legend.png'), 
            name_pdf.replace(name + '.pdf', kegg_pathway_map.pathway.title.replace('/','|') + '.png'))
        
        return 0
    
    def write_results(self, data, output_dir, output_type = 'excel'):
        if output_type == 'excel':
            data.to_excel('{}/KEGGCharter_results.xlsx'.format(output_dir))
        else:
            data.to_csv('{}/KEGGCharter_results.tsv'.format(output_dir), sep = '\t')
    
    def main(self):
        # Read input
        args = self.get_arguments()
        self.timed_message('Arguments valid.')
        
        if args.show_available_maps:
            exit(self.KEGG_metabolic_maps().to_string())
        else:
            if not hasattr(args, 'file'):
                sys.exit("Must specify input file with --file.")
        data = self.read_input(args.file)
        self.timed_message('Data successfully read.')
        
        # KEGG ID to KO -> if KO column is not set, it will get them with the KEGG API
        if not args.kos_column:
            kegg_ids = data[data['Cross-reference (KEGG)'].notnull()]['Cross-reference (KEGG)']
            kegg_ids = [ide.split(';')[0] for ide in kegg_ids]                      # sometimes Kegg IDs come as mfc:BRM9_0145;mfi:DSM1535_1468; from UniProt IDs mapping. Should be no problem since both IDs should always be same function
            self.timed_message('Converting {:d} KEGG IDs to KOs through the KEGG API.'.format(len(kegg_ids)))
            kos = self.keggid2ko(kegg_ids)
            data = pd.merge(data, kos, on = 'Cross-reference (KEGG)', how = 'left')
            kos_column = 'KO (KEGG Charter)'
        else:
            kos_column = args.kos_column

        # KO to EC number
        kos = data[data[kos_column].notnull()][kos_column].tolist()
        self.timed_message('Retrieving EC numbers from {} KOs.'.format(len(kos)))
        ecs = self.ko2ec(kos)
        data = pd.merge(data, ecs, on = 'KO (KEGG Charter)', how = 'left')
        self.timed_message('Merging UniProt and KEGG information on EC numbers')
        
        self.timed_message('Results saved to {}/KEGGCharter_results.{}'.format(
            args.output, 'tsv' if args.tsv else 'xlsx'))
        self.write_results(data, args.output, 
                           output_type = ('tsv' if args.tsv else 'excel'))
        
        # Begin dat chart magic
        metabolic_maps = args.metabolic_maps.split(',')
        self.timed_message('Creating KEGG Pathway representations for {} metabolic pathways.'.format(
                str(len(metabolic_maps))))
        
        not_loaded = list(); differential_no_kos = list()
        i = 1
        
        # Set colours for taxa if MG data is present
        if hasattr(args, 'metagenomic_columns'):
            args.metagenomic_columns = args.metagenomic_columns.split(',')
            if args.taxa_list is None:
                taxa = self.most_abundant_taxa(data, args.metagenomic_columns, 
                                   number_of_taxa = int(args.number_of_taxa),
                                   level_of_taxa = args.taxonomic_level.upper())
            else:
                taxa = args.taxa_list.split(',')
                
            colors = self.taxa_colors(ncolor = len(taxa))
            dic_colors = {taxa[i] : colors[i] for i in range(len(taxa))}
            
        if hasattr(args, 'metatranscriptomic_columns'):
            args.metatranscriptomic_columns = args.metatranscriptomic_columns.split(',')
        
        # For each metabolic map, will chart genomic potential and differential expression if MG and MT data are available, respectively
        for metabolic_map in metabolic_maps:
            kegg_pathway_map = KEGGPathwayMap(data, metabolic_map)              # if getting 404 here, is likely because map doesn't exist # TODO - put something here to inform users to report that problem
            print('[{}/{}] Creating representation for pathway: {}'.format(
                    str(i), str(len(metabolic_maps)), 
                    kegg_pathway_map.pathway.title))
            if args.metagenomic_columns:
                kegg_pathway_map = KEGGPathwayMap(data, metabolic_map)
                if args.uniprot_taxonomic_columns:
                    args.taxa_column = 'Taxonomic lineage ({})'.format(args.taxonomic_level.upper())
                self.genomic_potential_taxa(kegg_pathway_map, data,
                        args.metagenomic_columns, dic_colors, 
                        taxa_column = args.taxa_column,
                        output_basename = args.output + '/potential')           # TODO - log should be True, fix the other TODO
            if args.metatranscriptomic_columns:
                kegg_pathway_map = KEGGPathwayMap(data, metabolic_map)
                saul_good = self.differential_expression_sample(kegg_pathway_map, data, 
                        args.metatranscriptomic_columns, 
                        output_basename = args.output + '/differential',
                        log = False)
                if saul_good == 1:
                    differential_no_kos.append(metabolic_map)
            plt.close()
            i += 1
            
        not_loaded_filename = args.output + '/not_loaded.txt'
        no_kos_filename = args.output + '/no_kos.txt'
        print('{} maps could not be loaded. You can see which ones at {}'.format(
            len(not_loaded), not_loaded_filename))
        open(not_loaded_filename, 'w').write('\n'.join(not_loaded))
        print('{} maps had no KOs attributed to them. You can see which ones at {}'.format(
            len(differential_no_kos), no_kos_filename))
        open(no_kos_filename, 'w').write('\n'.join(differential_no_kos))
        
if __name__ == '__main__':
    KEGGCharter().main()