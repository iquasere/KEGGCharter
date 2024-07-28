import json
import csv
import xml.etree.ElementTree as ET
import os


def load_json_to_dict(filepath: str) -> dict[str, any]:
    """
    Loads a JSON file from a specified path and returns the contents as a dictionary.

    Parameters
    ----------
    filepath : str
        The path to the JSON file.

    Returns
    -------
    dict
        The contents of the JSON file as a dictionary.
    """
    with open(filepath, 'r') as file:
        data = json.load(file)
    return data


def load_kgml(caminho_kgml : str) -> None:
    """
    Loads a KGML file and returns the root of the XML document.

    Parameters
    ----------
    caminho_kgml : str
        The path to the KGML file.

    Returns
    -------
        The root of the XML document.
    """

    tree = ET.parse(caminho_kgml)
    root = tree.getroot()
    return root


def extract_pathway_title(caminho_kgml: str) -> str:
    """
    Extracts the pathway title from a KGML file.

    Parameters
    ----------
    caminho_kgml : str
        The path to the KGML file.

    Returns
    -------
    str
        The title of the pathway.
    """

    tree = ET.parse(caminho_kgml)
    root = tree.getroot()
    
    titulo_pathway = root.get('title', 'Title not found')
    return titulo_pathway


def extract_pathway_number(caminho_kgml: str) -> str:
    """
    Extracts the pathway number from a KGML file.

    Parameters
    ----------
    caminho_kgml : str
        The path to the KGML file.

    Returns
    -------
    str
        The pathway number.
    """

    tree = ET.parse(caminho_kgml)
    root = tree.getroot()
    
    numero_pathway = root.get('number', 'Number not found')
    return numero_pathway


def extract_graphical_elements_with_id(root, tipo: str = 'rectangle') -> list[dict[str, str]]:
    """
    Extracts graphic elements of a specific type from the root XML document and includes their IDs, 
    coordinates and a link associated with each ID.

    Parameters
    ----------
    root : xml.etree.ElementTree.Element
        The root of the XML document.

    tipo : str
        Type of graphic element to be extracted (the default is ‘rectangle’).

    Returns
    -------
    list[dict[str, str]]
        List of dictionaries containing data on graphic elements, including their IDs,
        the coordinates (x, y, width, height) and an associated link.
    """
    elementos_id = []
    
    for entry in root.findall(".//entry"):
        graphics = entry.find(".//graphics")
        
        if graphics is not None and graphics.get('type') == tipo:
            link_suffix = entry.get('link')
            if link_suffix and not link_suffix.startswith('http'):
                link = f"https://www.kegg.jp/dbget-bin/www_bget?{link_suffix}"
            else:
                link = link_suffix  
            
            dados_elemento = {
                'id': entry.get('id'),  
                'x': graphics.get('x'),
                'y': graphics.get('y'),
                'width': graphics.get('width'),
                'height': graphics.get('height'),
                'link': link
            }
            elementos_id.append(dados_elemento)
    
    return elementos_id


def enrich_graphic_elements(elementos_graficos_com_id: list[dict[str, str]], caminho_tsv: str, caminho_json_taxon: str, caminho_json_kos: str) -> list[dict[str, str]]:
    """
    It enriches the graphic elements with data from a TSV file and two JSON files.

    Parameters
    ----------
    elementos_graficos : list[dict[str, str]]
        List of dictionaries with data on graphic elements.

    caminho_tsv : str
        Path to the TSV file containing ID-to-name mappings (EC or KO).

    caminho_json_taxon : str
        Path to the JSON file containing ID-to-taxon mappings.

    caminho_json_kos : str
        Path to the JSON file containing ID mappings for KOs.

    Returns
    -------
    list[dict[str, str]]
        List of dictionaries with data on graphic elements enriched with additional information.
    """
    id_to_name = {}
    with open(caminho_tsv, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        next(reader)  
        for linha in reader:
            if linha:
                id_to_name[linha[0]] = linha[1]

    id_to_taxon = load_json_to_dict(caminho_json_taxon)
    id_to_kos = load_json_to_dict(caminho_json_kos)

    for elemento in elementos_graficos_com_id:
        elemento_id = elemento['id']
        elemento['EC_KO'] = id_to_name.get(elemento_id, ["Unknown"])
        elemento['taxon'] = id_to_taxon.get(elemento_id, ["Unknown"])
        elemento['kos'] = id_to_kos.get(elemento_id, ["Unknown"])

    return elementos_graficos_com_id

#######################################não faz nada#############


def create_coloured_boxes(id : str, colors_dict : dict[str, list[str]]) -> str:
    """
    Generates an HTML string representing several horizontally aligned sub-boxes,
    each coloured as specified in the dictionary for a specific ID.
    The sub-boxes are created with a black border for better visualisation.

    Parameters
    ----------
        id : str 
            The identifier that points to the list of colours in the dictionary.

        colors_dict : dict[str, list[str]] 
            A dictionary mapping identifiers to lists of colour strings.

    Returns
    -------
        str 
            An HTML string representing a row of coloured sub-boxes, each with a black border.
    """
    colors = colors_dict.get(str(id), ["#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF"])  
    boxes_html = ''.join([
        f'<div style="background-color: {color}; width: 25%; height: 50px; display: inline-block; border: 1px solid black;"></div>'
        for color in colors
    ])
    return f'<div style="width: 100%; height: 50px; display: flex;">{boxes_html}</div>'


def load_tsv_by_id(diretoria: str, id: int) -> dict:
    """
    Searches for a TSV file in the specified directory that matches the ID provided and converts it to a dictionary.
    
    Parameters
    ----------
        diretoria : str
            The path of the directory where the TSV files are stored.
        id : int
            The ID that identifies the specific TSV file.

    Returns
    -------
        dict
            Dictionary with the data structure extracted from the TSV file.

    Raises
    ------
        FileNotFoundError
            Se o arquivo com o ID especificado não for encontrado.
    """
    nome_arquivo = f"ko00680_{id}_info.json"  
    caminho_completo = os.path.join(diretoria, nome_arquivo)

    if not os.path.exists(caminho_completo):
        return None

    resultado = {}
    with open(caminho_completo, 'r', newline='', encoding='utf-8') as arquivo:
        reader = csv.DictReader(arquivo, delimiter='\t')
        for row in reader:
            taxa = row.pop('Taxonomic lineage (GENUS)')
            resultado[taxa] = {key: float(val) for key, val in row.items()}

    return resultado


def create_sample_boxes(taxa_expressions: dict, taxa_colors: dict) -> list:
    """
    Generates a list of HTML strings, each representing the differential expressions of a sample
    in the form of a coloured vertical box.

    Parameters
    ----------
        taxa_expressions : dict
            Dictionary with taxonomies and their expressions for each sample.
        taxa_colors : dict
            Dictionary with the corresponding colour for each taxonomy.

    Returns
    -------
        list
            List of HTML strings, each representing a box for a sample.
    """
    if not isinstance(taxa_colors, dict):
        raise ValueError("taxa_colors should be a colour dictionary")

    html = ''
    column = 1
    for amostra in ['MP1', 'MP2', 'MP3', 'MP4']:
        html_content = f'<div class="grid-item grid-item-center" style="grid-row-start: 4; grid-row-end: 4; grid-column-start: {column}; grid-column-end: {column + 1};">'
        html_content += f'<h3>{amostra}</h3>'
        total_expression = sum(taxa_expressions[taxa][amostra] for taxa in taxa_expressions)
        if total_expression == 0:
            html_content += '<div style="background-color: #FFFFFF; width: 50px; height: 300px; margin-bottom: 2px;">0</div>'
        else:
            for taxa, values in taxa_expressions.items():
                height = f"{(values[amostra] / total_expression) * 300}px"
                if values[amostra] / total_expression > 0:
                    color = taxa_colors.get(taxa, '#FFFFFF')
                    html_content += f'<div style="background-color: {color}; width: 50px; height: {height}; margin-bottom: 2px;">{round(values[amostra], 3)}</div>'
        html_content += '</div>'
        html += html_content
        column +=1

    return html


def create_detail_page(id : str, link : str, ec_ko : str, KO : list[str], Taxon : list[str], colors_dict : dict[str, list[str]], taxa_expressions : dict, taxa_colors : dict):
    """
    Creates a detailed HTML page for a specific ID based on the data provided and the colours defined in the dictionary.
    
    Parameters
    ----------
        id : str 
            The element identifier, used to create specific coloured boxes and to name the HTML file.
        link : str
            The URL to which the button on the HTML page should redirect.
        ec_ko : str 
            The EC or KO number that will be used as the page title.
        KO : list[str] 
            List of K numbers that will be displayed on the page.
        Taxon : list[str] 
            List of associated taxon names that will be displayed on the page.
        colors_dict : dict[str, list[str]] 
            Dictionary mapping identifiers to colour lists, used to colour boxes on the page.
        taxa_expressions : dict
            Dictionary with taxonomies and their expressions for each sample.
        taxa_colors : dict
            Dictionary with the corresponding colour for each taxonomy.

    Effects
    -------
        Creates an HTML file in the current directory with a name based on the pathway and ID provided.
    """
    numero_pathway = extract_pathway_number(caminho_kgml)
    titulo_pathway = extract_pathway_title(caminho_kgml)
    color_boxes = create_coloured_boxes(id, colors_dict)

    sample_boxes = create_sample_boxes(taxa_expressions, taxa_colors)  if taxa_expressions is not None else ''

    imagem_url_differential = "./KEGGCharter/first_time_running_KC/maps/differential_ko00680_legend.png"
    imagem_url_potential = "./KEGGCharter/first_time_running_KC/maps/potential_ko00680_legend.png"

    detalhe_html = f"""
    <!DOCTYPE html>
    <html lang='en'>
    <head>
        <meta charset='UTF-8'>
        <title>EC_Number {ec_ko}</title>
        <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 0;
            padding: 50px;
            background-color: #f4f4f4;
        }}
        .grid-container {{
            display: grid;
            grid-template-columns: repeat(5, 1fr);
            gap: 0px 10px;
        }}
        .grid-item-center {{
            display: flex;
            flex-direction: column; 
            align-items: center;    
            justify-content: center;  
        }}
        #botao {{
            grid-row-start: 1;
            grid-row-end: 2;
            grid-column-start: 5;
            grid-column-end: 5;
        }}
        #titulo {{
            grid-row-start: 1;
            grid-row-end: 2;
            grid-column-start: 1;
            grid-column-end: 5;
        }}
        #img_differential {{
            grid-row-start: 3;
            grid-row-end: 5;
            grid-column-start: 5;
            grid-column-end: 5;
        }}
        #box_differential {{
            grid-row-start: 3;
            grid-row-end: 3;
            grid-column-start: 1;
            grid-column-end: 5;
        }}
        #sample_boxes {{
            grid-row-start: 4;
            grid-row-end: 5;
            grid-column-start: 1;
            grid-column-end: 5;
        }}
        #expression_by_taxonomy {{
            grid-row-start: 4;
            grid-row-end: 5;
            grid-column-start: 1;
            grid-column-end: 2;
            font-weight: bold;
            writing-mode: vertical-rl;
            transform: rotate(180deg); 
        }}
        #total_expression {{
            grid-row-start: 2;
            grid-row-end: 3;
            grid-column-start: 1;
            grid-column-end: 2;
            font-weight: bold;
        }}
        </style>
    </head>
    <body>
        <div class="grid-container">
            <div class="grid-item" id = "botao" style="text-align: center;">
                <button onclick="location.href='{link}';">Link KEGG</button>
            </div>
            <div class="grid-item" id = "titulo">
                <h1 style="text-align: center;">EC Number - {ec_ko}</h1>
                <h3>MAP {numero_pathway} ({titulo_pathway})</h3>
                <p>Box ID - {id} <br />
                K numbers - {', '.join(KO)} <br />
                Taxonomy - {', '.join(Taxon)}</p>
            </div>
            <div class="grid-item grid-item-center" id = "img_differential" style="text-align: center;">
                <img src="{imagem_url_differential}" alt="Color Legend">
                <img src="{imagem_url_potential}" alt="Color Legend">
            </div>
            <div class="grid-item" id = "box_differential">
                {color_boxes} 
            </div>
                {sample_boxes}
            <div class="grid-item" id = "expression_by_taxonomy" style="text-align: center;">
                Expression by Taxonomy:
            </div>
            <div class="grid-item" id = "total_expression">
                Total Expression:
            </div>
        </div>
    </body>
    </html>
    """
    
    with open(f"Map{numero_pathway}_Box{id}.html", 'w') as file:
        file.write(detalhe_html)


def create_image_maps(taxa_colors : dict, coordenadas: list[dict[str, str]],
                     imagem: str = "./KEGGCharter/original_kegg_map.png", arquivo_saida: str = f"image_maps.html",
                     titulo_imagem: str = "Methane metabolism Map",
                     colors_dict: str = './KEGGCharter/info/ko00680_boxes2quant_colors.json') -> None:
    """
    Creates an HTML file with image maps based on specified coordinates.

    Parameters
    ----------
    taxa_colors : dict
        Dictionary with the corresponding colour for each taxonomy.
            
    coordenadas : list[dict[str, str]]
        Creates an HTML file with image maps based on specified coordinates.

    imagem : str
        Path to the background image.

    arquivo_saida : str
        Path of the output file where the HTML will be saved.

    titulo_imagem : str
        Title of the image to be included in the HTML.

    colors_dict : str

    Returns
    -------
    None
    """

    def coords_to_area(coords: dict[str, str]) -> str:
        """
        Converts coordinates of a graphic element into an area string for use in HTML image maps.

        Parameters
        ----------
        coords : dict[str, str]
            Dictionary with the element's 'x', 'y', 'width' and 'height' coordinates.

        Returns
        -------
        str
            String formatted for use in the 'coords' attribute of an HTML <area> tag.
        """

        x1 = int(coords['x']) * 2.1
        y1 = int(coords['y']) * 2.1
        x2 = x1 + int(coords['width']) * 2.2
        y2 = y1 + int(coords['height']) * 2.2
        return f"{x1},{y1},{x2},{y2}"

    numero_pathway = extract_pathway_number(caminho_kgml)
    html = f"<!DOCTYPE html>\n<html lang='en'>\n<head>\n    <meta charset='UTF-8'>\n    <title>{titulo_imagem}</title>\n</head>\n<body>\n    <h2>{titulo_imagem}</h2>\n    <img src='{imagem}' usemap='#potentialMap' alt='{titulo_imagem}'>\n    <map name='potentialMap'>\n"

    for coord in coordenadas:
        id = coord['id']
        link = coord['link']
        ec_ko = coord['EC_KO']
        KO = coord['kos']
        Taxon = coord['taxon']
        taxa_expressions = load_tsv_by_id(diretoria, id)
        create_detail_page(id, link, ec_ko, KO, Taxon, colors_dict, taxa_expressions, taxa_colors)
        area_str = coords_to_area(coord)
        html += f'        <area shape="rect" coords="{area_str}" href="Map{numero_pathway}_Box{id}.html" target="_blank" alt="Element ID {id}">\n'

    html += "</map>\n</body>\n</html>"

    with open(arquivo_saida, 'w') as file:
        file.write(html)

    print(f"HTML with image maps successfully generated and saved as {arquivo_saida}!")


if __name__ == '__main__':
    caminho_kgml = './resources_directory/kc_kgmls/ko00680.xml'
    caminho_tsv = './KEGGCharter/info/ko00680_box2name.tsv'
    caminho_json_taxon = './KEGGCharter/info/ko00680_boxes2taxon.json'
    caminho_json_kos = './KEGGCharter/info/ko00680_box2kos.json'
    caminho_json_colors = './KEGGCharter/info/ko00680_boxes2quant_colors.json'
    root = load_kgml(caminho_kgml)
    titulo_pathway = extract_pathway_title(caminho_kgml)
    retangulos = extract_graphical_elements_with_id(root)
    retangulos_enriquecidos = enrich_graphic_elements(retangulos, caminho_tsv, caminho_json_taxon,
                                                            caminho_json_kos)
    colors_dict = load_json_to_dict(caminho_json_colors)
    diretoria = 'KEGGCharter/info'
    taxa_colors = load_json_to_dict('./KEGGCharter/info/ko00680_taxa2colors.json')
    create_image_maps(taxa_colors, retangulos_enriquecidos,
                     imagem="./KEGGCharter/original_kegg_map.png", arquivo_saida=f"image_maps_{titulo_pathway}.html",
                     titulo_imagem=f"{titulo_pathway} Map", colors_dict=colors_dict)
