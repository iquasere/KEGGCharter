# KEGGCharter

A tool for representing genomic potential and transcriptomic expression into KEGG pathways.

## Features

KEGGCharter is a user-friendly implementation of KEGG API and Pathway functionalities. It allows for:
* Conversion of KEGG IDs to KEGG Orthologs (KO) and of KO to EC numbers
* Representation of the metabolic potential of the main taxa in KEGG metabolic maps (up to the top 10, each distinguished by its own colour)
* Representation of differential expression between samples in KEGG metabolic maps (the collective sum of each function will be represented)

## Installation

KEGGCharter can be easily installed with Bioconda
```
conda install -c conda-forge -c bioconda keggcharter
```

## Usage

KEGGCharter needs an input file, but that is all it needs!
```
usage: kegg_charter.py [-h] [-o OUTPUT] [--tsv] [-t THREADS]
                       [-mm METABOLIC_MAPS] [-gcol GENOMIC_COLUMNS]
                       [-tcol TRANSCRIPTOMIC_COLUMNS] [-tls TAXA_LIST]
                       [-not NUMBER_OF_TAXA] [-keggc KEGG_COLUMN]
                       [-koc KO_COLUMN] [-ecc EC_COLUMN] [--resume] [-v] -f
                       FILE -tc TAXA_COLUMN [--show-available-maps]

KEGGCharter - A tool for representing genomic potential and transcriptomic
expression into KEGG pathways

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output directory
  --tsv                 Results will be outputed in TSV format (and not
                        EXCEL).
  -mm METABOLIC_MAPS, --metabolic-maps METABOLIC_MAPS
                        IDs of metabolic maps to output
  -gcol GENOMIC_COLUMNS, --genomic-columns GENOMIC_COLUMNS
                        Names of columns with genomic identification
  -tcol TRANSCRIPTOMIC_COLUMNS, --transcriptomic-columns TRANSCRIPTOMIC_COLUMNS
                        Names of columns with transcriptomics quantification
  -tls TAXA_LIST, --taxa-list TAXA_LIST
                        List of taxa to represent in genomic potential charts
                        (comma separated)
  -not NUMBER_OF_TAXA, --number-of-taxa NUMBER_OF_TAXA
                        Number of taxa to represent in genomic potential
                        charts (comma separated)
  -keggc KEGG_COLUMN, --kegg-column KEGG_COLUMN
                        Column with KEGG IDs.
  -koc KO_COLUMN, --ko-column KO_COLUMN
                        Column with KOs.
  -ecc EC_COLUMN, --ec-column EC_COLUMN
                        Column with EC numbers.
  --resume              Data inputed has already been analyzed by KEGGCharter.
  -v, --version         show program's version number and exit
  -tc TAXA_COLUMN, --taxa-column TAXA_COLUMN
                        Column with the taxa designations to represent with
                        KEGGChart

required named arguments:
  -f FILE, --file FILE  TSV or EXCEL table with information to chart

Special functions:
  --show-available-maps
                        Outputs KEGG maps IDs and descriptions to the console
                        (so you may pick the ones you want!)
```

To run KEGGCharter, an input file must be supplied - see "Example" section - and the columns with genomic and/or transcriptomic information as well. Output directory is not mandatory, but may help find results.
```
python kegg_charter.py -f input_file.xlsx -o output_folder -mgc mg_column1,mg_column2 -mtc mt_column1,mt_column2 ...
```

## Example

An example input file is available when downloading the GitHub repository. Inserting the KEGG IDs and genomic and transcriptomic quantifications in this file will allow to use KEGGCharter with no errors... in principle.

## Outputs

KEGGCharter produces a table from the inputed data with two new columns - KO (KEGG Charter) and EC number (KEGG Charter) - containing the results of conversion of KEGG IDs to KOs and KOs to EC numbers, respectively. This file is saved as KEGGCharter_results in the output directory. 
KEGGCharter then represents this information in KEGG metabolic maps. If information is available as result of (meta)genomics analysis, KEGGCharter will localize the boxes whose functions are present in the organisms' genomes, mapping their genomic potential. If (meta)transcriptomics data is available, KEGGCharter will consider the sample as a whole, measuring gene expression and performing a multi-sample comparison for each function in the metabolic maps.
* maps with genomic information are identified with the prefix "potential_" from genomic potential (figure 1).

![ScreenShot](potential_Methane_metabolism.png)
Figure 1 - KEGG metabolic map of methane metabolism, with identified taxa for each function from a simulated dataset.

* maps with transcriptomic information are identified with the prefix "differential_" from differential expression (figure 2).

![ScreenShot](differential_Methane_metabolism.png)
Figure 2 - KEGG metabolic map of methane metabolism, with differential analysis of quantified expression for each function from a simulated dataset.