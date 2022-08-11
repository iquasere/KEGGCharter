# KEGGCharter

A tool for representing genomic potential and transcriptomic expression into KEGG pathways.

* [Features](https://github.com/iquasere/KEGGCharter#features)
* [Installation](https://github.com/iquasere/KEGGCharter#installation)
* [Running KEGGCharter](https://github.com/iquasere/KEGGCharter#running-keggcharter)
* [Testing KEGGCharter](https://github.com/iquasere/KEGGCharter#testing-keggcharter)
* [Outputs](https://github.com/iquasere/KEGGCharter#outputs)
* [Arguments for KEGGCharter](https://github.com/iquasere/KEGGCharter#arguments-for-keggcharter)
* [Referencing KEGGCharter](https://github.com/iquasere/KEGGCharter#referencing-keggcharter)


## Features

KEGGCharter is a user-friendly implementation of KEGG API and Pathway functionalities. It allows for:
* Conversion of KEGG IDs to KEGG Orthologs (KO) and of KO to EC numbers
* Representation of the metabolic potential of the main taxa in KEGG metabolic maps (up to the top 10, each distinguished by its own colour)
* Representation of differential expression between samples in KEGG metabolic maps (the collective sum of each function will be represented)

## Installation

KEGGCharter can be easily installed with Bioconda.
```
conda install -c conda-forge -c bioconda keggcharter
```

## Running KEGGCharter

To run KEGGCharter, an input file must be supplied - see ```Testing KEGGCharter``` section - and columns 
with genomic and/or transcriptomic information, as well as one column with either KEGG IDs, KOs or EC numbers, must be 
present in the file and specified through the command line.
```
python kegg_charter.py -f input_file.xlsx -o output_folder -mgc mg_column1,mg_column2 -mtc mt_column1,mt_column2 ...
```

## Testing KEGGCharter

An example input file is available [here](https://github.com/iquasere/KEGGCharter/blob/master/MOSCA_Entry_Report.xlsx). 
This is one output of [MOSCA](https://github.com/iquasere/MOSCA), which can be directly inputted to KEGGCharter to obtain
metabolic representations by running:
```
kegg_charter.py -f MOSCA_Entry_Report.xlsx -gcol mg -tcol mt_0.01a_normalized,mt_1a_normalized,mt_100a_normalized,mt_0.01b_normalized,mt_1b_normalized,mt_100b_normalized,mt_0.01c_normalized,mt_1c_normalized,mt_100c_normalized -keggc "Cross-reference (KEGG)" -o test_keggcharter -tc "Taxonomic lineage (GENUS)"
```
Just make sure ```MOSCA_Entry_Report.xlsx``` is in the present folder, or indicate the path to it. This command will create
representations for all 252 default maps of KEGGCharter. If you want to represent for less or more, run with the ```--metabolic-maps``` 
parameter to indicate to KEGGCharter what maps to run (comma separated).

### First time KEGGCharter runs it will take a long time

KEGGCharter needs KGMLs and EC numbers to boxes relations, which it will automatically retrieve for every map inputted. 
This might take some time, but you only need to run it once. 

Default directory for storing these files is the folder containing the ```kegg_charter.py``` script, but it can be customized
with the ```--resources-directory``` parameter.

## Outputs

KEGGCharter produces a table from the inputed data with two new columns - KO (KEGG Charter) and EC number (KEGG Charter) - containing the results of conversion of KEGG IDs to KOs and KOs to EC numbers, respectively. This file is saved as ```KEGGCharter_results``` in the output directory. 
KEGGCharter then represents this information in KEGG metabolic maps. If information is available as result of (meta)genomics analysis, KEGGCharter will localize the boxes whose functions are present in the organisms' genomes, mapping their genomic potential. If (meta)transcriptomics data is available, KEGGCharter will consider the sample as a whole, measuring gene expression and performing a multi-sample comparison for each function in the metabolic maps.
* maps with genomic information are identified with the prefix "potential_" from genomic potential (figure 1).

![ScreenShot](potential_Methane_metabolism.png)
Figure 1 - KEGG metabolic map of methane metabolism, with identified taxa for each function from a simulated dataset.

* maps with transcriptomic information are identified with the prefix "differential_" from differential expression (figure 2).

![ScreenShot](differential_Methane_metabolism.png)
Figure 2 - KEGG metabolic map of methane metabolism, with differential analysis of quantified expression for each function from a simulated dataset.

## Arguments for KEGGCharter

KEGGCharter provides several options for customizing its workflow.
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

### Mock imputation of quantification and taxonomy

Sometimes, not all information required for KEGGCharter will be available. 
In this cases, KEGGCharter may use mock imputations of quantification and/or taxonomy.
To input mock quantification, run with the ```--input-quantification``` parameter. This will attribute a quantification 
of 1 to every protein in the input dataset.
This replaces the ```--genomic-columns``` parameter.
To input mock taxonomy, run with the ```--input-taxonomy [mock taxonomy]``` parameter, where ```[mock taxonomy]``` 
should be replaced with the value to be presented in the maps. This will attribute that taxonomic classification to 
every protein in the input dataset, which might be useful to, for example, represent "metagenome" in the genomic 
potential maps.
This replaces the ```--taxonomic-columns``` parameter.

```

## Referencing KEGGCharter

If you use KEGGCharter, please cite its [publication](https://www.sciencedirect.com/science/article/pii/S2001037022001179).