
# Biosurfer: Modeling connections between genomic, transcriptomic, and proteomic layers to characterize isoform function

"Surf" the biological network, from genome to transcriptome to proteome and back to gain insights into human disease biology.


**Contents**

- [Installation](#installation)
- [Features](#features)
- [References](#references)

## Installation
 

#### Building Requirements

* Python 3.9 or higher 
* Python packages (numpy, more-itertools, intervaltree, biopython, attrs, tqdm)
* Database (sqlalchemy >=1.4)
* Vizualization (matplotlib, brokenaxes)


#### Toy dataset

#### Local building (without installation)

Clone the project repository and create a new conda environment if needed.

    git clone https://github.com/sheynkman-lab/biosurfer
    cd biosurfer
    pip install --editable .
    
The dependency on [`graph-tool`](https://graph-tool.skewed.de/) currently requires a separate installation step. Run `conda install -c conda-forge graph-tool`.

#### Usage
Download dataset from Zenodo 

## 2. Features

#### Data sources and preprocessing
- External reference databases
    - GENCODE -- genes, transcripts, ORFs, proteins
    - Ensembl Pfam mappings (via BioMart) -- protein features
 

#### Software implementation details

- SQLAlchemy ORM creates an internal database using the input files and associates related Python objects with one another.
- describing the hybrid alignment algorithm?
- briefly describe implementation of plotting package?

#### Analysis-specific methods
- describe pblock similarity metrics
- sQTL functional analysis -- describe metrics?... (it doesn’t feel like we really used them much)
- describe statistical methodology for proteome-wide analysis of frameshift regions

    
## 3. References

* D’Antonio,M. et al. (2015) ASPicDB: a database web tool for alternative splicing analysis. Methods Mol. Biol., 1269, 365–378.
* Frankish,A. et al. (2019) GENCODE reference annotation for the human and mouse genomes. Nucleic Acids Res., 47, D766–D773.
* de la Fuente,L. et al. (2020)  21, 119.
* Gal-Oz,S.T. et al. (2021) DoChaP: thtappAS: a comprehensive computational framework for the analysis of the functional impact of differential splicing. Genome Biol.,e domain change presenter. Nucleic Acids Res., 49, W162–W168.
* Louadi,Z. et al. (2021) DIGGER: exploring the functional role of alternative splicing in protein interactions. Nucleic Acids Res., 49, D309–D318.
