
# Biosurfer

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)  [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7004071.svg)](https://doi.org/10.5281/zenodo.7004071)

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

#### Local building (without installation)


Clone the project repository and create a [new conda environment](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands) if needed.

```
# Clone the repository
git clone https://github.com/sheynkman-lab/biosurfer
    
# Move to the folder
cd biosurfer
    
# Run setup 
pip install --editable .
``` 

## Usage
* Download the [toy gencode data](https://zenodo.org/record/7004071) from Zenodo into the project directory.

* Run the following command to create a database and run hybrid alignment

```
biosurfer \
    --db data/gencode/gencode.v41.annotation.gtf data/gencode/gencode.v41.pc_transcripts.fa data/gencode/gencode.v41.pc_translations.fa data/gencode/grch38-protein-features.tsv data/gencode/pfamA.tsv data/gencode/prosite.dat gencode_v41 \
    --alignment \
    --o output/alignment-analysis/gencode_v41 
```
```
Usage: biosurfer [OPTIONS] [FILENAME]... DB_NAME OUTPUT

Options:
  --verbose    Will print verbose messages.
  --db         Creates database for the provided genocode files.
  --alignment  Run hybrid alignment script.
  --o          Directory to write output to.
  --help       Show this message and exit. 
```

### Input

Biosurfer takes the following gencode files as input. The data from these files are loaded into a local database to run scripts.
NOTE: The files must be in the same order as provided below.
1. Gene annotation file (GTF)
2. Transcript reference sequence file (FASTA)
3. Translation reference sequence file (FASTA)
4. grch38 protein feature file (TSV)
5. Protein Family mapping file (TSV)
6. PROSITE pattern data file    


    
## 3. References

* D’Antonio,M. et al. (2015) ASPicDB: a database web tool for alternative splicing analysis. Methods Mol. Biol., 1269, 365–378.
* Frankish,A. et al. (2019) GENCODE reference annotation for the human and mouse genomes. Nucleic Acids Res., 47, D766–D773.
* de la Fuente,L. et al. (2020)  21, 119.
* Gal-Oz,S.T. et al. (2021) DoChaP: thtappAS: a comprehensive computational framework for the analysis of the functional impact of differential splicing. Genome Biol.,e domain change presenter. Nucleic Acids Res., 49, W162–W168.
* Louadi,Z. et al. (2021) DIGGER: exploring the functional role of alternative splicing in protein interactions. Nucleic Acids Res., 49, D309–D318.
