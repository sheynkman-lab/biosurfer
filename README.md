
# Biosurfer

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)  [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7297008.svg)](https://doi.org/10.5281/zenodo.7297008)


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


Clone the project repository (using SSH if need be) and create a [new conda environment](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands) if needed. 

```
# Clone the repository
git clone https://github.com/sheynkman-lab/biosurfer
    
# Move to the folder
cd biosurfer
    
# Run setup 
pip install --editable .
``` 

## Usage

### Biosurfer command line options:
```
Usage: biosurfer [OPTIONS] COMMAND1 [ARGS]... [COMMAND2 [ARGS]...]...

Options:
  --help  Show this message and exit.

Commands:
  hybrid_alignment  This script runs hybrid alignment on the provided...
  load_db           Loads transcript and protein isoform information from...
  plot              Plot isoforms from a single gene, specified by...
```
* Download the [toy gencode data](https://zenodo.org/record/7297008/files/biosurfer_gencode_toy_data.zip?download=1) from Zenodo into the project directory.

### 1. Load database:

```
Usage: biosurfer load_db [OPTIONS]

  Loads transcript and protein isoform information from provided files into a
  Biosurfer database. A new database is created if the target database does
  not exist.

Options:
  -v, --verbose              Will print verbose messages
  -d, --db_name TEXT         Database name  [required]
  --source [GENCODE|PacBio]  Source of input data  [required]
  --gtf PATH                 Path to gtf file  [required]
  --tx_fasta PATH            Path to transcript sequence fasta file
                             [required]
  --tl_fasta PATH            Path to protein sequence fasta file  [required]
  --sqanti PATH              Path to SQANTI classification tsv file (only for
                             PacBio isoforms)
  --help                     Show this message and exit.
```

 #### Load database using GENCODE reference (toy version) 
 
```bash
biosurfer load_db --source=GENCODE --gtf biosurfer_gencode_toy_data/gencode.v38.toy.gtf --tx_fasta biosurfer_gencode_toy_data/gencode.v38.toy.transcripts.fa --tl_fasta biosurfer_gencode_toy_data/gencode.v38.toy.translations.fa --db_name gencode_toy
``` 
Running GENCODE files without ```--ref``` will 
 #### Load database using PacBio data without reference (WTC11 data)
 
```bash
biosurfer load_db --source=PacBio --gtf biosurfer_wtc11_data/wtc11_with_cds.gtf --tx_fasta biosurfer_wtc11_data/wtc11_corrected.fasta  --tl_fasta biosurfer_wtc11_data/wtc11_orf_refined.fasta --sqanti biosurfer_wtc11_data/wtc11_classification.txt --db_name wtc11_db
``` 

### 2. Hybrid alignment
* Run hybdrid alignment script on the created database. Create a directory to store the output files.

```shell
biosurfer hybrid_alignment -d gencode_toy -o output/gencode_toy -- gencode
```

```
Usage: biosurfer hybrid_alignment [OPTIONS]

  This script runs hybrid alignment on the provided database.

Options:
  -v, --verbose           Print verbose messages
  -d, --db_name TEXT      Database name  [required]
  -o, --output DIRECTORY  Directory for output files
  --gencode               Also compare all GENCODE isoforms of a gene against
                          its anchor isoform
  --anchors FILE          TSV file with gene names in column 1 and anchor
                          isoform IDs in column 2
  --help                  Show this message and exit.
```
> Please note that in the code, the terms *`anchor`* and *`other`* correspond to the *`reference`* and *`alternative`* isoforms mentioned in the manuscript.

### 3. Visualize protein isoforms
* To visualization isoforms of *CRYBG2* gene, run the following snippet.

```shell
biosurfer plot -d gencode_toy --gene CRYBG2
```

```
Usage: biosurfer plot [OPTIONS] [TRANSCRIPT_IDS]...

  Plot isoforms from a single gene, specified by TRANSCRIPT_IDS.

Options:
  -v, --verbose           Print verbose messages
  -o, --output DIRECTORY  Directory in which to save plots
  -d, --db_name TEXT      Database name  [required]
  --gene TEXT             Name of gene for which to plot all isoforms;
                          overrides TRANSCRIPT_IDS
  --help                  Show this message and exit.
```


