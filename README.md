
# Biosurfer

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)  [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7182809.svg)](https://doi.org/10.5281/zenodo.7182809)



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

### Biosurfer command line:
```
Usage: biosurfer [OPTIONS] COMMAND1 [ARGS]... [COMMAND2 [ARGS]...]...
          
Options:
  --help  Show this message and exit.

Commands:
  hybrid_alignment  This script runs hybrid alignment on the...
  load_db           Loads transcript and protein isoform...
  plot              Plot isoforms from a single gene...
```
* Download the [toy gencode data](https://zenodo.org/record/7182809) from Zenodo into the project directory.

### Load database:

```
Usage: biosurfer load_db [OPTIONS]

  Loads transcript and protein isoform information from
  provided files into a Biosurfer database. A new database is
  created if the target database does not exist.

Options:
  -v, --verbose              Will print verbose messages
  -d, --db_name TEXT         Database name  [required]
  --source [GENCODE|PacBio]  Source of input data  [required]
  --gtf PATH                 Path to gtf file  [required]
  --tx_fasta PATH            Path to transcript sequence fasta
                             file  [required]
  --tl_fasta PATH            Path to protein sequence fasta
                             file  [required]
  --sqanti PATH              Path to SQANTI classification tsv
                             file
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

### Hybrid alignment
* Run hybdrid alignment script on the created database. Create a directory to store the output files.

```shell
biosurfer hybrid_alignment -d gencode_toy -o output/gencode_toy --summary
```

```
Usage: biosurfer hybrid_alignment [OPTIONS]

  This script runs hybrid alignment on the provided database.

Options:
  -v, --verbose           Will print verbose messages.
  -d, --db_name TEXT      Database name  [required]
  -o, --output DIRECTORY
  --summary               Prints summary statistics and plots
                          for hybrid alignment.
  --help                  Show this message and exit.
```

### Plots
* To plot the hybrid alignment result, run the following snippet.

```shell
biosurfer plot --db_name gencode_toy --gene CRYBG2 -o output/gencode_toy
```

```
Usage: biosurfer plot [OPTIONS] [TRANSCRIPT_IDS]...

  Plot isoforms from a single gene, specified by TRANSCRIPT_IDS.

Options:
  -v, --verbose           Print verbose messages
  -o, --output DIRECTORY  Directory in which to save plots
  -d, --db_name TEXT      Database name  [required]
  --gene TEXT             Name of gene for which to plot all
                          isoforms; overrides TRANSCRIPT_IDS
  --help                  Show this message and exit.
```

### Input

Biosurfer takes the following gencode files as input. The data from these files are loaded into a local database to run scripts.
NOTE: The files must be in the same order as provided below.
1. Gene annotation file (GTF)
2. Transcript reference sequence file (FASTA)
3. Translation reference sequence file (FASTA)

