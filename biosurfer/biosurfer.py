from operator import attrgetter
from pathlib import Path

import click
from more_itertools import partition

from biosurfer.analysis.alignment_analysis_gencode import run_hybrid_alignment
from biosurfer.core.alignments import ProteinAlignment
from biosurfer.core.constants import APPRIS
from biosurfer.core.database import Database
from biosurfer.core.helpers import (get_ids_from_gencode_fasta,
                                    get_ids_from_lrp_fasta,
                                    get_ids_from_pacbio_fasta, skip_gencode,
                                    skip_par_y)
from biosurfer.core.models.biomolecules import Gene, Transcript
from biosurfer.plots.plotting import IsoformPlot
from biosurfer.analysis.plot_biosurfer import run_plot

@click.group(chain=True)
def cli():
    """
    \b
--------------------------------------
         Welcome to Biosurfer!
--------------------------------------
    """

@cli.command("load_db")
@click.option('-v', '--verbose', is_flag=True, help="Will print verbose messages")
@click.option('-d', '--db_name', required=True, help="Database name")
@click.option('--source', type=click.Choice(['GENCODE', 'PacBio'], case_sensitive=False),  required=True, help="Source of input data")
@click.option('--gtf', required=True, type=click.Path(exists=True), help='Path to gtf file')
@click.option('--tx_fasta', required=True, type=click.Path(exists=True, path_type=Path), help='Path to transcript sequence fasta file')
@click.option('--tl_fasta', required=True, type=click.Path(exists=True, path_type=Path), help='Path to protein sequence fasta file')
@click.option('--sqanti', type=click.Path(exists=True, path_type=Path), help='Path to SQANTI classification tsv file')
def run_populate_database(verbose: bool, db_name: str, source :str, gtf: Path, tx_fasta: Path, tl_fasta: Path, sqanti: Path):
    """Loads transcript and protein isoform information from provided files into a Biosurfer database.
    A new database is created if the target database does not exist."""

    db = Database(db_name)
    if source == "GENCODE":
        click.echo('----- Loading database with reference ', err=True)
        db.load_gencode_gtf(gtf)
        db.load_transcript_fasta(tx_fasta, get_ids_from_gencode_fasta, skip_par_y)
        db.load_translation_fasta(tl_fasta, get_ids_from_gencode_fasta, skip_par_y)

    elif source == "PacBio":
        click.echo('----- Loading database without reference ', err=True)
        db.load_pacbio_gtf(gtf)
        db.load_transcript_fasta(tx_fasta, get_ids_from_pacbio_fasta)
        db.load_translation_fasta(tl_fasta, get_ids_from_lrp_fasta, skip_gencode)
        if sqanti:
            db.load_sqanti_classifications(sqanti)

# Run with 2>/dev/null to avoid printing error messages on terminal
@cli.command("hybrid_alignment")
@click.option('-v', '--verbose', is_flag=True, help="Will print verbose messages.")
@click.option('-d', '--db_name', required=True, nargs=1, help='Database name')
@click.option('-o', '--output', type=click.Path(exists=True, file_okay=False, writable=True, path_type=Path))
# @click.option('--summary', is_flag=True, help="Prints summary statistics and plots for hybrid alignment.")
def run_hybrid_al(verbose, db_name, output):
    """ This script runs hybrid alignment on the provided database. """
    click.echo('')
    click.echo('----- Running hybrid alignment: ', err=True)
    click.echo('')
    if not output:
        output = Path('.')
    run_hybrid_alignment(db_name, output)

@cli.command("plot")
@click.option('-v', '--verbose', is_flag=True, help="Print verbose messages")
@click.option('-o', '--output', type=click.Path(exists=True, file_okay=False, writable=True, path_type=Path), help='Directory in which to save plots')
@click.option('-d', '--db_name', required=True, nargs=1, help='Database name')
@click.option('--gene', type=str, help='Name of gene for which to plot all isoforms; overrides TRANSCRIPT_IDS')
@click.argument('transcript_ids', nargs=-1)
def plot_isoforms(verbose: bool, output: Path, gene: str, db_name: str, transcript_ids: tuple[str]):
    """Plot isoforms from a single gene, specified by TRANSCRIPT_IDS."""
    run_plot(output, gene, db_name, transcript_ids)

    