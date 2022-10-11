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
# @click.option('--ref', is_flag=True, help='Load reference isoforms')
@click.option('--source', type=click.Choice(['GENCODE', 'PacBio'], case_sensitive=False),  required=True, help="Source of input data")
@click.option('--gtf', required=True, type=click.Path(exists=True), help='Path to gtf file')
@click.option('--tx_fasta', required=True, type=click.Path(exists=True, path_type=Path), help='Path to transcript sequence fasta file')
@click.option('--tl_fasta', required=True, type=click.Path(exists=True, path_type=Path), help='Path to protein sequence fasta file')
@click.option('--sqanti', type=click.Path(exists=True, path_type=Path), help='Path to SQANTI classification tsv file')
def run_populate_database(verbose: bool, db_name: str, source, gtf: Path, tx_fasta: Path, tl_fasta: Path, sqanti: Path):
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

@cli.command("hybrid_alignment")
@click.option('-v', '--verbose', is_flag=True, help="Will print verbose messages.")
@click.option('-d', '--db_name', required=True, nargs=1, help='Database name')
@click.option('-o', '--output', type=click.Path(exists=True, file_okay=False, writable=True, path_type=Path))
@click.option('--summary', is_flag=True, help="Prints summary statistics and plots for hybrid alignmentßß.")
def run_hybrid_al(verbose, db_name, output, summary):
    """ This script runs hybrid alignment on the provided database. """
    click.echo('')
    click.echo('----- Running hybrid alignment: ', err=True)
    click.echo('')
    if summary:
        click.echo('----- with stats: ', err=True)
    else: 
        click.echo('----- without stats: ', err=True)
    if not output:
        output = Path('.')
    run_hybrid_alignment(db_name, output, summary)

@cli.command("plot")
@click.option('-v', '--verbose', is_flag=True, help="Print verbose messages")
@click.option('-o', '--output', type=click.Path(exists=True, file_okay=False, writable=True, path_type=Path), help='Directory in which to save plots')
@click.option('-d', '--db_name', required=True, nargs=1, help='Database name')
@click.option('--gene', type=str, help='Name of gene for which to plot all isoforms; overrides TRANSCRIPT_IDS')
@click.argument('transcript_ids', nargs=-1)
def plot_isoforms(verbose: bool, output: Path, gene: str, db_name: str, transcript_ids: tuple[str]):
    """Plot isoforms from a single gene, specified by TRANSCRIPT_IDS."""
    if not output:
        output = Path('.')
    db = Database(db_name)
    with db.get_session() as s:
        if verbose:
            click.echo(f'Loading transcripts from database...')
        
        if gene:
            gene_obj = Gene.from_name(s, gene)
            if gene_obj is None:
                click.echo(f'Gene "{gene}" not found in database', err=True)
                transcripts = dict()
                anchor = None
            else:
                transcripts = {tx.accession: tx for tx in gene_obj.transcripts}
                anchor = max(transcripts.values(), key=lambda tx: getattr(tx, 'appris', APPRIS.NONE))
            others = [tx for tx in transcripts.values() if tx is not anchor]
        else:
            transcripts: dict[str, Transcript] = {tx.accession: tx for tx in Transcript.from_accessions(s, transcript_ids).values()}
            not_found, found = partition(lambda tx_id: tx_id in transcripts, transcript_ids)
            for tx_id in not_found:
                click.echo(f'Transcript ID "{tx_id}" not found in database', err=True)
            if transcript_ids:
                anchor = transcripts.get(transcript_ids[0], None)
            else:
                if verbose:
                    click.echo('No isoforms provided')
                anchor = None
            others = [tx for tx in map(transcripts.get, found) if tx is not anchor]

        if anchor:
            if verbose:
                click.echo(f'Reference isoform: {anchor}')
            gene = anchor.gene.name

            alns: dict[Transcript, ProteinAlignment] = dict()
            for other in others:
                if anchor.protein is None or other.protein is None:
                    alns[other] = None
                else:
                    try:
                        alns[other] = ProteinAlignment.from_proteins(anchor.protein, other.protein)
                    except ValueError:
                        click.echo(f'Could not plot isoform {other}', err=True)
            
            filename = f'{db_name}-{gene}.png'
            plot = IsoformPlot([anchor] + list(alns.keys()))
            plot.draw_all_isoforms()
            plot.draw_frameshifts()
            for other, aln in alns.items():
                if aln:
                    plot.draw_protein_alignment_blocks(aln.blocks, anchor.protein, other.protein)
            plot.draw_legend()
            filepath = str(output/filename)
            plot.savefig(filepath)
            if verbose:
                click.echo(f'Saved {filepath}')
