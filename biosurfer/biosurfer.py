from operator import attrgetter
from pathlib import Path
from more_itertools import partition
import click

from biosurfer.analysis.alignment_analysis_gencode import run_hybrid_alignment
from biosurfer.analysis.load_gencode_database import check_database
from biosurfer.analysis.plot_biosurfer import run_plot
from biosurfer.core.alignments import ProteinAlignment
from biosurfer.core.constants import APPRIS
from biosurfer.core.database import Database
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

@cli.command("create_database")

@click.option('--verbose', is_flag=True, help="Will print verbose messages.")
@click.argument('filename', nargs=6, type=click.Path(exists=True))
@click.option('--db', is_flag=True, required=True, help="Creates database for the provided genocode files.")
@click.argument('db_name')

def run_populate_database(verbose, filename, db, db_name):
    """ This script creates a database for the provided gencode files. """
    if db and db_name:
        click.echo('----- Input files:', err=True)
        click.echo(click.format_filename(filename[0]))
        click.echo(click.format_filename(filename[1]))
        click.echo(click.format_filename(filename[2]))
        click.echo(click.format_filename(filename[3]))
        click.echo(click.format_filename(filename[4]))
        click.echo(click.format_filename(filename[5]))
        click.echo('')
        click.echo('----- Creating database: ', err=True)
        check_database(click.format_filename(filename[0]),click.format_filename(filename[1]),click.format_filename(filename[2]),click.format_filename(filename[3]), click.format_filename(filename[4]),click.format_filename(filename[5]), db_name)

@cli.command("hybrid_alignment")
@click.option('--verbose', is_flag=True, help="Will print verbose messages.")
@click.option('--o', is_flag=True, help="Directory to write output to.")
@click.argument('output_path', type=click.Path(exists=True))
@click.argument('db_name')

def run_hybrid_al(verbose, db_name, o, output_path):
    """ This script runs hybrid alignment on the provided database. """
    if o and db_name:
        click.echo('')
        click.echo('----- Running hybrid alignment: ', err=True)
        click.echo('')
        run_hybrid_alignment(db_name, output_path)


@cli.command("plot")
@click.option('-v', '--verbose', is_flag=True, help="Print verbose messages")
@click.option('-o', '--output', type=click.Path(exists=True, file_okay=False, writable=True, path_type=Path), help='Directory in which to save plots')
@click.option('-d', '--db_name', required=True, nargs=1, help='Database name')
@click.option('--gene', type=str, help='Name of gene for which to plot all isoforms; overrides TRANSCRIPT_IDS')
@click.argument('transcript_ids', nargs=-1)
def plot_isoforms(verbose: bool, output: Path, gene: str, db_name: str, transcript_ids: tuple[str]):
    """Plot isoforms from a single gene, specified by TRANSCRIPT_IDS."""
    if output is None:
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
