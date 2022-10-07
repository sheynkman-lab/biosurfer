from operator import attrgetter
from pathlib import Path
from typing import Iterable
import click
from biosurfer.analysis.load_gencode_database import check_database
from biosurfer.analysis.alignment_analysis_gencode import run_hybrid_alignment
from biosurfer.analysis.plot_biosurfer import run_plot
from biosurfer.core.alignments import ProteinAlignment
from biosurfer.core.database import Database
from biosurfer.core.models.biomolecules import Gene
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
@click.option('-v', '--verbose', is_flag=True, help="Will print verbose messages.")
@click.option('-o', '--output', type=click.Path(exists=True, file_okay=False, writable=True), help='Directory in which to save plots.')
@click.argument('db_name', nargs=1)
@click.argument('gene_list', nargs=-1)
def plot_hal(verbose: bool, output: Path, db_name: str, gene_list: Iterable[str]):
    """Generate plots for a list of genes/transcripts."""
    if output is None:
        output = Path('.')
    
    db = Database(db_name)
    with db.get_session() as s:
        if verbose:
            click.echo(f'Loading genes from database...')
        genes = Gene.from_names(s, gene_list)
        not_found = sorted(set(gene_list) - set(genes))
        for gene in not_found:
            click.echo(f'Gene "{gene}" not found in database', err=True)
        
        for gene_name, gene in genes.items():
            txs = sorted(gene.transcripts, key=attrgetter('appris'), reverse=True)
            anchor = txs[0]
            if anchor.protein is None:
                click.echo(f'Cannot compare isoforms of gene {gene} against non-coding reference', err=True)
                continue

            alns = dict()
            for other in txs[1:]:
                if other.protein is None:
                    alns[other] = None
                else:
                    try:
                        alns[other] = ProteinAlignment.from_proteins(anchor.protein, other.protein)
                    except ValueError:
                        click.echo(f'Could not plot isoform {other}', err=True)
            filename = f'{gene_name}.png'
            plot = IsoformPlot([anchor] + list(alns.keys()))
            plot.draw_all_isoforms()
            plot.draw_frameshifts()
            plot.draw_legend()
            filepath = str(output/filename)
            plot.savefig(filepath)
            if verbose:
                click.echo(f'Saved {filepath}')
