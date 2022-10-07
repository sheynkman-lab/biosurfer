import click
from biosurfer.analysis.load_gencode_database import check_database
from biosurfer.analysis.alignment_analysis_gencode import run_hybrid_alignment
from biosurfer.analysis.plot_biosurfer import run_plot




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
@click.option('--summary', is_flag=True, help="Prints summary statistics and plots for hybrid alignmentßß.")

def run_hybrid_al(verbose, db_name, o, output_path, summary):
    """ This script runs hybrid alignment on the provided database. """
    if o and db_name:
        click.echo('')
        click.echo('----- Running hybrid alignment: ', err=True)
        click.echo('')
        if summary:
            click.echo('----- with stats: ', err=True)
            run_hybrid_alignment(db_name, output_path, True)
        else: 
            click.echo('----- without stats: ', err=True)
            run_hybrid_alignment(db_name, output_path, False)


@cli.command("plot")
@click.option('--verbose', is_flag=True, help="Will print verbose messages.")
@click.option('--hal', is_flag=True, help="Plot hybrid alignment result.")

def plot_hal(verbose, hal):
    """ This script enables plotting functionality. """
    if hal:
        click.echo('')
        click.echo('----- Plotting hybrid alignment result: ', err=True)
        click.echo('')
        run_plot()