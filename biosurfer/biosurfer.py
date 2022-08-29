import click
from example_scripts.load_gencode_database import check_database
from example_scripts.alignment_analysis_gencode import run_hybrid_alignment



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
@click.option('--db', is_flag=True, help="Creates database for the provided genocode files.")
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
@click.argument('db_name')
@click.option('--o', is_flag=True, help="Directory to write output to.")
@click.argument('output')

def run_hybrid_al(verbose, db_name, o, output):
    """ This script runs hybrid alignment on the provided database. """
    if o and db_name:
        click.echo('')
        click.echo('----- Running hybrid alignment: ', err=True)
        click.echo('')
        run_hybrid_alignment(db_name, output)