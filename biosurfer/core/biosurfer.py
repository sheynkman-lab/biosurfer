import click
from example_scripts.load_gencode_database import check_database
from example_scripts.alignment_analysis_gencode import run_hybrid_alignment


@click.command()

@click.option('--verbose', is_flag=True, help="Will print verbose messages.")
#@click.argument('db_name', nargs=1, required=True, type=click.STRING)
@click.argument('filename', nargs=6, type=click.Path(exists=True), required=False)
@click.option('--db', is_flag=True, help="Creates database for the provided genocode files.")
@click.option('--alignment', is_flag=True, help="Run hybrid alignment script.")
@click.argument('db_name')

def cli(verbose, filename, db, alignment, db_name):
    """ Command line interface function to retrieve user inputs
    Args:
        verbose: Option to show more details on terminal
        filename: List to store all input file names
        db: Option to input data files
        alignment: Option to run hybrid alignment
        db_name: User input database name

    Returns:
      Nothing
    """
    if verbose:
        click.echo("We are in the verbose mode.")
    click.echo("""
--------------------------------------
         Welcome to Biosurfer!
--------------------------------------
""", err=True)
    click.echo()
    # TODO: Sanity check for Gencode file - from Gloria @Liz : subflags for all files
    # TODO: Break down creating database - Check Kalisto tool for example  
    if db:
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

    if alignment:
        click.echo('')
        click.echo('----- Running hybrid alignment: ', err=True)
        click.echo('')
        run_hybrid_alignment(db_name)

    # parser = argparse.ArgumentParser(description='Import gencode files.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # #parser.add_argument('filepath', type=argparse.FileType('r'), nargs='+', help='Path of a file or a folder of files.')
    # parser.add_argument("-i", "--input", help="File or directory to import", required=True)
    # parser.add_argument("-v", "--verbose", help="Display info status messages", action="store_true")
    # parser.add_argument("-q", "--quiet", help="Suppress most output", action="store_true")
    # parser.add_argument("--debug", help="Set logging to debug", action="store_true")
    # args = parser.parse_args()

    # logpath = None
    # if args.log:
    #     logpath = os.path.abspath(args.log)
    #     if os.path.isdir(logpath):
    #         logpath = os.path.join(logpath, "biosurfer.log")
    # else:
    #     logpath = ("biosurfer.log")

    # logger = logging.getLogger("Analysis") # create logger
    # sh = logging.StreamHandler()
    # fh = logging.FileHandler(logpath)
    # formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d,%H:%M:%S')
    # sh.setFormatter(formatter)
    # fh.setFormatter(formatter)

    # # set level based on args
    # if args.debug:
    #     sh.setLevel(logging.DEBUG)
    #     fh.setLevel(logging.DEBUG)
    # elif args.verbose:
    #     sh.setLevel(logging.INFO)
    #     fh.setLevel(logging.INFO)
    # elif args.quiet:
    #     sh.setLevel(logging.ERROR)
    #     fh.setLevel(logging.ERROR)
    # else:
    #     sh.setLevel(logging.WARNING)
    #     fh.setLevel(logging.WARNING)

    # logger.addHandler(sh) # add handler to logger
    # logger.addHandler(fh)

    # success = 0

    # if os.path.isdir(INPUT):
    #     for f in os.listdir(INPUT):
    #         if f.endswith(".gtf"):
    #             logging.info("** Importing {} **".format(f))
    #             gencode_gtf = os.path.join(INPUT, f)
    #             success = 1
    # elif os.path.isfile(INPUT):
    #     if INPUT.endswith(".gb") or f.endswith(".gbk"):
    #         logging.info("** Importing {} **".format(INPUT))
    #         mongo_import_genbank(INPUT, DB, "genes")
    #         success = 1