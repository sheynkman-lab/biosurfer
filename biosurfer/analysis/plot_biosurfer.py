# %%
from pathlib import Path
from more_itertools import partition
from biosurfer.core.alignments import ProteinAlignment
from biosurfer.core.constants import APPRIS
from biosurfer.core.database import Database
from biosurfer.core.helpers import (get_ids_from_gencode_fasta,
                                    get_ids_from_lrp_fasta,
                                    get_ids_from_pacbio_fasta, skip_gencode,
                                    skip_par_y)
from biosurfer.core.models.biomolecules import Gene, Transcript
from biosurfer.plots.plotting import IsoformPlot

#%%
def run_plot(output: Path, gene: str, db_name: str, transcript_ids: tuple[str]):
    """ Main plot function to invoke plotting for different pipelines/scripts.
    Args:
      Nothing
    Returns:
      Nothing
    """
    if not output:
      output = Path('.')
    db = Database(db_name)
    with db.get_session() as s:
        print(f'Loading transcripts from database...')
        
        if gene:
            gene_obj = Gene.from_name(s, gene)
            if gene_obj is None:
                print(f'Gene "{gene}" not found in database')
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
                print(f'Transcript ID "{tx_id}" not found in database')
            if transcript_ids:
                anchor = transcripts.get(transcript_ids[0], None)
            else:
                print('No isoforms provided')
                anchor = None
            others = [tx for tx in map(transcripts.get, found) if tx is not anchor]

        if anchor:
            print(f'Reference isoform: {anchor}')
            gene = anchor.gene.name

            alns: dict[Transcript, ProteinAlignment] = dict()
            for other in others:
                if anchor.protein is None or other.protein is None:
                    alns[other] = None
                else:
                    try:
                        alns[other] = ProteinAlignment.from_proteins(anchor.protein, other.protein)
                    except ValueError:
                        print(f'Could not plot isoform {other}')
            
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
            print(f'Saved {filepath}')

 


if __name__ == "__main__":
    run_plot(output, gene, db_name, transcript_ids)