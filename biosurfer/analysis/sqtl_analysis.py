from itertools import filterfalse, tee
from operator import attrgetter
from typing import Iterable
from warnings import filterwarnings

from IPython.display import display
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib._api.deprecation import MatplotlibDeprecationWarning

from biosurfer.core.alignments import pairwise_align_protein_sets, export_annotated_pblocks_to_tsv
from biosurfer.core.models import Chromosome, Gene, Junction, Transcript
from biosurfer.plots.plotting import IsoformPlot

filterwarnings("ignore", category=MatplotlibDeprecationWarning)

def split_transcripts_on_junction_usage(junction: 'Junction', transcripts: Iterable['Transcript']):
    def uses_junction(transcript):
        return junction in transcript.junctions
    tx1, tx2 = tee(transcripts)
    transcripts_using = filter(uses_junction, tx1)
    transcripts_not_using = filterfalse(uses_junction, tx2)
    return set(transcripts_using), set(transcripts_not_using)

if __name__ == '__main__':
    working_dir = '../../data/bone'
    sqtls = pd.read_csv(f'{working_dir}/fibroblast_coloc_sqtls.csv')
    sqtls = sqtls[sqtls['gene_type'] == 'protein_coding']
    genes = {gene.name: gene for gene in Gene.query.where(Gene.name.in_(sqtls['gene_name']))}

    for index, row in sqtls.iterrows():
        gene = genes[row['gene_name']]
        
        junc = Junction(row['event_start'], row['event_end'], Chromosome.from_name(row['chr']), gene.strand)
        using, not_using = split_transcripts_on_junction_usage(junc, gene.transcripts)
        if not all((using, not_using)):
            continue
        using = sorted(using, key=attrgetter('appris'))
        not_using = sorted(not_using, key=attrgetter('appris'))
        print(junc)
        print(f'\tusing: {using}')
        print(f'\tnot using: {not_using}')

        pblocks_containing_junc = []
        alns = pairwise_align_protein_sets((tx.orfs[0].protein for tx in not_using), (tx.orfs[0].protein for tx in using))
        for aln in alns:
            aln.annotate()
            up_exon, down_exon = aln.other.transcript.get_exons_from_junction(junc)
            pblocks_containing_junc.extend(pblock for pblock in aln.protein_blocks if {up_exon, down_exon} & pblock.other_exons)
        export_annotated_pblocks_to_tsv(f'{working_dir}/{gene.name}_{junc.donor}_{junc.acceptor}.tsv', pblocks_containing_junc)

        isoplot = IsoformPlot(using + not_using)
        isoplot.draw_all_isoforms()
        isoplot.draw_frameshifts()
        isoplot.draw_background_rect(start=junc.donor, stop=junc.acceptor, facecolor='#ffffb7')
        isoplot.fig.set_size_inches(9, 0.5*len(gene.transcripts))
        fig_path = f'{working_dir}/{gene.name}_{junc.donor}_{junc.acceptor}.png'
        plt.savefig(fig_path, facecolor='w', transparent=False, dpi=300, bbox_inches='tight')
        print('\tsaved '+fig_path)
        plt.close(isoplot.fig)
