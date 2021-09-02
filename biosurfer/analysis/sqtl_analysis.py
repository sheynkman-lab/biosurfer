from itertools import filterfalse, tee
from operator import attrgetter
from os.path import isfile
from typing import Iterable
from warnings import filterwarnings

import matplotlib.pyplot as plt
import pandas as pd
from biosurfer.core.alignments import (export_annotated_pblocks_to_tsv,
                                       pairwise_align_protein_sets)
from biosurfer.core.models import Chromosome, Gene, Junction, Transcript
from biosurfer.plots.plotting import IsoformPlot
from IPython.display import display
from matplotlib._api.deprecation import MatplotlibDeprecationWarning

filterwarnings("ignore", category=MatplotlibDeprecationWarning)


def split_transcripts_on_junction_usage(junction: 'Junction', transcripts: Iterable['Transcript']):
    def contains_both_splice_sites(transcript):
        return (transcript.start <= junction.donor <= transcript.stop and
                transcript.start <= junction.acceptor <= transcript.stop)
    def uses_junction(transcript):
        return junction in transcript.junctions
    tx1, tx2 = tee(filter(contains_both_splice_sites, transcripts))
    transcripts_using = filter(uses_junction, tx1)
    transcripts_not_using = filterfalse(uses_junction, tx2)
    return set(transcripts_using), set(transcripts_not_using)

if __name__ == '__main__':
    data_dir = '../../data/bone'
    output_dir = '../../output/bone'
    sqtls = pd.read_csv(f'{data_dir}/fibroblast_coloc_sqtls.csv')
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
            def pblock_is_related_to_junc(pblock):
                return up_exon in pblock.other_exons or down_exon in pblock.other_exons
            pblocks_containing_junc.extend(filter(pblock_is_related_to_junc, aln.protein_blocks))
        export_annotated_pblocks_to_tsv(f'{output_dir}/{gene.name}_{junc.donor}_{junc.acceptor}.tsv', pblocks_containing_junc)

        fig_path = f'{output_dir}/{gene.name}_{junc.donor}_{junc.acceptor}.png'
        if not isfile(fig_path):
            isoplot = IsoformPlot(using + not_using)
            isoplot.draw_all_isoforms()
            isoplot.draw_frameshifts()
            isoplot.draw_background_rect(start=junc.donor, stop=junc.acceptor, facecolor='#ffffb7')
            isoplot.fig.set_size_inches(9, 0.5*len(gene.transcripts))
            plt.savefig(fig_path, facecolor='w', transparent=False, dpi=300, bbox_inches='tight')
            print('\tsaved '+fig_path)
            plt.close(isoplot.fig)
