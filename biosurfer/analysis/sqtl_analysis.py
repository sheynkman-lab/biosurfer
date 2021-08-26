from itertools import filterfalse, tee
from operator import attrgetter
from typing import Iterable

from IPython.display import display

from biosurfer.core.alignments import pairwise_align_protein_sets, export_annotated_pblocks_to_tsv
from biosurfer.core.models import Chromosome, Gene, Junction, Transcript
from biosurfer.plots.plotting import IsoformPlot


def split_transcripts_on_junction_usage(junction: 'Junction', transcripts: Iterable['Transcript']):
    def uses_junction(transcript):
        return junction in transcript.junctions
    tx1, tx2 = tee(transcripts)
    transcripts_using = filter(uses_junction, tx1)
    transcripts_not_using = filterfalse(uses_junction, tx2)
    return set(transcripts_using), set(transcripts_not_using)

if __name__ == '__main__':
    gene = Gene.from_name('TANGO2')
    junc = gene.transcripts[1].junctions[3]
    # junc = Junction(69072023, 69072627, gene.chromosome, gene.strand)
    using, not_using = split_transcripts_on_junction_usage(junc, gene.transcripts)
    print(junc)
    print(f'using: {using}')
    print(f'not using: {not_using}')

    pblocks_containing_junc = []
    alns = pairwise_align_protein_sets((tx.orfs[0].protein for tx in not_using), (tx.orfs[0].protein for tx in using))
    for aln in alns:
        aln.annotate()
        up_exon, down_exon = aln.other.transcript.get_exons_from_junction(junc)
        pblocks_containing_junc.extend(pblock for pblock in aln.protein_blocks if {up_exon, down_exon} & pblock.other_exons)
    export_annotated_pblocks_to_tsv('test_sQTL_annotations.tsv', pblocks_containing_junc)
    print('done')

    isoplot = IsoformPlot(sorted(gene.transcripts, key=attrgetter('appris')))
    isoplot.draw_all_isoforms()
    isoplot.draw_frameshifts()
    isoplot.fig.set_size_inches(9, 0.5*len(gene.transcripts))
