# %%
from itertools import groupby
from operator import attrgetter

import matplotlib.pyplot as plt
import seaborn as sns
from biosurfer.core.database import Database
from biosurfer.core.helpers import ExceptionLogger
from biosurfer.core.models import (ORF, Gene, Protein, ProteinFeature,
                                   Transcript)
from biosurfer.plots.plotting import IsoformPlot
from matplotlib.patches import Patch
from sqlalchemy import select

db = Database('sqlite:///../biosurfer/core/gencode.sqlite3')
session = db.get_session()
Gene.session = session

gene_list = (
    "ABCC10",
    "ABCF2",
    "AFTPH",
    "ARHGAP17",
    "ASPSCR1",
    "ATG7",
    "BMP4",
    "BOD1",
    "CALCOCO1",
    "CAPZB",
    "CD44",
    "CDK10",
    "CDK5RAP3",
    "CEP89",
    "CLCN7",
    "CLYBL",
    "CRELD2",
    "CRTC3",
    "CSNK1G3",
    "CTDSP2",
    "CYBRD1",
    "DHRS1",
    "DPP8",
    "DUSP3",
    "EEF1AKMT2",
    "EIF3L",
    "ELMOD3",
    "EPB41L2",
    "EPS15L1",
    "ETFA",
    "EXOSC10",
    "FAM118A",
    "FHL3",
    "FHOD1",
    "FIBP",
    "FIP1L1",
    "FLCN",
    "GAK",
    "GANAB",
    "GGT7",
    "GMDS",
    "GNB2",
    "GPS1",
    "GSS",
    "HLA-A",
    "HM13",
    "HNRNPUL1",
    "HYOU1",
    "ISYNA1",
    "JAG1",
    "KAT2A",
    "KLC1",
    "LETMD1",
    "LGALS3BP",
    "LMF2",
    "LSS",
    "LUC7L",
    "MACROD1",
    "MED16",
    "MLPH",
    "MMD",
    "MON1A",
    "MRC2",
    "MRM3",
    "NAGLU",
    "NCAPH2",
    "NCLN",
    "NECAB3",
    "NIBAN2",
    "NOC2L",
    "NPRL3",
    "NQO1",
    "NT5C2",
    "NUCB1",
    "OS9",
    "PCBP2",
    "PCNP",
    "PGM3",
    "PGS1",
    "PHKG2",
    "PICK1",
    "PNKD",
    "POLDIP2",
    "POLR2E",
    "PPP6R2",
    "PRKAR1A",
    "RANGAP1",
    "RECQL5",
    "RELA",
    "RIC8B",
    "RNPS1",
    "RPS6KA1",
    "RSAD1",
    "SCYL1",
    "SLC12A9",
    "SLC9A3R2",
    "SMARCD3",
    "SNRNP70",
    "SPAG9",
    "SPATA20",
    "SPG7",
    "SPPL2B",
    "SPRED2",
    "ST7L",
    "STIMATE",
    "SUN2",
    "TCEA3",
    "TGFBI",
    "TINF2",
    "TMED10",
    "TMEM175",
    "TMEM263",
    "TMEM94",
    "TNFAIP1",
    "TPCN2",
    "TPM2",
    "TRIP6",
    "TRMT61B",
    "TRPC4AP",
    "TSC2",
    "TSC22D4",
    "UBE2L3",
    "UBTF",
    "UNKL",
    "VAPA",
    "WLS",
    "YDJC",
    "YIPF3",
    "ZNF263",
    "ZNF304",
    "ZNF408",
    "ZNF584",
    "ZNF75A",
    "ZNF800",
    "ZSCAN32"
)

# %%
transcripts = set(session.execute(
    select(Transcript).
    select_from(ProteinFeature).
    join(ProteinFeature.protein).
    join(Protein.orf).
    join(ORF.transcript).
    join(Transcript.gene).
    where(Gene.name.in_(gene_list))
).scalars())

# %%
for gene, transcripts in groupby(sorted(transcripts, key=attrgetter('name')), key=attrgetter('gene')):
    isoplot = IsoformPlot(transcripts)
    with ExceptionLogger(f'{gene}'):
        domain_names = list({domain.name for tx in isoplot.transcripts for domain in tx.protein.features})
        cmap = sns.color_palette('pastel', len(domain_names), as_cmap=True)
        domain_colors = dict(zip(domain_names, cmap))
        isoplot.draw_all_isoforms()
        isoplot.draw_frameshifts()
        for track, tx in enumerate(isoplot.transcripts):
            if not tx.orfs:
                continue
            for domain in tx.protein.features:
                for i, (_, domain_exon) in enumerate(groupby(domain.residues, key=attrgetter('primary_exon'))):
                    domain_exon = list(domain_exon)
                    start = domain_exon[0].codon[1].coordinate
                    stop = domain_exon[-1].codon[1].coordinate
                    isoplot.draw_region(
                        track,
                        start = start,
                        stop = stop,
                        edgecolor = 'none',
                        facecolor = domain_colors[domain.name],
                        zorder = 1.8,
                        label = domain.name
                    )
                    # isoplot.draw_text((start+stop)//2, track, text=domain.name, color='w', ha='center', size='xx-small')
        isoplot._bax.big_ax.legend(
            handles = [Patch(color=color, label=name) for name, color in domain_colors.items()],
            ncol = 1,
            loc = 'center left',
            # mode = 'expand',
            bbox_to_anchor = (1.05, 0.5)
        )
        isoplot.fig.set_size_inches(10, 0.5*len(isoplot.transcripts))
        plt.savefig(f'../output/domains/{gene.name}_domains.png', facecolor='w', transparent=False, dpi=200, bbox_inches='tight')
    plt.close(isoplot.fig)
# %%
