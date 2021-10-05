# %%
import matplotlib.pyplot as plt
from biosurfer.core.database import Database
from biosurfer.core.helpers import ExceptionLogger
from biosurfer.core.models import (ORF, GencodeTranscript, Gene, PacBioTranscript, Protein, ProteinFeature,
                                   Transcript)
from biosurfer.plots.plotting import IsoformPlot
from sqlalchemy import select

db = Database('bone')
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
genes = set(session.execute(
    select(Gene).
    select_from(ProteinFeature).
    join(ProteinFeature.protein).
    join(Protein.orf).
    join(ORF.transcript).
    join(Transcript.gene).
    where(Gene.name.in_(gene_list))
).scalars())

# %%
for gene in genes:
    gc_transcripts = [tx for tx in gene.transcripts if isinstance(tx, GencodeTranscript) if tx.protein]
    pb_transcripts = [tx for tx in gene.transcripts if isinstance(tx, PacBioTranscript) if tx.protein]
    isoplot = IsoformPlot(gc_transcripts + pb_transcripts, track_spacing=1.0)
    fig_path = f'../output/domains/{gene.name}_domains.svg'
    with ExceptionLogger(f'{gene.name}'):
        isoplot.draw_all_isoforms()
        isoplot.draw_frameshifts()
        isoplot.draw_domains()
        isoplot.draw_legend()
        isoplot.fig.set_size_inches(10, 0.75*len(isoplot.transcripts))
        print(f'Saving {fig_path}...')
        plt.savefig(fig_path, facecolor='w', transparent=False, dpi=200, bbox_inches='tight')
    plt.close(isoplot.fig)
# %%
