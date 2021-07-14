#%%
# from exon import Exon
# from orf import ORF
# from gene import Gene
# pretend that ORF is all exons
# currently very simple gtf with transcript and exon features only
from models import Gene, ORF, Exon, session, Base, engine
def read_gtf_line(line: str) -> list:
    """Read and parse a single gtf line

    Args:
        line (str): unbroken line of a gtf file

    Returns:
        list: gtf attributes
            chromosome : str
            source : str
            feature : str
            start : int
            stop : int
            score : str
            strand : str
            phase : str
            attributes: dict

    """
    chromosome, source, feature, start, stop, score, strand, phase, attributes = line.split('\t')
    start = int(start)
    stop = int(stop)
    attributes = attributes.split(';')[:-1]
    attributes = [att.strip(' ').split(' ') for att in attributes]
    attributes = {att[0]: att[1].strip('"') for att in attributes}
    return [chromosome, source, feature, start, stop, score, strand, phase, attributes]

def load_data_from_gtf(gtf_file: str) -> None:
    Base.metadata.create_all(engine)
    genes = {}
    transcripts = {}
    exons = {}
    with open(gtf_file) as gtf:
        for line in gtf:
            chromosome, source, feature, start, stop, score, strand, phase, attributes = read_gtf_line(line)
            if feature == 'transcript':
                if attributes['gene_id'] not in genes.keys():
                    gene = Gene()
                    gene.name = attributes['gene_id']
                    gene.strand = strand
                    genes[attributes['gene_id']] = gene
                    session.add(gene)
                orf = ORF()
                orf.name = attributes['transcript_id']
                genes[attributes['gene_id']].orfs.append(orf)
                transcripts[attributes['transcript_id']] = orf
                session.add(orf)
                
            elif feature == 'exon':
                if (chromosome, start, stop) in exons.keys():
                    exon = exons[(chromosome, start, stop)]
                else:
                    exon = Exon()
                    exon.start = start
                    exon.stop = stop
                exon.orfs.append(transcripts[attributes['transcript_id']])
                session.add(exon)
                exons[(chromosome, start, stop)] = exon
    session.commit() #Attempt to commit all the records
    



load_data_from_gtf('/Users/bj8th/Documents/Sheynkman-Lab/Data/test/jurkat_chr22_corrected.gtf')

# %%
