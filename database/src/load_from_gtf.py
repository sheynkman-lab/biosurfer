#%%
# from exon import Exon
# from orf import ORF
# from gene import Gene
# pretend that ORF is all exons
# currently very simple gtf with transcript and exon features only
from models import Chromosome, Gene, Transcript, Exon
from database import db_session, Base, engine
from sqlalchemy import select
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
    chromosomes = {}
    with open(gtf_file) as gtf:
        for line in gtf:
            chr, source, feature, start, stop, score, strand, phase, attributes = read_gtf_line(line)
            if feature == 'gene':
                if chr not in chromosomes.keys():
                  chromosome = Chromosome()
                  chromosome.name = chr
                  chromosomes[chr] = chromosome
                else:
                    chromosome = chromosomes[chr]

                gene = Gene()
                gene.name = attributes['gene_id']
                gene.strand = strand
                genes[attributes['gene_id']] = gene
                db_session.add(gene)

                chromosome.genes.append(gene)
                db_session.add(chromosome)
            elif feature == 'transcript':
                transcript = Transcript()
                transcript.name = attributes['transcript_id']
                genes[attributes['gene_id']].transcripts.append(transcript)
                transcripts[attributes['transcript_id']] = transcript
                db_session.add(transcript)
                
            elif feature == 'exon':
                if (chromosome, start, stop) in exons.keys():
                    exon = exons[(chromosome, start, stop)]
                else:
                    exon = Exon()
                    exon.start = start
                    exon.stop = stop
                exon.transcripts.append(transcripts[attributes['transcript_id']])
                db_session.add(exon)
                exons[(chromosome, start, stop)] = exon
    db_session.commit() #Attempt to commit all the records
    


import time
start = time.time()
load_data_from_gtf('/Users/bj8th/Documents/Sheynkman-Lab/Data/test/gencode.v35.annotation.chr22.gtf')
end = time.time()
print(f"Time to load gtf file\t{end - start}")
# %%
# statement = select(Transcript)
# result = db_session.execute(statement).all()
# transcript = result.pop()[0]

# print(transcript.name)
# print(transcript.exons)
# print(transcript.gene)
# print(transcript.gene.chromosome.name)
# print(transcript.gene.transcripts)
# print(transcript.length)
# # print(transcript.gene.chromosome.genes)
# # %%
# exon = transcript.exons[0]
# exon.length
# exon.gene
# # %%
