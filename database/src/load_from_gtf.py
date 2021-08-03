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
    chromosomes = {}
    with open(gtf_file) as gtf:
        for line in gtf:
            if line.startswith("#"): 
                continue
            chr, source, feature, start, stop, score, strand, phase, attributes = read_gtf_line(line)
            if feature == 'gene':
                if chr not in chromosomes.keys():
                  chromosome = Chromosome()
                  chromosome.name = chr
                  chromosomes[chr] = chromosome
                else:
                    chromosome = chromosomes[chr]

                gene = Gene()
                gene.accession = attributes['gene_id']
                gene.name = attributes['gene_name']
                gene.strand = strand
                genes[attributes['gene_id']] = gene
                db_session.add(gene)

                chromosome.genes.append(gene)
                db_session.add(chromosome)
            elif feature == 'transcript':
                transcript = Transcript()
                transcript.accession = attributes['transcript_id']
                transcript.name = attributes['transcript_name']
                genes[attributes['gene_id']].transcripts.append(transcript)
                transcripts[attributes['transcript_id']] = transcript
                db_session.add(transcript)
                
            elif feature == 'exon':
                exon = Exon()
                exon.accession = attributes['exon_id']
                exon.start = start
                exon.stop = stop
                exon.transcript = transcripts[attributes['transcript_id']]
                db_session.add(exon)
    db_session.commit() #Attempt to commit all the records
    

#%%
import time
start = time.time()
load_data_from_gtf('/home/redox/sheynkman-lab/biosurfer/data/biosurfer_demo_data/gencode.v38.annotation.gtf.toy')
end = time.time()
print(f"Time to load gtf file\t{end - start}")


# %%
