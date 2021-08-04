#%%
import logging
import time

from Bio import SeqIO
from sqlalchemy import select
from sqlalchemy.exc import NoResultFound

from database import Base, db_session, engine
from models import (ORF, Chromosome, GencodeExon, GencodeTranscript, Gene,
                    Protein, Transcript)


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
                transcript = GencodeTranscript(
                    accession = attributes['transcript_id'],
                    name = attributes['transcript_name'],
                    # TODO: implement these tags
                    appris = None,
                    start_nf = None,
                    end_nf = None
                )
                # transcript.accession = attributes['transcript_id']
                # transcript.name = attributes['transcript_name']
                genes[attributes['gene_id']].transcripts.append(transcript)
                transcripts[attributes['transcript_id']] = transcript
                db_session.add(transcript)
                
            elif feature == 'exon':
                exon = GencodeExon(
                    accession = attributes['exon_id'],
                    start = start,
                    stop = stop
                )
                exon.accession = attributes['exon_id']
                exon.start = start
                exon.stop = stop
                exon.transcript = transcripts[attributes['transcript_id']]
                db_session.add(exon)
    
    # TODO: calculate exon.start_tx and exon.stop_tx for transcript in transcripts for exon in transcript

    db_session.commit() #Attempt to commit all the records
    
def load_transcript_fasta(transcript_fasta):
    for record in SeqIO.parse(transcript_fasta, 'fasta'):
        fields = record.id.split('|')
        transcript_name = fields[0]
        orf_coords = [field[4:] for field in fields if field.startswith('CDS:')][0]
        orf_start = int(orf_coords.split('-')[0])
        orf_end = int(orf_coords.split('-')[1])
        sequence = str(record.seq)

        try:
            statement = select(Transcript).filter(Transcript.accession == transcript_name)
            result = db_session.execute(statement).one()
            transcript = result[0]

            orf = ORF()
            orf.transcript_start, orf.transcript_stop = orf_start, orf_end
            orf.transcript = transcript
            # TODO: calculate genomic start and end? or get from CDS?
            db_session.add(orf)

            prior_length = 0
            for exon in transcript.exons:
                exon_sequence = sequence[prior_length:prior_length + exon.length]
                if exon.sequence != exon_sequence:
                    logging.info(f"updating exon {exon} in transcript {transcript}\nPrior\t{exon.sequence}\nNew\t{exon_sequence}")
                    exon.sequence = exon_sequence
                prior_length = prior_length + exon.length
            db_session.commit() #Attempt to commit all the records
        except NoResultFound:
            logging.warning(f'could not get transcript {transcript_name} from database')
        
    # db_session.commit() #Attempt to commit all the records

def load_translation_fasta(translation_fasta):
    for record in SeqIO.parse(translation_fasta, 'fasta'):
        fields = record.id.split('|')
        protein_id, transcript_id = fields[:2]
        sequence = str(record.seq)
        
        try:
            statement = select(Transcript).filter(Transcript.accession == transcript_id)
            result = db_session.execute(statement).one()
            transcript = result[0]

            protein = Protein()
            protein.sequence = sequence
            protein.orf = transcript.orfs[0]
            db_session.add(protein)
        except NoResultFound:
            logging.warning(f'could not get transcript {transcript_id} from database')
        
    db_session.commit()

path = '/home/redox/sheynkman-lab/biosurfer/data/biosurfer_demo_data/'
gtf_file = 'gencode.v38.annotation.gtf.toy'
tx_file = 'gencode.v38.pc_transcripts.fa.toy'
tl_file = 'gencode.v38.pc_translations.fa.toy'

#%%
start = time.time()
load_data_from_gtf(path + gtf_file)
end = time.time()
print(f"Time to load gtf file\t{end - start:0.3g}")

#%%
start = time.time()
load_transcript_fasta(path + tx_file)
end = time.time()
print(f"time to load transcript fasta\t{end - start:0.3g} seconds")
print(f"time to load transcript fasta\t{(end - start)/60:0.3g} minutes")

#%%
start = time.time()
load_translation_fasta(path + tl_file)
end = time.time()
print(f"time to load translation fasta\t{end - start:0.3g} seconds")
# %%
