#%%
from Bio import SeqIO
from sqlalchemy import select
from database import Base, db_session, engine
from models import Transcript, Exon, Nucleotide
import time
import logging
def load_transcript_fasta(transcript_fasta):
    for record in SeqIO.parse(transcript_fasta, 'fasta'):
        transcript_name = record.id.split('|')[0]
        sequence = str(record.seq)
        try:
            statement = select(Transcript).filter(Transcript.ensembl == transcript_name)
            result = db_session.execute(statement).one()
            transcript = result[0]
            prior_length = 0
            
            for exon in transcript.exons:
                exon_sequence = sequence[prior_length:prior_length + exon.length]
                if exon.sequence != exon_sequence:
                    logging.info(f"updating exon {exon} in transcript {transcript}\nPrior\t{exon.sequence}\nNew\t{exon_sequence}")
                    exon.sequence = exon_sequence
                prior_length = prior_length + exon.length
            db_session.commit() #Attempt to commit all the records
        except:
            logging.warning(f"TRANSCRIPT NOT FOUND IN DATABASE:\t{transcript_name}")
        
    # db_session.commit() #Attempt to commit all the records



start = time.time()
load_transcript_fasta('/home/redox/sheynkman-lab/biosurfer/data/biosurfer_demo_data/gencode.v38.pc_transcripts.fa.toy')
end = time.time()

print(f"time to load fasta\t{end - start} seconds")
print(f"time to load fasta\t{(end - start)/60} minutes")

