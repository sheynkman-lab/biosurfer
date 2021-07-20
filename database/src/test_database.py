#%%
from database import  db_session
from models import Transcript
from Bio import SeqIO
from sqlalchemy import select
# from models import Transcript

statement = select(Transcript).where(Transcript.name == 'ENST00000436958.5')
transcript = db_session.execute(statement).one()[0]
# transcript.sequence
# transcript = result.pop()[0]
#%%
import time
start = time.time()
transcript_fasta = '/Users/bj8th/Documents/Sheynkman-Lab/Data/test/gencode.v35.pc_transcripts.chr22.fa'
transcript_sequences = {}
for record in SeqIO.parse(transcript_fasta, 'fasta'):
    transcript_name = record.id.split('|')[0]
    sequence = str(record.seq)
    transcript_sequences[transcript_name] = sequence

    statement = select(Transcript).where(Transcript.name == transcript_name)
    transcript = db_session.execute(statement).one()[0]
    if transcript.sequence != sequence:
        print(transcript)
    


# %%
import time
start = time.time()
statement = select(Transcript).where(Transcript.name == 'ENST00000262607.3')
transcript = db_session.execute(statement).one()[0]
# statement = select(Transcript)
# result = db_session.execute(statement).all()
# transcript = result[8][0]

transcript.sequence
# transcript.sequence == transcript_sequences[transcript.name]
stop = time.time()
print(stop-start)
# %%
