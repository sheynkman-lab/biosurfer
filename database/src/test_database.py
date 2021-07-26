#%%
from database import  db_session
from models import Transcript
from Bio import SeqIO
from sqlalchemy import select
# from models import Transcript
import time
start = time.time()
statement = select(Transcript).where(Transcript.name == 'ENST00000436958.5')
transcript = db_session.execute(statement).one()[0]
stop1 = time.time()
seq = transcript.sequence
stop2 = time.time()
nucs = transcript.nucleotides
stop3 = time.time()

print(f"load transcript\t{stop1-start}")
print(f"get sequence\t{stop2-stop1}")
print(f"get nucloetides\t{stop3-stop2}")
# transcript.sequence
# transcript = result.pop()[0]
#%%
import time
start = time.time()
transcript_fasta = '/Users/bj8th/Documents/Sheynkman-Lab/Data/test/gencode.v35.pc_transcripts.chr22.small.fa'
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
# statement = select(Transcript).where(Transcript.name == 'ENST00000262607.3')
# transcript = db_session.execute(statement).all()[0]
statement = select(Transcript)
result = db_session.execute(statement).all()

transcript = result[11][0]
print(len(transcript.sequence))

# transcript.sequence == transcript_sequences[transcript.name]
stop = time.time()
print(stop-start)

# %%
for i in range(1000):
    transcript = result[i][0]
    if len(transcript.exons) > 1 and len(transcript.sequence) > 0:
         print((transcript))
         break

# %%


