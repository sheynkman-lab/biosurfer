#%%
from database import db_session
from models import Gene, Transcript
from Bio import SeqIO
from sqlalchemy import select

# %%
import time
start = time.time()
# statement = select(Transcript).where(Transcript.name == 'ENST00000262607.3')
# transcript = db_session.execute(statement).all()[0]
statement = select(Gene)
result = db_session.execute(statement).all()
transcript = result[11][0]
# transcript.sequence == transcript_sequences[transcript.name]
stop = time.time()
print(f'loaded {len(transcript.sequence)} transcript sequences in {stop-start:0.3g} s')

# %%


