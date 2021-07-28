#%%
from database import db_session
from models import Gene, Transcript, Exon, ORF, Protein
from sqlalchemy import select

#%%
statement = select(Protein).join(ORF).join(Transcript)
proteins = [row[0] for row in db_session.execute(statement).all()]
# %%
isoform = db_session.query(Transcript).filter(Transcript.name == 'PAX5-201').one()

# %%
