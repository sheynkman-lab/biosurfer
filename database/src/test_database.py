#%%
from database import db_session
from models import Gene, Transcript, Exon, ORF, Protein
from sqlalchemy import select

#%%
statement = select(Protein).join(ORF).join(Transcript)
result = db_session.execute(statement).all()
# %%
