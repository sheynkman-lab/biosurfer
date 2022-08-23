#%%
from time import time

from biosurfer.core.database import Database
from biosurfer.core.helpers import get_ids_from_gencode_fasta, skip_par_y
import os

#%%
def check_database(gencode_gtf, gencode_tx, gencode_tl, gencode_doms, pfam_dom_info, prosite_patterns, db_name):    

    #check if database already exits
    path = os.getcwd() +'/'
    db_path = path + 'databases'
    db_path_list = os.listdir(db_path)
    print(db_path_list)

    if (db_name + '.sqlite3') in db_path_list:
        print('\n Database already exists. Loading ' + db_name + ' ...')
        load_gencode(gencode_gtf, gencode_tx, gencode_tl, gencode_doms, pfam_dom_info, prosite_patterns, db_name, path)
    else:
        print('\n Creating a new database ' + db_name + ' ...')
        load_gencode(gencode_gtf, gencode_tx, gencode_tl, gencode_doms, pfam_dom_info, prosite_patterns, db_name, path)    

#%%
def load_gencode(gencode_gtf, gencode_tx, gencode_tl, gencode_doms, pfam_dom_info, prosite_patterns, db_name, path):
    
    db = Database(db_name)
    gencode_gtf = 'data/gencode/' + os.path.basename(gencode_gtf)
    gencode_tx = 'data/gencode/' + os.path.basename(gencode_tx)
    gencode_tl = 'data/gencode/' + os.path.basename(gencode_tl)
    gencode_doms = 'data/gencode/' + os.path.basename(gencode_doms)
    pfam_dom_info = 'data/gencode/' + os.path.basename(pfam_dom_info)
    prosite_patterns = 'data/gencode/' + os.path.basename(prosite_patterns)

    start = time()
    # db.recreate_tables()
    #TODO: Join paths using pathlib :  os.path.join()
    #TODO: Script doesnt read gtf file to create database
    db.load_gencode_gtf(path + gencode_gtf, overwrite=True)
    db.load_transcript_fasta(path + gencode_tx, get_ids_from_gencode_fasta, skip_par_y)
    db.load_translation_fasta(path + gencode_tl, get_ids_from_gencode_fasta, skip_par_y, overwrite=True)
    db.load_domains(path + pfam_dom_info)
    db.load_patterns(path + prosite_patterns)
    db.load_feature_mappings(path + gencode_doms, overwrite=False)
    end = time()
    print(f'Total time: {end - start:.3g} s')

#%%
if __name__ == '__main__':
    check_database(gencode_gtf, gencode_tx, gencode_tl, gencode_doms, pfam_dom_info, prosite_patterns, db_name)
