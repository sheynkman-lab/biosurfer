import csv
import os
from operator import itemgetter
from typing import Callable

from Bio import SeqIO
from biosurfer.core.constants import APPRIS, SQANTI, Strand
from biosurfer.core.helpers import (FastaHeaderFields, bulk_update_and_insert,
                                    read_gtf_line)
from biosurfer.core.models import (ORF, Base, Chromosome, Exon, GencodeExon,
                                   GencodeTranscript, Gene, PacBioExon,
                                   PacBioTranscript, Protein, Transcript)
from sqlalchemy import create_engine
from sqlalchemy.orm import scoped_session, sessionmaker
from tqdm import tqdm

pkg_dir = os.path.dirname(os.path.abspath(__file__))
DB_GENCODE = 'sqlite:///' + os.path.join(pkg_dir, 'gencode.sqlite3')
DB_BONE = 'sqlite:///' + os.path.join(pkg_dir, 'bone.sqlite3')
DB_MEMORY = 'sqlite://'

CHUNK_SIZE = 10000
SQANTI_DICT = {
    'full-splice_match': SQANTI.FSM,
    'incomplete-splice_match': SQANTI.ISM,
    'novel_in_catalog': SQANTI.NIC,
    'novel_not_in_catalog': SQANTI.NNC
}


class Database:
    registry = {}

    def __new__(cls, path):
        if path in Database.registry:
            return Database.registry[path]
        else:
            obj = super().__new__(cls)
            Database.registry[path] = obj
            return obj

    def __init__(self, path, sessionfactory=None):
        self.path = path
        self._engine = None
        if sessionfactory is None:
            self._sessionmaker = scoped_session(sessionmaker(autocommit=False, autoflush=False, bind=self.engine))
        else:
            self._sessionmaker = sessionfactory
    
    @property
    def engine(self):
        if self._engine is None:
            self._engine = create_engine(self.path)
        return self._engine

    def get_session(self, **kwargs):
        return self._sessionmaker(**kwargs)

    def load_gencode_gtf(self, gtf_file: str, overwrite=False) -> None:
        if overwrite:
            Base.metadata.drop_all(self.engine)
            print(f'Dropped all tables from {self.engine.url}')
        Base.metadata.create_all(self.engine)

        with self.get_session() as session:
            with session.begin():
                existing_genes = {row.accession for row in session.query(Gene.accession).all()}
                existing_transcripts = {row.accession for row in session.query(Transcript.accession).all()}
                existing_exons = {(row.accession, row.transcript_id) for row in session.query(Exon.accession, Exon.transcript_id).all()}
                existing_orfs = {(row.transcript_id, row.transcript_start, row.transcript_stop) for row in session.query(ORF.transcript_id, ORF.transcript_start, ORF.transcript_stop).all()}

            chromosomes = {}
            genes_to_update, genes_to_insert = [], []
            transcripts_to_update, transcripts_to_insert = [], []
            exons_to_update, exons_to_insert = [], []
            orfs_to_update, orfs_to_insert = [], []

            minus_transcripts = set()
            transcripts_to_exons = {}
            transcripts_to_cdss = {}

            with open(gtf_file) as gtf:
                t = tqdm(gtf, desc=f'Reading {gtf_file} as GENCODE', unit='lines')
                for i, line in enumerate(t, start=1):
                    if line.startswith("#"): 
                        continue
                    chr, source, feature, start, stop, score, strand, phase, attributes, tags = read_gtf_line(line)
                    gene_id = attributes['gene_id']
                    gene_name = attributes['gene_name']
                    if attributes['gene_type'] != 'protein_coding':
                        continue
                    
                    if feature == 'gene':
                        if chr not in chromosomes:
                            with session.begin():
                                chromosome = session.merge(Chromosome(name=chr))
                            chromosomes[chr] = chromosome
                        else:
                            chromosome = chromosomes[chr]

                        gene = {
                            'accession': gene_id,
                            'name': gene_name,
                            'strand': Strand.from_symbol(strand),
                            'chromosome_id': chr
                        }
                        if gene_id in existing_genes:
                            genes_to_update.append(gene)
                        else:
                            genes_to_insert.append(gene)
                            
                    elif feature == 'transcript':
                        transcript_id = attributes['transcript_id']
                        transcript_name = attributes['transcript_name']
                        appris = APPRIS.NONE
                        start_nf, end_nf = False, False
                        for tag in tags:
                            if 'appris_principal' in tag:
                                appris = APPRIS.PRINCIPAL
                            if 'appris_alternative' in tag:
                                appris = APPRIS.ALTERNATIVE
                            start_nf = start_nf or 'start_NF' in tag
                            end_nf = end_nf or 'end_NF' in tag
                        transcript = {
                            'accession': transcript_id,
                            'name': transcript_name,
                            'type': 'gencodetranscript',
                            'gene_id': gene_id,
                            'strand': Strand.from_symbol(strand),
                            'appris': appris,
                            'start_nf': start_nf,
                            'end_nf': end_nf
                        }
                        if transcript_id in existing_transcripts:
                            transcripts_to_update.append(transcript)
                        else:
                            transcripts_to_insert.append(transcript)
                        if Strand.from_symbol(strand) is Strand.MINUS:
                            minus_transcripts.add(transcript_id)
                        
                    elif feature == 'exon':
                        transcript_id = attributes['transcript_id']
                        exon_id = attributes['exon_id']
                        exon = {
                            'accession': exon_id,
                            'type': 'gencodeexon',
                            'start': start,
                            'stop': stop,
                            'transcript_id': transcript_id
                        }
                        if (exon_id, transcript_id) in existing_exons:
                            exons_to_update.append(exon)
                        else:
                            exons_to_insert.append(exon)
                        if transcript_id in transcripts_to_exons:
                            transcripts_to_exons[transcript_id].append(exon)
                        else:
                            transcripts_to_exons[transcript_id] = [exon,]
                    
                    elif feature == 'CDS':
                        transcript_id = attributes['transcript_id']
                        protein_id = attributes['protein_id']
                        cds = (start, stop, protein_id)
                        if transcript_id in transcripts_to_cdss:
                            transcripts_to_cdss[transcript_id].append(cds)
                        else:
                            transcripts_to_cdss[transcript_id] = [cds,]
                    
                    if i % CHUNK_SIZE == 0:
                        bulk_update_and_insert(session, Gene, genes_to_update, genes_to_insert)
                        bulk_update_and_insert(session, GencodeTranscript, transcripts_to_update, transcripts_to_insert)
                bulk_update_and_insert(session, Gene, genes_to_update, genes_to_insert)
                bulk_update_and_insert(session, GencodeTranscript, transcripts_to_update, transcripts_to_insert)

            # calculate the coordinates of each exon relative to the sequence of its parent transcript
            t = tqdm(
                transcripts_to_exons.items(),
                desc = 'Calculating transcript-relative exon coords',
                total = len(transcripts_to_exons),
                unit = 'transcripts'
            )
            for i, (transcript_id, exon_list) in enumerate(t):
                exon_list.sort(key=itemgetter('start'), reverse=transcript_id in minus_transcripts)
                tx_idx = 1
                for i, exon in enumerate(exon_list, start=1):
                    exon['position'] = i
                    exon_length = exon['stop'] - exon['start'] + 1
                    exon['transcript_start'] = tx_idx
                    exon['transcript_stop'] = tx_idx + exon_length - 1
                    tx_idx += exon_length

                if i % CHUNK_SIZE == 0:
                    bulk_update_and_insert(session, GencodeExon, exons_to_update, exons_to_insert)
            bulk_update_and_insert(session, GencodeExon, exons_to_update, exons_to_insert)

            # assemble CDS intervals into ORFs
            t = tqdm(
                transcripts_to_cdss.items(),
                desc = 'Calculating transcript-relative ORF coords',
                total = len(transcripts_to_cdss),
                unit = 'transcripts'
            )
            for i, (transcript_id, cds_list) in enumerate(t):
                assert len({cds[2] for cds in cds_list}) == 1
                exon_list = transcripts_to_exons[transcript_id]
                cds_list.sort(key=itemgetter(0), reverse=transcript_id in minus_transcripts)
                first_cds, last_cds = cds_list[0], cds_list[-1]
                # assuming that the first and last CDSs are the ORF boundaries -- won't work when dealing with multiple ORFs
                orf_start = first_cds[0]
                orf_stop = last_cds[1]
                if transcript_id in minus_transcripts:
                    orf_start, orf_stop = orf_stop, orf_start
                first_exon = next((exon for exon in exon_list if exon['start'] <= first_cds[0] and first_cds[1] <= exon['stop']), None)
                last_exon = next((exon for exon in reversed(exon_list) if exon['start'] <= last_cds[0] and last_cds[1] <= exon['stop']), None)
                # find ORF start/end relative to exons
                if transcript_id not in minus_transcripts:
                    first_offset = first_cds[0] - first_exon['start']
                    last_offset = last_exon['stop'] - last_cds[1]
                else:
                    first_offset = first_exon['stop'] - first_cds[1]
                    last_offset = last_cds[0] - last_exon['start']
                # convert to transcript-relative coords
                orf_tx_start = first_exon['transcript_start'] + first_offset
                orf_tx_stop = last_exon['transcript_stop'] - last_offset
                orf_tx_stop += 3  # account for stop codon
                orf = {
                    'transcript_id': transcript_id,
                    'transcript_start': orf_tx_start,
                    'transcript_stop': orf_tx_stop,
                    'start': orf_start,
                    'stop': orf_stop
                }
                if (transcript_id, orf_tx_start, orf_tx_stop) in existing_orfs:
                    orfs_to_update.append(orf)
                else:
                    orfs_to_insert.append(orf)
                
                if i % CHUNK_SIZE == 0:
                    bulk_update_and_insert(session, ORF, orfs_to_update, orfs_to_insert)
            bulk_update_and_insert(session, ORF, orfs_to_update, orfs_to_insert)

    def load_pacbio_gtf(self, gtf_file: str, overwrite=False) -> None:
        if overwrite:
            Base.metadata.drop_all(self.engine)
            print(f'Dropped all tables from {self.engine.url}')
        Base.metadata.create_all(self.engine)

        with self.get_session() as session:
            with session.begin():
                existing_genes = {row.name: row.accession for row in session.query(Gene.name, Gene.accession).all()}
                existing_transcripts = {row.accession for row in session.query(Transcript.accession).all()}
                existing_exons = {(row.transcript_id, row.position) for row in session.query(Exon.transcript_id, Exon.position).all()}
                existing_orfs = {(row.transcript_id, row.transcript_start, row.transcript_stop) for row in session.query(ORF.transcript_id, ORF.transcript_start, ORF.transcript_stop).all()}

            chromosomes = {}
            genes_to_update, genes_to_insert = [], []
            transcripts_to_update, transcripts_to_insert = [], []
            exons_to_update, exons_to_insert = [], []
            orfs_to_update, orfs_to_insert = [], []

            minus_transcripts = set()
            transcripts_to_exons = {}
            transcripts_to_cdss = {}

            with open(gtf_file) as gtf:
                t = tqdm(gtf, desc=f'Reading {gtf_file} as PacBio', unit='lines')
                for i, line in enumerate(t, start=1):
                    if line.startswith("#"): 
                        continue
                    chr, source, feature, start, stop, score, strand, phase, attributes, tags = read_gtf_line(line)
                    if chr not in chromosomes:
                        with session.begin():
                            chromosome = session.merge(Chromosome(name=chr))
                        chromosomes[chr] = chromosome
                    gene_name = attributes['gene_id']
                    if gene_name not in existing_genes:
                        existing_genes[gene_name] = gene_name
                        gene = {
                            'accession': gene_name,
                            'name': gene_name,
                            'strand': Strand.from_symbol(strand),
                            'chromosome_id': chr
                        }
                        genes_to_insert.append(gene)
                    gene_id = existing_genes[gene_name]

                    if feature == 'transcript':
                        transcript_id = attributes['transcript_id'].split('|')[1]
                        transcript_name = gene_name + '|' + transcript_id
                        transcript = {
                            'accession': transcript_id,
                            'name': transcript_name,
                            'type': 'pacbiotranscript',
                            'gene_id': gene_id,
                            'strand': Strand.from_symbol(strand),
                        }
                        if transcript_id in existing_transcripts:
                            transcripts_to_update.append(transcript)
                        else:
                            transcripts_to_insert.append(transcript)
                        if Strand.from_symbol(strand) is Strand.MINUS:
                            minus_transcripts.add(transcript_id)
                        
                    elif feature == 'exon':
                        transcript_id = attributes['transcript_id'].split('|')[1]
                        exon = {
                            'type': 'pacbioexon',
                            'start': start,
                            'stop': stop,
                            'transcript_id': transcript_id
                        }
                        if transcript_id in transcripts_to_exons:
                            transcripts_to_exons[transcript_id].append(exon)
                        else:
                            transcripts_to_exons[transcript_id] = [exon,]
                    
                    elif feature == 'CDS':
                        transcript_id = attributes['transcript_id'].split('|')[1]
                        cds = (start, stop)
                        if transcript_id in transcripts_to_cdss:
                            transcripts_to_cdss[transcript_id].append(cds)
                        else:
                            transcripts_to_cdss[transcript_id] = [cds,]
                    
                    if i % CHUNK_SIZE == 0:
                        bulk_update_and_insert(session, Gene, genes_to_update, genes_to_insert)
                        bulk_update_and_insert(session, PacBioTranscript, transcripts_to_update, transcripts_to_insert)
                bulk_update_and_insert(session, Gene, genes_to_update, genes_to_insert)
                bulk_update_and_insert(session, PacBioTranscript, transcripts_to_update, transcripts_to_insert)

            # calculate the coordinates of each exon relative to the sequence of its parent transcript
            t = tqdm(
                transcripts_to_exons.items(),
                desc = 'Calculating transcript-relative exon coords',
                total = len(transcripts_to_exons),
                unit = 'transcripts'
            )
            for i, (transcript_id, exon_list) in enumerate(t):
                exon_list.sort(key=itemgetter('start'), reverse=transcript_id in minus_transcripts)
                tx_idx = 1
                for i, exon in enumerate(exon_list, start=1):
                    exon['position'] = i
                    exon['accession'] = transcript_id + f':EXON{i}'
                    if (transcript_id, exon['position']) in existing_exons:
                        exons_to_update.append(exon)
                    else:
                        exons_to_insert.append(exon)
                    exon_length = exon['stop'] - exon['start'] + 1
                    exon['transcript_start'] = tx_idx
                    exon['transcript_stop'] = tx_idx + exon_length - 1
                    tx_idx += exon_length

                if i % CHUNK_SIZE == 0:
                    bulk_update_and_insert(session, PacBioExon, exons_to_update, exons_to_insert)
            bulk_update_and_insert(session, PacBioExon, exons_to_update, exons_to_insert)

            # assemble CDS intervals into ORFs
            t = tqdm(
                transcripts_to_cdss.items(),
                desc = 'Calculating transcript-relative ORF coords',
                total = len(transcripts_to_cdss),
                unit = 'transcripts'
            )
            for i, (transcript_id, cds_list) in enumerate(t):
                exon_list = transcripts_to_exons[transcript_id]
                cds_list.sort(key=itemgetter(0), reverse=transcript_id in minus_transcripts)
                first_cds, last_cds = cds_list[0], cds_list[-1]
                # assuming that the first and last CDSs are the ORF boundaries -- won't work when dealing with multiple ORFs
                orf_start = first_cds[0]
                orf_stop = last_cds[1]
                if transcript_id in minus_transcripts:
                    orf_start, orf_stop = orf_stop, orf_start
                first_exon = next((exon for exon in exon_list if exon['start'] <= first_cds[0] and first_cds[1] <= exon['stop']), None)
                last_exon = next((exon for exon in reversed(exon_list) if exon['start'] <= last_cds[0] and last_cds[1] <= exon['stop']), None)
                # find ORF start/end relative to exons
                if transcript_id not in minus_transcripts:
                    first_offset = first_cds[0] - first_exon['start']
                    last_offset = last_exon['stop'] - last_cds[1]
                else:
                    first_offset = first_exon['stop'] - first_cds[1]
                    last_offset = last_cds[0] - last_exon['start']
                # convert to transcript-relative coords
                orf_tx_start = first_exon['transcript_start'] + first_offset
                orf_tx_stop = last_exon['transcript_stop'] - last_offset
                orf_tx_stop += 3  # account for stop codon
                orf = {
                    'transcript_id': transcript_id,
                    'transcript_start': orf_tx_start,
                    'transcript_stop': orf_tx_stop,
                    'start': orf_start,
                    'stop': orf_stop
                }
                if (transcript_id, orf_tx_start, orf_tx_stop) in existing_orfs:
                    orfs_to_update.append(orf)
                else:
                    orfs_to_insert.append(orf)

                if i % CHUNK_SIZE == 0:
                    bulk_update_and_insert(session, ORF, orfs_to_update, orfs_to_insert)
            bulk_update_and_insert(session, ORF, orfs_to_update, orfs_to_insert)

    def load_transcript_fasta(self, transcript_fasta: str, id_extractor: Callable[[str], 'FastaHeaderFields'], id_filter: Callable[[str], bool] = lambda x: False):
        with self.get_session() as session:
            with session.begin():
                existing_transcripts = {row.accession for row in session.query(Transcript.accession)}
            transcripts_to_update = []
            t = tqdm(SeqIO.parse(transcript_fasta, 'fasta'), desc='Reading transcripts fasta', unit='seq')
            for record in t:
                if id_filter(record.id):
                    continue
                ids = id_extractor(record.id)
                transcript_id = ids.transcript_id
                sequence = str(record.seq)
                if transcript_id in existing_transcripts:
                    transcript = {'accession': transcript_id, 'sequence': sequence}
                    transcripts_to_update.append(transcript)
                if len(transcripts_to_update) == CHUNK_SIZE:
                    bulk_update_and_insert(session, Transcript, transcripts_to_update, [])
            bulk_update_and_insert(session, Transcript, transcripts_to_update, [])

    def load_translation_fasta(self, translation_fasta: str, id_extractor: Callable[[str], 'FastaHeaderFields'], id_filter: Callable[[str], bool] = lambda x: False):
        with self.get_session() as session:
            with session.begin():
                existing_proteins = {row.accession for row in session.query(Protein.accession)}
                existing_orfs = {}
                for transcript_id, start, stop in session.query(ORF.transcript_id, ORF.transcript_start, ORF.transcript_stop):
                    if transcript_id in existing_orfs:
                        existing_orfs[transcript_id].append((start, stop))
                    else:
                        existing_orfs[transcript_id] = [(start, stop)]

            orfs_to_update = []
            proteins_to_update = []
            proteins_to_insert = []
            t = tqdm(SeqIO.parse(translation_fasta, 'fasta'), desc='Reading translations fasta', unit='seqs')
            for i, record in enumerate(t, start=1):
                if id_filter(record.id):
                    continue
                ids = id_extractor(record.id)
                transcript_id = ids.transcript_id
                protein_id = ids.protein_id 
                sequence = str(record.seq)
                seq_length = len(sequence)

                protein = {'accession': protein_id, 'sequence': sequence}
                if protein_id in existing_proteins:
                    proteins_to_update.append(protein)
                else:
                    proteins_to_insert.append(protein)
                if transcript_id in existing_orfs:
                    orfs = existing_orfs[transcript_id]
                    for start, stop in orfs:
                        if stop - start + 1 == (seq_length + 1)*3:
                            orfs_to_update.append({
                                'transcript_id': transcript_id,
                                'transcript_start': start,
                                'transcript_stop': stop,
                                'protein_id': protein_id
                            })
                            break
                
                if i % CHUNK_SIZE == 0:
                    bulk_update_and_insert(session, ORF, orfs_to_update, [])
                    bulk_update_and_insert(session, Protein, proteins_to_update, proteins_to_insert)
            bulk_update_and_insert(session, ORF, orfs_to_update, [])
            bulk_update_and_insert(session, Protein, proteins_to_update, proteins_to_insert)  

    def load_sqanti_classifications(self, sqanti_file: str):
        with self.get_session() as session:
            with session.begin():
                existing_gencode_transcripts = {row.accession for row in session.query(GencodeTranscript.accession)}
                existing_pacbio_transcripts = {row.accession for row in session.query(PacBioTranscript.accession)}
            transcripts_to_update = []
            with open(sqanti_file) as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in tqdm(reader, desc='Reading SQANTI classifications', unit='lines'):
                    if row['isoform'] in existing_pacbio_transcripts:
                        tx = {
                            'accession': row['isoform'],
                            'sqanti': SQANTI_DICT.get(row['structural_category'], SQANTI.OTHER)
                        }
                        associated_tx = row['associated_transcript']
                        if tx['sqanti'] in {SQANTI.FSM, SQANTI.ISM} and associated_tx in existing_gencode_transcripts:
                            tx['gencode_id'] = associated_tx
                        transcripts_to_update.append(tx)
                    
                    if len(transcripts_to_update) == CHUNK_SIZE:
                        bulk_update_and_insert(session, PacBioTranscript, transcripts_to_update, [])
                bulk_update_and_insert(session, PacBioTranscript, transcripts_to_update, [])
