import csv
import os
from operator import itemgetter
from typing import Callable

from Bio import SeqIO
from biosurfer.core.constants import APPRIS, SQANTI, STOP_CODONS, AminoAcid, Strand
from biosurfer.core.helpers import (FastaHeaderFields, bulk_update_and_insert,
                                    read_gtf_line)
from biosurfer.core.models import (ORF, Base, Chromosome, Exon,
                                   GencodeTranscript, Gene, PacBioTranscript,
                                   Protein, ProteinFeature, Transcript)
from more_itertools import chunked
from sqlalchemy import create_engine, select
from sqlalchemy.orm import scoped_session, sessionmaker
from tqdm import tqdm

pkg_dir = os.path.dirname(os.path.abspath(__file__))
DB_GENCODE = 'sqlite:///' + os.path.join(pkg_dir, 'gencode.sqlite3')
DB_BONE = 'sqlite:///' + os.path.join(pkg_dir, 'bone.sqlite3')
DB_MEMORY = 'sqlite://'

CHUNK_SIZE = 5000
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
                existing_exons = {(row.transcript_id, row.position) for row in session.query(Exon.transcript_id, Exon.position).all()}
                existing_orfs = {(row.transcript_id, row.position) for row in session.query(ORF.transcript_id, ORF.position).all()}

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
                    if any('_NF' in tag for tag in tags):  # biosurfer does not handle start_NF and end_NF transcripts very well
                        continue
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
                            start_nf, end_nf = False, False
                            # start_nf = start_nf or 'start_NF' in tag
                            # end_nf = end_nf or 'end_NF' in tag
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
            _process_and_upsert_exons(session, transcripts_to_exons, minus_transcripts, existing_exons, exons_to_update, exons_to_insert)

            # assemble CDS intervals into ORFs
            _process_and_upsert_orfs(session, transcripts_to_cdss, transcripts_to_exons, minus_transcripts, existing_orfs, orfs_to_update, orfs_to_insert)

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
                existing_orfs = {(row.transcript_id, row.position) for row in session.query(ORF.transcript_id, ORF.position).all()}

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
            _process_and_upsert_exons(session, transcripts_to_exons, minus_transcripts, existing_exons, exons_to_update, exons_to_insert)

            # assemble CDS intervals into ORFs
            _process_and_upsert_orfs(session, transcripts_to_cdss, transcripts_to_exons, minus_transcripts, existing_orfs, orfs_to_update, orfs_to_insert)

    def load_transcript_fasta(self, transcript_fasta: str, id_extractor: Callable[[str], 'FastaHeaderFields'], id_filter: Callable[[str], bool] = lambda x: False):
        with self.get_session() as session:
            with session.begin():
                existing_transcripts = {row.accession for row in session.query(Transcript.accession)}
                existing_orfs = {
                    row.transcript_id: (row.position, row.transcript_start, row.transcript_stop)
                    for row in session.execute(select(ORF.position, ORF.transcript_id, ORF.transcript_start, ORF.transcript_stop))
                }
            transcripts_to_update = []
            orfs_to_update = []
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
                if transcript_id in existing_orfs:
                    position, tx_start, tx_stop = existing_orfs[transcript_id]
                    # if orf is followed by a stop codon in transcript sequence, modify tx_stop to include the stop codon
                    if sequence[tx_stop:tx_stop+3] in STOP_CODONS:
                        orfs_to_update.append({
                            'transcript_id': transcript_id,
                            'position': position,
                            'transcript_start': tx_start,
                            'transcript_stop': tx_stop + 3,
                            'has_stop_codon': True
                        })

                if len(transcripts_to_update) == CHUNK_SIZE:
                    bulk_update_and_insert(session, Transcript, transcripts_to_update, [])
                    bulk_update_and_insert(session, ORF, orfs_to_update, [])
            bulk_update_and_insert(session, Transcript, transcripts_to_update, [])
            bulk_update_and_insert(session, ORF, orfs_to_update, [])

    def load_translation_fasta(self, translation_fasta: str, id_extractor: Callable[[str], 'FastaHeaderFields'], id_filter: Callable[[str], bool] = lambda x: False):
        with self.get_session() as session:
            with session.begin():
                existing_proteins = {row.accession for row in session.query(Protein.accession)}
                existing_orfs = {}
                for transcript_id, position, start, stop, has_stop_codon in session.query(ORF.transcript_id, ORF.position, ORF.transcript_start, ORF.transcript_stop, ORF.has_stop_codon):
                    if transcript_id in existing_orfs:
                        existing_orfs[transcript_id].append((position, start, stop, has_stop_codon))
                    else:
                        existing_orfs[transcript_id] = [(position, start, stop, has_stop_codon)]

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
                    for position, start, stop, has_stop_codon in orfs:
                        if stop - start + 1 == (seq_length + int(has_stop_codon))*3:
                            orfs_to_update.append({
                                'transcript_id': transcript_id,
                                'position': position,
                                'transcript_start': start,
                                'transcript_stop': stop,
                                'protein_id': protein_id
                            })
                            if has_stop_codon:
                                protein['sequence'] = protein['sequence'] + AminoAcid.STOP.value
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
    
    def load_domain_mappings(self, domain_mapping_file: str, domain_name_file: str):
        with open(domain_name_file) as f:
            reader = csv.DictReader(f, delimiter='\t')
            t = tqdm(reader, desc='Reading domain names', unit='accessions')
            domain_accession_to_name = {row['acc']: row['name'] for row in t}
        with self.get_session() as session:
            with session.begin():
                transcript_name_to_protein_acc = {
                    name: accession
                    for name, accession in session.execute(
                        select(Transcript.name, Protein.accession).
                        join(Protein.orf).
                        join(ORF.transcript)
                    )
                }
            with open(domain_mapping_file) as f:
                reader = csv.DictReader(f, delimiter='\t')
                t = tqdm(reader, desc='Reading domain mappings', unit='mappings')
                domains_to_insert = [
                    {
                        'type': 'feature',
                        'accession': row['Pfam_acc'],
                        'name': domain_accession_to_name.get(row['Pfam_acc'], None),
                        'protein_id': transcript_name_to_protein_acc[row['transcript_name(iso_acc)']],
                        'protein_start': row['start'],
                        'protein_stop': row['end']
                    }
                    for row in t if row['transcript_name(iso_acc)'] in transcript_name_to_protein_acc
                ]
            bulk_update_and_insert(session, ProteinFeature, [], domains_to_insert)        


def _process_and_upsert_exons(session, transcripts_to_exons, minus_transcripts, existing_exons, exons_to_update, exons_to_insert):
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
            if (transcript_id, exon['position']) in existing_exons:
                exons_to_update.append(exon)
            else:
                exons_to_insert.append(exon)
            if 'accession' not in exon:
                exon['accession'] = transcript_id + f':EXON{i}'
            exon_length = exon['stop'] - exon['start'] + 1
            exon['transcript_start'] = tx_idx
            exon['transcript_stop'] = tx_idx + exon_length - 1
            tx_idx += exon_length
    t = tqdm(
        exons_to_update,
        desc = 'Updating existing exons',
        total = len(exons_to_update),
        unit = 'exons'
    )
    for exon_update_chunk in chunked(t, CHUNK_SIZE):
        bulk_update_and_insert(session, Exon, list(exon_update_chunk), [])
    t = tqdm(
        exons_to_insert,
        desc = 'Inserting new exons',
        total = len(exons_to_insert),
        unit = 'exons'
    )
    for exon_insert_chunk in chunked(t, CHUNK_SIZE):
        bulk_update_and_insert(session, Exon, [], list(exon_insert_chunk))


def _process_and_upsert_orfs(session, transcripts_to_cdss, transcripts_to_exons, minus_transcripts, existing_orfs, orfs_to_update, orfs_to_insert):
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
        orf_length = orf_tx_stop - orf_tx_start + 1
        if orf_length % 3 != 0:
            continue  # ORFs with nt lengths indivisible by 3 should not be considered
        orf = {
            'transcript_id': transcript_id,
            'position': 1,
            'transcript_start': orf_tx_start,
            'transcript_stop': orf_tx_stop,
            'has_stop_codon': False,
            # 'start': orf_start,
            # 'stop': orf_stop,
        }
        if (transcript_id, orf['position']) in existing_orfs:
            orfs_to_update.append(orf)
        else:
            orfs_to_insert.append(orf)

        if i % CHUNK_SIZE == 0:
            bulk_update_and_insert(session, ORF, orfs_to_update, orfs_to_insert)
    bulk_update_and_insert(session, ORF, orfs_to_update, orfs_to_insert)
