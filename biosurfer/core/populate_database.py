#%%
import logging
import time
from dataclasses import dataclass
from operator import itemgetter
from typing import Any, Callable, Iterable, Tuple, Type

from Bio import SeqIO
from biosurfer.core.constants import APPRIS, Strand
from biosurfer.core.database import Base, db_session, engine
from biosurfer.core.models import (ORF, Chromosome, Exon, GencodeExon,
                                   GencodeTranscript, Gene, PacBioExon, PacBioTranscript, Protein,
                                   Transcript)


@dataclass
class FastaHeaderFields:
    transcript_id: str = None
    protein_id: str = None


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
        tags: list

    """
    chromosome, source, feature, start, stop, score, strand, phase, attributes = line.split('\t')
    start = int(start)
    stop = int(stop)
    attributes = attributes.split(';')[:-1]
    attributes = [att.strip(' ').split(' ') for att in attributes]
    tags = [att[1].strip('"') for att in attributes if att[0] == 'tag']
    attributes = {att[0]: att[1].strip('"') for att in attributes if att[0] != 'tag'}
    return chromosome, source, feature, start, stop, score, strand, phase, attributes, tags

def load_gencode_gtf(gtf_file: str, overwrite=False) -> None:
    if overwrite:
        Base.metadata.drop_all(engine)
        print(f'dropped all tables from {engine.url}')
    Base.metadata.create_all(engine)

    existing_genes = {row.accession for row in db_session.query(Gene.accession).all()}
    existing_transcripts = {row.accession for row in db_session.query(Transcript.accession).all()}
    existing_exons = {(row.accession, row.transcript_id) for row in db_session.query(Exon.accession, Exon.transcript_id).all()}
    existing_orfs = {(row.transcript_id, row.transcript_start, row.transcript_stop) for row in db_session.query(ORF.transcript_id, ORF.transcript_start, ORF.transcript_stop).all()}

    chromosomes = {}
    genes_to_update, genes_to_insert = [], []
    transcripts_to_update, transcripts_to_insert = [], []
    exons_to_update, exons_to_insert = [], []
    orfs_to_update, orfs_to_insert = [], []

    minus_transcripts = set()
    transcripts_to_exons = {}
    transcripts_to_cdss = {}

    def update_and_insert(mapper, mappings_to_update, mappings_to_insert):
        db_session.bulk_update_mappings(mapper, mappings_to_update)
        db_session.bulk_insert_mappings(mapper, mappings_to_insert)
        mappings_to_update[:] = []
        mappings_to_insert[:] = []

    with open(gtf_file) as gtf:
        for i, line in enumerate(gtf, start=1):
            if line.startswith("#"): 
                continue
            chr, source, feature, start, stop, score, strand, phase, attributes, tags = read_gtf_line(line)
            gene_id = attributes['gene_id']
            gene_name = attributes['gene_name']
            
            if feature == 'gene':
                if chr not in chromosomes:
                    chromosome = db_session.merge(Chromosome(name=chr))
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
            
            if i % 20000 == 0:
                print(f'read {i//1000}k lines...')
                update_and_insert(Gene, genes_to_update, genes_to_insert)
                update_and_insert(GencodeTranscript, transcripts_to_update, transcripts_to_insert)
        update_and_insert(Gene, genes_to_update, genes_to_insert)
        update_and_insert(GencodeTranscript, transcripts_to_update, transcripts_to_insert)
    db_session.commit()

    print('calculating transcript-relative exon coordinates...')
    # calculate the coordinates of each exon relative to the sequence of its parent transcript
    for transcript_id, exon_list in transcripts_to_exons.items():
        exon_list.sort(key=itemgetter('start'), reverse=transcript_id in minus_transcripts)
        tx_idx = 1
        for i, exon in enumerate(exon_list, start=1):
            exon['position'] = i
            exon_length = exon['stop'] - exon['start'] + 1
            exon['transcript_start'] = tx_idx
            exon['transcript_stop'] = tx_idx + exon_length - 1
            tx_idx += exon_length
    db_session.bulk_update_mappings(GencodeExon, exons_to_update)
    db_session.bulk_insert_mappings(GencodeExon, exons_to_insert)
    db_session.commit()

    print('calculating transcript-relative ORF coordinates...')
    # assemble CDS intervals into ORFs
    for transcript_id, cds_list in transcripts_to_cdss.items():
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
    db_session.bulk_update_mappings(ORF, orfs_to_update)
    db_session.bulk_insert_mappings(ORF, orfs_to_insert)
    db_session.commit()

def load_pacbio_gtf(gtf_file: str, overwrite=False) -> None:
    if overwrite:
        Base.metadata.drop_all(engine)
        print(f'dropped all tables from {engine.url}')
    Base.metadata.create_all(engine)

    existing_genes = {row.name: row.accession for row in db_session.query(Gene.name, Gene.accession).all()}
    existing_transcripts = {row.accession for row in db_session.query(Transcript.accession).all()}
    existing_exons = {(row.transcript_id, row.position) for row in db_session.query(Exon.transcript_id, Exon.position).all()}
    existing_orfs = {(row.transcript_id, row.transcript_start, row.transcript_stop) for row in db_session.query(ORF.transcript_id, ORF.transcript_start, ORF.transcript_stop).all()}

    chromosomes = {}
    genes_to_update, genes_to_insert = [], []
    transcripts_to_update, transcripts_to_insert = [], []
    exons_to_update, exons_to_insert = [], []
    orfs_to_update, orfs_to_insert = [], []

    minus_transcripts = set()
    transcripts_to_exons = {}
    transcripts_to_cdss = {}

    def update_and_insert(mapper, mappings_to_update, mappings_to_insert):
        db_session.bulk_update_mappings(mapper, mappings_to_update)
        db_session.bulk_insert_mappings(mapper, mappings_to_insert)
        mappings_to_update[:] = []
        mappings_to_insert[:] = []

    with open(gtf_file) as gtf:
        for i, line in enumerate(gtf, start=1):
            if line.startswith("#"): 
                continue
            chr, source, feature, start, stop, score, strand, phase, attributes, tags = read_gtf_line(line)
            if chr not in chromosomes:
                chromosome = db_session.merge(Chromosome(name=chr))
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
            
            if i % 20000 == 0:
                print(f'read {i//1000}k lines...')
                update_and_insert(Gene, genes_to_update, genes_to_insert)
                update_and_insert(PacBioTranscript, transcripts_to_update, transcripts_to_insert)
        update_and_insert(Gene, genes_to_update, genes_to_insert)
        update_and_insert(PacBioTranscript, transcripts_to_update, transcripts_to_insert)
    db_session.commit()

    print('calculating transcript-relative exon coordinates...')
    # calculate the coordinates of each exon relative to the sequence of its parent transcript
    for transcript_id, exon_list in transcripts_to_exons.items():
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
    db_session.bulk_update_mappings(PacBioExon, exons_to_update)
    db_session.bulk_insert_mappings(PacBioExon, exons_to_insert)
    db_session.commit()

    print('calculating transcript-relative ORF coordinates...')
    # assemble CDS intervals into ORFs
    for transcript_id, cds_list in transcripts_to_cdss.items():
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
    db_session.bulk_update_mappings(ORF, orfs_to_update)
    db_session.bulk_insert_mappings(ORF, orfs_to_insert)
    db_session.commit()

def load_transcript_fasta(transcript_fasta: str, id_extractor: Callable[[str], FastaHeaderFields], id_filter: Callable[[str], bool] = lambda x: False):
    existing_transcripts = {row.accession for row in db_session.query(Transcript.accession)}

    transcripts_to_update = []
    for i, record in enumerate(SeqIO.parse(transcript_fasta, 'fasta'), start=1):
        if id_filter(record.id):
            continue
        ids = id_extractor(record.id)
        transcript_id = ids.transcript_id
        sequence = str(record.seq)

        if transcript_id in existing_transcripts:
            transcript = {'accession': transcript_id, 'sequence': sequence}
            transcripts_to_update.append(transcript)
        if i % 20000 == 0:
            print(f'read {i//1000}k entries...')
            db_session.bulk_update_mappings(Transcript, transcripts_to_update)
            transcripts_to_update[:] = []
    db_session.bulk_update_mappings(Transcript, transcripts_to_update)
    db_session.commit() #Attempt to commit all the records

def load_translation_fasta(translation_fasta: str, id_extractor: Callable[[str], FastaHeaderFields], id_filter: Callable[[str], bool] = lambda x: False):
    existing_proteins = {row.accession for row in db_session.query(Protein.accession)}
    existing_orfs = {}
    for transcript_id, start, stop in db_session.query(ORF.transcript_id, ORF.transcript_start, ORF.transcript_stop):
        if transcript_id in existing_orfs:
            existing_orfs[transcript_id].append((start, stop))
        else:
            existing_orfs[transcript_id] = [(start, stop)]

    orfs_to_update = []
    proteins_to_update = []
    proteins_to_insert = []
    for i, record in enumerate(SeqIO.parse(translation_fasta, 'fasta'), start=1):
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
        if i % 20000 == 0:
            print(f'read {i//1000}k entries...')
            db_session.bulk_update_mappings(ORF, orfs_to_update)
            db_session.bulk_update_mappings(Protein, proteins_to_update)
            db_session.bulk_insert_mappings(Protein, proteins_to_insert)
            orfs_to_update[:] = []
            proteins_to_update[:] = []
            proteins_to_insert[:] = []
    db_session.bulk_update_mappings(ORF, orfs_to_update)
    db_session.bulk_update_mappings(Protein, proteins_to_update)
    db_session.bulk_insert_mappings(Protein, proteins_to_insert)     
    db_session.commit()

path = '/home/redox/sheynkman-lab/biosurfer/data/'
gencode_gtf = 'biosurfer_demo_data/gencode.v38.basic.annotation.gtf'
gencode_tx = 'biosurfer_demo_data/gencode.v38.pc_transcripts.fa'
gencode_tl = 'biosurfer_demo_data/gencode.v38.pc_translations.fa'
pacbio_gtf = 'bone/bone_cds_high_confidence.gtf'
pacbio_tx = 'bone/filtered_bone_corrected.fasta'
pacbio_tl = 'bone/bone_hybrid.fasta'

def get_ids_from_gencode_fasta(header: str):
    fields = [field for field in header.split('|') if field and not field.startswith(('UTR', 'CDS'))]
    transcript_id = next((field for field in fields if field.startswith('ENST')))
    protein_id = next((field for field in fields if field.startswith('ENSP')), None)
    return FastaHeaderFields(transcript_id, protein_id)

def get_ids_from_pacbio_fasta(header: str):
    return FastaHeaderFields(header, None)

def get_ids_from_hybrid_fasta(header: str):
    fields = header.split('|')
    return FastaHeaderFields(fields[1], fields[1] + ':PROT1')

def skip_par_y(header: str):  # these have duplicate ENSEMBL accessions and that makes SQLAlchemy very sad
    return 'PAR_Y' in header

def skip_gencode(header: str):
    return header.startswith('gc')

#%%
start = time.time()
load_gencode_gtf(path + gencode_gtf, overwrite=True)
end = time.time()
print(f"Time to load gtf file\t{end - start:0.3g}s")

start = time.time()
load_transcript_fasta(path + gencode_tx, get_ids_from_gencode_fasta, skip_par_y)
end = time.time()
print(f"time to load transcript fasta\t{end - start:0.3g}s")

start = time.time()
load_translation_fasta(path + gencode_tl, get_ids_from_gencode_fasta, skip_par_y)
end = time.time()
print(f"time to load translation fasta\t{end - start:0.3g}s")

#%%
start = time.time()
load_pacbio_gtf(path + pacbio_gtf, overwrite=False)
end = time.time()
print(f"Time to load gtf file\t{end - start:0.3g}s")

start = time.time()
load_transcript_fasta(path + pacbio_tx, get_ids_from_pacbio_fasta)
end = time.time()
print(f"time to load transcript fasta\t{end - start:0.3g}s")

start = time.time()
load_translation_fasta(path + pacbio_tl, get_ids_from_hybrid_fasta, skip_gencode)
end = time.time()
print(f"time to load translation fasta\t{end - start:0.3g}s")
# %%
