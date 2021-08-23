#%%
import logging
import time
from operator import attrgetter

from Bio import SeqIO
from inscripta.biocantor.location.location_impl import SingleInterval, Strand
from sqlalchemy import and_, select
from sqlalchemy.exc import NoResultFound
from sqlalchemy.sql import exists

from constants import APPRIS
from models import (ORF, Base, Chromosome, Exon, GencodeExon,
                    GencodeTranscript, Gene, Protein, Transcript, db_session,
                    engine)


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

    """
    chromosome, source, feature, start, stop, score, strand, phase, attributes = line.split('\t')
    start = int(start)
    stop = int(stop)
    attributes = attributes.split(';')[:-1]
    attributes = [att.strip(' ').split(' ') for att in attributes]
    tags = [att[1].strip('"') for att in attributes if att[0] == 'tag']
    attributes = {att[0]: att[1].strip('"') for att in attributes if att[0] != 'tag'}
    return chromosome, source, feature, start, stop, score, strand, phase, attributes, tags

def load_data_from_gtf(gtf_file: str) -> None:
    Base.metadata.create_all(engine)

    existing_genes = {row.accession for row in db_session.query(Gene.accession).all()}
    existing_transcripts = {row.accession for row in db_session.query(Transcript.accession).all()}
    existing_exons = {(row.accession, row.transcript_id) for row in db_session.query(Exon.accession, Exon.transcript_id).all()}

    chromosomes = {}
    genes_to_update, genes_to_insert = [], []
    transcripts_to_update, transcripts_to_insert = [], []
    exons_to_update, exons_to_insert = [], []

    transcripts_to_exons = {}
    minus_transcripts = set()

    def update_and_insert(mapper, mappings_to_update, mappings_to_insert):
        db_session.bulk_update_mappings(mapper, mappings_to_update)
        db_session.bulk_insert_mappings(mapper, mappings_to_insert)
        mappings_to_update[:] = []
        mappings_to_insert[:] = []

    with open(gtf_file) as gtf:
        for i, line in enumerate(gtf):
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
            
            if i % 10000 == 0:
                print(f'read {i//1000}k lines...')
                update_and_insert(Gene, genes_to_update, genes_to_insert)
                update_and_insert(GencodeTranscript, transcripts_to_update, transcripts_to_insert)
        update_and_insert(Gene, genes_to_update, genes_to_insert)
        update_and_insert(GencodeTranscript, transcripts_to_update, transcripts_to_insert)
    
    print('calculating transcript-relative exon coordinates...')
    # calculate the coordinates of each exon relative to the sequence of its parent transcript
    for transcript_id, exon_list in transcripts_to_exons.items():
        strand = Strand.MINUS if transcript_id in minus_transcripts else Strand.PLUS
        exon_to_genomic_loc = [SingleInterval(exon['start'] - 1, exon['stop'], strand) for exon in exon_list]
        transcript_genomic_loc = exon_to_genomic_loc[0]
        for exon_genomic_loc in exon_to_genomic_loc[1:]:
            transcript_genomic_loc = transcript_genomic_loc.union(exon_genomic_loc)
        for i, exon in enumerate(exon_list):
            # TODO: is it faster to just loop through transcript's exons and count off lengths manually?
            exon_transcript_loc = exon_to_genomic_loc[i].location_relative_to(transcript_genomic_loc)
            exon['transcript_start'] = exon_transcript_loc.start + 1
            exon['transcript_stop'] = exon_transcript_loc.end
        # transcript.exons.sort(key=attrgetter('transcript_start'))
    db_session.bulk_update_mappings(GencodeExon, exons_to_update)
    db_session.bulk_insert_mappings(GencodeExon, exons_to_insert)

    db_session.commit()

def load_transcript_fasta(transcript_fasta):
    existing_transcripts = {row.accession for row in db_session.query(Transcript.accession)}
    existing_orfs = set(db_session.query(ORF.transcript_id, ORF.transcript_start, ORF.transcript_stop))

    transcripts_to_update = []
    orfs_to_insert = []
    for i, record in enumerate(SeqIO.parse(transcript_fasta, 'fasta')):
        fields = record.id.split('|')
        transcript_id = fields[0]
        orf_coords = [field[4:] for field in fields if field.startswith('CDS:')][0]
        orf_start = int(orf_coords.split('-')[0])
        orf_stop = int(orf_coords.split('-')[1])
        sequence = str(record.seq)

        if transcript_id in existing_transcripts:
            transcript = {'accession': transcript_id, 'sequence': sequence}
            transcripts_to_update.append(transcript)
            orf = {
                'transcript_id': transcript_id,
                'transcript_start': orf_start,
                'transcript_stop': orf_stop
            }
            if (transcript_id, orf_start, orf_stop) not in existing_orfs:
                orfs_to_insert.append(orf)
        if i % 10000 == 0:
            print(f'read {i//1000}k entries...')
            db_session.bulk_update_mappings(Transcript, transcripts_to_update)
            db_session.bulk_insert_mappings(ORF, orfs_to_insert)
            transcripts_to_update[:] = []
            orfs_to_insert[:] = []
    db_session.bulk_update_mappings(Transcript, transcripts_to_update)
    db_session.bulk_insert_mappings(ORF, orfs_to_insert)
    db_session.commit() #Attempt to commit all the records

def load_translation_fasta(translation_fasta):
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
    for i, record in enumerate(SeqIO.parse(translation_fasta, 'fasta')):
        fields = record.id.split('|')
        protein_id, transcript_id = fields[:2]
        if 'PAR_Y' in transcript_id:
            continue  # these have duplicate ENSEMBL accessions and that makes SQLAlchemy very sad
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
        if i % 10000 == 0:
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

path = '/home/redox/sheynkman-lab/biosurfer/data/biosurfer_demo_data/'
gtf_file = 'gencode.v38.basic.annotation.gtf'
tx_file = 'gencode.v38.pc_transcripts.fa'
tl_file = 'gencode.v38.pc_translations.fa'
# gtf_file = 'gencode.v38.annotation.gtf.toy'
# tx_file = 'gencode.v38.pc_transcripts.fa.toy'
# tl_file = 'gencode.v38.pc_translations.fa.toy'

#%%
start = time.time()
load_data_from_gtf(path + gtf_file)
end = time.time()
print(f"Time to load gtf file\t{end - start:0.3g}s")

#%%
start = time.time()
load_transcript_fasta(path + tx_file)
end = time.time()
print(f"time to load transcript fasta\t{end - start:0.3g}s / {(end - start)/60:0.3g}min")

#%%
start = time.time()
load_translation_fasta(path + tl_file)
end = time.time()
print(f"time to load translation fasta\t{end - start:0.3g}s / {(end - start)/60:0.3g}min")
# %%
