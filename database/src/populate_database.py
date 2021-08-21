#%%
import logging
import time
from operator import attrgetter

from Bio import SeqIO
from inscripta.biocantor.location.location_impl import SingleInterval
from sqlalchemy import select
from sqlalchemy.exc import NoResultFound

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
    genes = {}
    transcripts = {}
    chromosomes = {}
    with open(gtf_file) as gtf:
        for i, line in enumerate(gtf):
            if i % 5000 == 0:
                print(f'read {i//1000}k lines...')
            if line.startswith("#"): 
                continue
            chr, source, feature, start, stop, score, strand, phase, attributes, tags = read_gtf_line(line)
            gene_id = attributes['gene_id']
            gene_name = attributes['gene_name']
            
            if feature == 'gene':
                if chr not in chromosomes:
                    chromosome = Chromosome.from_name(chr)
                    if chromosome is None:
                        chromosome = Chromosome()
                        chromosome.name = chr
                    chromosomes[chr] = chromosome
                    db_session.add(chromosome)
                else:
                    chromosome = chromosomes[chr]

                if Gene.from_accession(gene_id) is None:
                    gene = Gene()
                    gene.accession = gene_id
                    gene.name = gene_name
                    # gene.strand = strand
                    genes[gene_id] = gene
                    db_session.add(gene)
                    chromosome.genes.append(gene)
                    
            elif feature == 'transcript':
                transcript_id = attributes['transcript_id']
                transcript_name = attributes['transcript_name']
                if Transcript.from_accession(transcript_id) is None:
                    appris = APPRIS.NONE
                    start_nf, end_nf = False, False
                    for tag in tags:
                        if 'appris_principal' in tag:
                            appris = APPRIS.PRINCIPAL
                        if 'appris_alternative' in tag:
                            appris = APPRIS.ALTERNATIVE
                        start_nf = start_nf or 'start_NF' in tag
                        end_nf = end_nf or 'end_NF' in tag
                    transcript = GencodeTranscript(
                        accession = transcript_id,
                        name = transcript_name,
                        strand = strand,
                        appris = appris,
                        start_nf = start_nf,
                        end_nf = end_nf
                    )
                    genes[gene_id].transcripts.append(transcript)
                    transcripts[transcript_id] = transcript
                    db_session.add(transcript)
                
            elif feature == 'exon':
                exon_id = attributes['exon_id']
                if Exon.from_accession(exon_id) is None:
                    exon = GencodeExon(
                        accession = exon_id,
                        start = start,
                        stop = stop,
                        transcript = transcripts[attributes['transcript_id']]
                    )
                    db_session.add(exon)
    
    # calculate the coordinates of each exon relative to the sequence of its parent transcript
    for transcript in transcripts.values():
        exon_to_genomic_loc = [SingleInterval(exon.start - 1, exon.stop, transcript.strand) for exon in transcript.exons]
        transcript_genomic_loc = exon_to_genomic_loc[0]
        for exon_genomic_loc in exon_to_genomic_loc:
            transcript_genomic_loc = transcript_genomic_loc.union(exon_genomic_loc)
        for i, exon in enumerate(transcript.exons):
            # TODO: is it faster to just loop through transcript's exons and count off lengths manually?
            exon_transcript_loc = exon_to_genomic_loc[i].location_relative_to(transcript_genomic_loc)
            exon.transcript_start = exon_transcript_loc.start + 1
            exon.transcript_stop = exon_transcript_loc.end
        transcript.exons.sort(key=attrgetter('transcript_start'))

    db_session.commit() #Attempt to commit all the records
    
def load_transcript_fasta(transcript_fasta):
    for i, record in enumerate(SeqIO.parse(transcript_fasta, 'fasta')):
        if i % 5000 == 0:
            print(f'read {i//1000}k entries...')
        fields = record.id.split('|')
        transcript_name = fields[0]
        orf_coords = [field[4:] for field in fields if field.startswith('CDS:')][0]
        orf_start = int(orf_coords.split('-')[0])
        orf_end = int(orf_coords.split('-')[1])
        sequence = str(record.seq)

        try:
            statement = select(Transcript).filter(Transcript.accession == transcript_name)
            result = db_session.execute(statement).one()
            transcript = result[0]
        except NoResultFound:
            pass
            # logging.warning(f'could not get transcript {transcript_name} from database')
        else:
            transcript.sequence = sequence
            # check if ORF already exists
            def has_same_range(orf):
                return (orf.transcript_start, orf.transcript_stop) == (orf_start, orf_end)
            if not any(has_same_range(orf) for orf in transcript.orfs):
                orf = ORF()
                orf.transcript = transcript
                orf.transcript_start, orf.transcript_stop = orf_start, orf_end
                orf.start = orf.nucleotides[0].coordinate
                orf.stop = orf.nucleotides[-1].coordinate
                if orf.start > orf.stop:
                    orf.start, orf.stop = orf.stop, orf.start
                db_session.add(orf)
            db_session.commit() #Attempt to commit all the records
    # db_session.commit() #Attempt to commit all the records

def load_translation_fasta(translation_fasta):
    for i, record in enumerate(SeqIO.parse(translation_fasta, 'fasta')):
        if i % 5000 == 0:
            print(f'read {i//1000}k entries...')
        fields = record.id.split('|')
        protein_id, transcript_id = fields[:2]
        sequence = str(record.seq)
        
        try:
            statement = select(Transcript).filter(Transcript.accession == transcript_id)
            result = db_session.execute(statement).one()
            transcript = result[0]
        except NoResultFound:
            pass
            # logging.warning(f'could not get transcript {transcript_id} from database')
        else:
            def has_same_sequence(protein):
                return protein.sequence == sequence
            if not any(has_same_sequence(orf.protein) for orf in transcript.orfs if orf.protein):
                protein = Protein()
                protein.sequence = sequence
                protein.orf = transcript.orfs[0]
                db_session.add(protein)
    db_session.commit()

path = '/home/redox/sheynkman-lab/biosurfer/data/biosurfer_demo_data/'
gtf_file = 'chr19.gtf'
tx_file = 'gencode.v38.pc_transcripts.fa'
tl_file = 'gencode.v38.pc_translations.fa'
# gtf_file = 'gencode.v38.annotation.gtf.toy'
# tx_file = 'gencode.v38.pc_transcripts.fa.toy'
# tl_file = 'gencode.v38.pc_translations.fa.toy'

#%%
start = time.time()
load_data_from_gtf(path + 'chr19.gtf')
end = time.time()
print(f"Time to load gtf file\t{end - start:0.3g}")

#%%
start = time.time()
load_data_from_gtf(path + 'chr22.gtf')
end = time.time()
print(f"Time to load gtf file\t{end - start:0.3g}")

#%%
start = time.time()
load_transcript_fasta(path + tx_file)
end = time.time()
print(f"time to load transcript fasta\t{end - start:0.3g} seconds")
print(f"time to load transcript fasta\t{(end - start)/60:0.3g} minutes")

#%%
start = time.time()
load_translation_fasta(path + tl_file)
end = time.time()
print(f"time to load translation fasta\t{end - start:0.3g} seconds")
# %%
