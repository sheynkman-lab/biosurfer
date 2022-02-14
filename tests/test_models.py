from itertools import chain, product

import hypothesis.strategies as st
import pytest
from Bio.Seq import Seq
from biosurfer.core.constants import (START_CODON, STOP_CODONS, AminoAcid,
                                      Strand)
from biosurfer.core.models.biomolecules import (ORF, Exon, Gene, Protein,
                                                Transcript)
from biosurfer.core.models.features import ProteinFeature
from hypothesis import given, note
from hypothesis.control import assume
from sqlalchemy.sql.expression import select

NT_ALPHABET = 'ACGT'

codons = (''.join(triplet) for triplet in product(NT_ALPHABET, repeat=3))
STOP_CODONS = tuple(STOP_CODONS)
non_stop_codons = tuple(codon for codon in codons if codon not in STOP_CODONS)

def nt_seqs(min_size=0, max_size=None):
    return st.text(alphabet=NT_ALPHABET, min_size=min_size, max_size=max_size)

@st.composite
def orf_seqs_missing_stop(draw, min_aa_length=1):
    codon_iter = draw(
        st.iterables(
            st.sampled_from(non_stop_codons),
            min_size = min_aa_length-1,
        )
    )
    seq = START_CODON + ''.join(codon_iter)
    return seq

@st.composite
def orf_seqs(draw, min_aa_length=1):
    return draw(orf_seqs_missing_stop(min_aa_length)) + draw(st.sampled_from(STOP_CODONS))

@st.composite
def transcript(draw, coding=False, min_size=1):
    if coding:
        utr5_seq = draw(nt_seqs())
        orf_seq = draw(orf_seqs())
        utr3_seq = draw(nt_seqs())
        tx_seq = utr5_seq + orf_seq + utr3_seq

        protein_seq = str(Seq(orf_seq).translate())
        protein = Protein(
            sequence = protein_seq
        )
        orf_tx_start = len(utr5_seq) + 1
        orf_tx_stop = len(utr5_seq) + len(orf_seq)
        orf = ORF(
            transcript_start = orf_tx_start,
            transcript_stop = orf_tx_stop,
            protein = protein
        )
    else:
        tx_seq = draw(nt_seqs(min_size=min_size))
        protein = None
        orf = None
    
    length = len(tx_seq)

    strand = draw(st.sampled_from(Strand))
    exons = []
    if length < 6:
        exon_tx_stops = []
    else:
        exon_tx_stops = draw(
            st.lists(
                st.integers(
                    min_value = 4,
                    max_value = length - 2
                ),
                unique = True,
                max_size = len(tx_seq)//2
            )
        )
    exon_tx_stops.append(length)
    exon_tx_stops.sort()

    intron_lengths = draw(
        st.lists(
            st.integers(min_value=31),
            min_size = len(exon_tx_stops)
        )
    )

    exon_tx_start = 1
    pos = draw(st.integers(min_value=length + 2*sum(intron_lengths)))
    for i, exon_tx_stop in enumerate(exon_tx_stops):
        exon_length = exon_tx_stop - exon_tx_start
        exon = Exon(
            position = i + 1,
            transcript_start = exon_tx_start,
            transcript_stop = exon_tx_stop
        )
        if strand is Strand.MINUS:
            exon.stop = pos
            exon.start = pos - exon_length
            pos = exon.start - intron_lengths[i]
        else:
            exon.start = pos
            exon.stop = pos + exon_length
            pos = exon.stop + intron_lengths[i]

        exon_tx_start = exon_tx_stop + 1
        exons.append(exon)
    
    tx = Transcript(
        name = 'TEST-1',
        gene = Gene(name='TEST', chromosome_id='chrZ'),
        strand = strand,
        sequence = tx_seq,
        exons = exons,
        orfs = [orf] if orf else []
    )
    return tx

@pytest.fixture(scope='module')
def transcript_getter(session):
    db_transcripts = session.execute(select(Transcript)).scalars().all()
    def _get_transcript(data):
        return data.draw(
            st.one_of(
                transcript(),
                st.sampled_from(db_transcripts)
            )
        )
    return _get_transcript

@pytest.fixture(scope='module')
def coding_transcript_getter(session):
    db_coding_transcripts = session.execute(
        select(Transcript).
        select_from(ORF).
        join(ORF.transcript)
    ).scalars().all()
    def _get_coding_transcript(data):
        return data.draw(
            st.one_of(
                transcript(coding=True),
                st.sampled_from(db_coding_transcripts)
            )
        )
    return _get_coding_transcript

@pytest.fixture(scope='module')
def feature_getter(session):
    db_features = session.execute(select(ProteinFeature)).scalars().all()
    def _get_feature(data):
        return data.draw(
            st.sampled_from(db_features)
        )
    return _get_feature

########################
### TESTS BEGIN HERE ###
########################

@given(data=st.data())
def test_transcript_nucleotides_match_sequence(data, transcript_getter):
    transcript = transcript_getter(data)
    assert ''.join(str(nt.base) for nt in transcript.nucleotides) == transcript.sequence

@given(data=st.data())
def test_exon_lengths_correct(data, transcript_getter):
    transcript = transcript_getter(data)
    assert sum(exon.length for exon in transcript.exons) == transcript.length

@given(data=st.data())
def test_exon_sequences_correct(data, transcript_getter):
    transcript = transcript_getter(data)
    assert ''.join(exon.sequence for exon in transcript.exons) == transcript.sequence

@given(data=st.data())
def test_exon_coords_correct(data, transcript_getter):
    transcript = transcript_getter(data)
    for exon in transcript.exons:
        assert exon.stop - exon.start == exon.transcript_stop - exon.transcript_start

@given(data=st.data())
def test_junction_coords_correct(data, transcript_getter):
    transcript = transcript_getter(data)
    for junction in transcript.junctions:
        upstream_exon, downstream_exon = transcript.get_exons_from_junction(junction)
        donor_nt = transcript.get_nucleotide_from_coordinate((junction.donor.coordinate))
        acceptor_nt = transcript.get_nucleotide_from_coordinate((junction.acceptor.coordinate))
        note(f'junction: {junction}')
        note(f'upstream exon: {upstream_exon}, {upstream_exon.start}-{upstream_exon.stop}')
        note(f'downstream exon: {downstream_exon}, {downstream_exon.start}-{downstream_exon.stop}')
        assert bool(junction) and donor_nt.exon is upstream_exon and acceptor_nt.exon is downstream_exon

@given(data=st.data())
def test_orf_utr_lengths_correct(data, coding_transcript_getter):
    transcript = coding_transcript_getter(data)
    for orf in transcript.orfs:
        utr5_length = orf.utr5.length if orf.utr5 else 0
        utr3_length = orf.utr3.length if orf.utr3 else 0
        assert utr5_length + orf.length + utr3_length == transcript.length

@given(data=st.data())
def test_orf_utr_sequences_correct(data, coding_transcript_getter):
    transcript = coding_transcript_getter(data)
    for orf in transcript.orfs:
        utr5_sequence = orf.utr5.sequence if orf.utr5 else ''
        utr3_sequence = orf.utr3.sequence if orf.utr3 else ''
        assert utr5_sequence + orf.sequence + utr3_sequence == transcript.sequence

@given(data=st.data())
def test_nucleotide_has_exon(data, transcript_getter):
    transcript = transcript_getter(data)
    for exon in transcript.exons:
        for nt in exon.nucleotides:
            assert nt.exon is exon

@given(data=st.data())
def test_protein_residues_match_sequence(data, coding_transcript_getter):
    transcript = coding_transcript_getter(data)
    protein = transcript.protein
    assert ''.join(str(res.amino_acid) for res in protein.residues) == protein.sequence

@given(data=st.data())
def test_nucleotide_has_residue(data, coding_transcript_getter):
    transcript = coding_transcript_getter(data)
    orf = transcript.primary_orf
    note(f'tx seq: {transcript.sequence}')
    note(f'orf seq: {orf.sequence}')
    note(f'protein seq: {orf.protein.sequence}')
    assert list(chain.from_iterable((res, res, res) for res in orf.protein.residues)) == [nt.residue for nt in orf.nucleotides]

@given(data=st.data())
def test_noncoding_nucleotide_has_no_residue(data, transcript_getter):
    transcript = transcript_getter(data)
    assume(len(transcript.orfs) == 0)
    note(f'tx seq: {transcript.sequence}')
    for nt in transcript.nucleotides:
        assert nt.residue is None

@given(data=st.data())
def test_residue_has_nucleotide(data, coding_transcript_getter):
    transcript = coding_transcript_getter(data)
    orf = transcript.primary_orf
    note(f'tx seq: {transcript.sequence}')
    note(f'orf seq: {orf.sequence}')
    note(f'protein seq: {orf.protein.sequence}')
    assert list(chain.from_iterable(res.codon for res in orf.protein.residues)) == orf.nucleotides

@given(data=st.data())
def test_feature_residues_correct(data, feature_getter):
    feature = feature_getter(data)
    protein = feature.protein
    note(f'protein seq: {protein.sequence}')
    note(f'feature seq: {feature.sequence}')
    assert [(res.position, res.amino_acid) for res in feature.residues] == [(feature.protein_start + i, AminoAcid(aa)) for i, aa in enumerate(feature.sequence)]

# TODO: test junctions
