from itertools import chain, product
from typing import TYPE_CHECKING, Callable, Dict

from Bio.Seq import Seq
from biosurfer.core.constants import Strand
from biosurfer.core.models import (ORF, Exon, Gene, Protein, ProteinFeature,
                                   Transcript)

import hypothesis.strategies as st
from hypothesis import given, note

NT_ALPHABET = 'ACGT'

codons = (''.join(triplet) for triplet in product(NT_ALPHABET, repeat=3))
start_codon = 'ATG'
stop_codons = ('TGA', 'TAA', 'TAG')
non_stop_codons = tuple(codon for codon in codons if codon not in stop_codons)

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
    seq = start_codon + ''.join(codon_iter)
    return seq

@st.composite
def orf_seqs(draw, min_aa_length=1):
    return draw(orf_seqs_missing_stop(min_aa_length)) + draw(st.sampled_from(stop_codons))

@st.composite
def coding_transcript(draw, add_exons=True, add_protein=True):
    utr5_seq = draw(nt_seqs())
    orf_seq = draw(orf_seqs())
    utr3_seq = draw(nt_seqs())
    tx_seq = utr5_seq + orf_seq + utr3_seq
    length = len(tx_seq)

    if add_protein:
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
        protein = None
        orf = None

    strand = draw(st.sampled_from(Strand))
    exons = []
    if add_exons:
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
            st.iterables(
                st.integers(min_value=31),
                min_size = len(exon_tx_stops)
            )
        )

        exon_tx_start = 1
        exon_genome_start = draw(st.integers(min_value=0))
        for i, exon_tx_stop in enumerate(exon_tx_stops, start=1):
            exon = Exon(
                position = i,
                transcript_start = exon_tx_start,
                transcript_stop = exon_tx_stop,
                start = exon_genome_start,
                stop = exon_genome_start + exon_tx_stop - exon_tx_start
            )
            exons.append(exon)
            exon_tx_start = exon_tx_stop + 1
            exon_genome_start = exon.stop + next(intron_lengths)

    
    tx = Transcript(
        name = 'TEST-1',
        strand = strand,
        sequence = tx_seq,
        exons = exons,
        orfs = [orf] if orf else []
    )
    return tx

@given(coding_transcript(add_protein=False))
def test_exon_lengths_correct(transcript):
    assert sum(exon.length for exon in transcript.exons) == transcript.length

@given(coding_transcript(add_protein=False))
def test_exon_sequences_correct(transcript):
    assert ''.join(exon.sequence for exon in transcript.exons) == transcript.sequence

@given(coding_transcript(add_protein=False))
def test_nucleotide_has_exon(transcript):
    assert all(nt.exon is exon for exon in transcript.exons for nt in exon.nucleotides)

@given(coding_transcript(add_exons=False))
def test_nucleotide_has_residue(transcript):
    orf = transcript.primary_orf
    note(f'tx seq: {transcript.sequence}')
    note(f'orf seq: {orf.sequence}')
    note(f'protein seq: {orf.protein.sequence}')
    assert list(chain.from_iterable((res, res, res) for res in orf.protein.residues)) == [nt.residue for nt in orf.nucleotides]

@given(coding_transcript(add_exons=False))
def test_residue_has_nucleotide(transcript):
    orf = transcript.primary_orf
    note(f'tx seq: {transcript.sequence}')
    note(f'orf seq: {orf.sequence}')
    note(f'protein seq: {orf.protein.sequence}')
    assert list(chain.from_iterable(res.codon for res in orf.protein.residues)) == orf.nucleotides
