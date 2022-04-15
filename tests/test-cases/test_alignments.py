import pytest
from biosurfer.core.alignments import CodonAlignment, CodonAlignCat, ProteinAlignment, SeqAlignCat
from biosurfer.core.models.biomolecules import Transcript
from biosurfer.core.models.nonpersistent import Position
from hypothesis import given, note, strategies as st
from more_itertools import one

def test_codon_alignment_blocks(session, alignment_case):
    anchor: 'Transcript' = Transcript.from_name(session, alignment_case['anchor'])
    other: 'Transcript' = Transcript.from_name(session, alignment_case['other'])
    chr = anchor.gene.chromosome_id
    strand = anchor.strand
    aln = CodonAlignment.from_proteins(anchor.protein, other.protein)
    for block in aln.blocks:
        print(block)
    for block in aln.blocks:
        if block.category is CodonAlignCat.MATCH:
            assert [anchor.protein.residues[i].codon for i in block.anchor_range] == [other.protein.residues[i].codon for i in block.other_range]
        elif block.category in {CodonAlignCat.FRAME_AHEAD, CodonAlignCat.FRAME_BEHIND}:
            assert len(block.anchor_range) == len(block.other_range)
            pr_coords = zip(block.anchor_range, block.other_range)
            if block.category is CodonAlignCat.FRAME_AHEAD:
                for a, o in pr_coords:
                    assert anchor.protein.residues[a].codon[:2] == other.protein.residues[o].codon[-2:]
            else:
                for a, o in pr_coords:
                    assert anchor.protein.residues[a].codon[-2:] == other.protein.residues[o].codon[:2]
        elif block.category in {CodonAlignCat.EDGE, CodonAlignCat.COMPLEX}:
            assert len(block.anchor_range) == len(block.other_range) == 1
            anchor_res = anchor.protein.residues[block.anchor_range[0]]
            other_res = other.protein.residues[block.other_range[0]]
            if block.category is CodonAlignCat.EDGE:
                difference = set(nt.coordinate for nt in anchor_res.codon) ^ set(nt.coordinate for nt in other_res.codon)
                assert difference == {anchor_res.codon[0].coordinate, other_res.codon[0].coordinate} or difference == {anchor_res.codon[-1].coordinate, other_res.codon[-1].coordinate}
        elif block.category is CodonAlignCat.DELETION:
            assert not block.other_range
            for i in block.anchor_range:
                with pytest.raises(KeyError):
                    for nt in anchor.protein.residues[i].codon:
                        other.get_transcript_coord_from_genome_coord(Position(chr, strand, nt.coordinate))
        elif block.category is CodonAlignCat.INSERTION:
            assert not block.anchor_range
            for i in block.other_range:
                with pytest.raises(KeyError):
                    for nt in other.protein.residues[i].codon:
                        anchor.get_transcript_coord_from_genome_coord(Position(chr, strand, nt.coordinate))
        elif block.category is CodonAlignCat.UNTRANSLATED:
            assert not block.other_range and (block.other_range.start == 0 or block.other_range.start == other.protein.length)
        elif block.category is CodonAlignCat.TRANSLATED:
            assert not block.anchor_range and (block.anchor_range.start == 0 or block.anchor_range.start == anchor.protein.length)

def test_protein_alignment_blocks(session, alignment_case):
    anchor: 'Transcript' = Transcript.from_name(session, alignment_case['anchor'])
    other: 'Transcript' = Transcript.from_name(session, alignment_case['other'])
    aln = ProteinAlignment.from_proteins(anchor.protein, other.protein)
    for block in aln.blocks:
        print(block)
    for block in aln.blocks:
        if block.category is SeqAlignCat.MATCH:
            assert anchor.protein.sequence[block.anchor_range.start:block.anchor_range.stop] == other.protein.sequence[block.other_range.start:block.other_range.stop]
        elif block.category is SeqAlignCat.SUBSTITUTION:
            assert anchor.protein.sequence[block.anchor_range.start:block.anchor_range.stop] != other.protein.sequence[block.other_range.start:block.other_range.stop]
        elif block.category is SeqAlignCat.DELETION:
            assert not block.other_range
        elif block.category is SeqAlignCat.INSERTION:
            assert not block.anchor_range
        else:
            raise ValueError(block.category)

def test_alignment_coordinate_projection(session, alignment_case):
    anchor: 'Transcript' = Transcript.from_name(session, alignment_case['anchor'])
    other: 'Transcript' = Transcript.from_name(session, alignment_case['other'])
    aln = CodonAlignment.from_proteins(anchor.protein, other.protein)
    print(aln.blocks)
    for anchor_coord in range(anchor.protein.length):
        other_coord = aln.project_coordinate(anchor_coord)
        if other_coord is None:
            assert one(aln.anchor_blocks.at(anchor_coord)).data.category in {CodonAlignCat.DELETION, CodonAlignCat.UNTRANSLATED}
        else:
            assert anchor_coord == aln.project_coordinate(other_coord, from_anchor=False)

@given(data=st.data())
def test_alignment_range_projection(data, session, alignment_case):
    anchor: 'Transcript' = Transcript.from_name(session, alignment_case['anchor'])
    other: 'Transcript' = Transcript.from_name(session, alignment_case['other'])
    aln = CodonAlignment.from_proteins(anchor.protein, other.protein)
    note(aln.blocks)
    anchor_start = data.draw(st.integers(min_value=0, max_value=anchor.protein.length - 1))
    anchor_stop = data.draw(st.integers(min_value=anchor_start + 1, max_value=anchor.protein.length))
    other_range = aln.project_range(anchor_start, anchor_stop)
    if other_range:
        note(f'{other_range = }')
        anchor_range = aln.project_range(other_range.start, other_range.stop, from_anchor=False)
        assert anchor_start <= anchor_range.start and anchor_range.stop <= anchor_stop
        assert other_range.start == 0 or aln.project_coordinate(other_range.start - 1, from_anchor=False) not in range(anchor_start, anchor_stop)
        assert other_range.stop == other.protein.length or aln.project_coordinate(other_range.stop, from_anchor=False) not in range(anchor_start, anchor_stop)

def test_alignment_feature_projection(session, alignment_case):
    anchor: 'Transcript' = Transcript.from_name(session, alignment_case['anchor'])
    other: 'Transcript' = Transcript.from_name(session, alignment_case['other'])
    aln = CodonAlignment.from_proteins(anchor.protein, other.protein)
    print(aln.blocks)
    for feature in anchor.protein.features:
        proj_feat = aln.project_feature(feature)
        if proj_feat:
            print(f'{feature}\t{feature.sequence}')
            print(f'{proj_feat}\t{proj_feat.sequence}')
            for residue in (residue for residue in proj_feat.residues if residue not in proj_feat.altered_residues):
                anchor_residue_pos = aln.project_coordinate(residue.position - 1, from_anchor=False)
                assert residue.amino_acid is anchor.protein.residues[anchor_residue_pos].amino_acid
        else:
            assert not any(interval.data.other_range for interval in aln.anchor_blocks.overlap(feature.protein_start - 1, feature.protein_stop))
