import pytest
from biosurfer.core.alignments import Alignment
from biosurfer.core.alignments_new import ProteinAlignment, CodonAlignCat
from biosurfer.core.models.biomolecules import Transcript
from biosurfer.core.models.nonpersistent import Position

def test_codon_alignment_blocks(session, alignment_case):
    print(session.get_bind())
    anchor: 'Transcript' = Transcript.from_name(session, alignment_case['anchor'])
    other: 'Transcript' = Transcript.from_name(session, alignment_case['other'])
    chr = anchor.gene.chromosome_id
    strand = anchor.strand
    aln = ProteinAlignment.from_proteins(anchor.protein, other.protein)
    for block in aln.blocks:
        assert block.anchor_range or block.other_range
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
        elif block.category in {CodonAlignCat.EDGE_MATCH, CodonAlignCat.EDGE_MISMATCH}:
            assert len(block.anchor_range) == len(block.other_range) == 1
            anchor_edge_res = anchor.protein.residues[block.anchor_range[0]]
            other_edge_res = other.protein.residues[block.other_range[0]]
            assert anchor_edge_res.codon[2:] == other_edge_res.codon[2:] or anchor_edge_res.codon[:2] == other_edge_res.codon[:2]
            assert (anchor_edge_res.amino_acid == other_edge_res.amino_acid) ^ (block.category is CodonAlignCat.EDGE_MISMATCH)
        elif block.category is CodonAlignCat.COMPLEX:
            assert len(block.anchor_range) == len(block.other_range) == 1
            anchor_res = anchor.protein.residues[block.anchor_range[0]]
            other_res = other.protein.residues[block.other_range[0]]
            overlap = set(anchor_res.codon) & set(other_res.codon)
            assert overlap == {anchor_res.codon[1]} or overlap == {other_res.codon[1]}
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

# @pytest.mark.skip(reason='old Alignment is deprecated')
# def test_alignment_full(session, alignment_case):
#     print(session.get_bind())
#     anchor = alignment_case['anchor']
#     other = alignment_case['other']
#     expected = alignment_case['full']
#     txs = Transcript.from_names(session, (anchor, other))
#     aln = Alignment(txs[anchor].protein, txs[other].protein)
#     assert aln.full == expected

# @pytest.mark.skip(reason='feature cases may be outdated')
# def test_feature_alignment(session, feature_case):
#     print(session.get_bind())
#     anchor = feature_case['anchor']
#     other = feature_case['other']
#     txs = Transcript.from_names(session, (anchor, other))
#     aln = Alignment(txs[anchor].protein, txs[other].protein)
#     for feat_aln, expected in zip((aln.project_feature(feat)[1] for feat in aln.anchor.features), feature_case['expected']):
#         assert feat_aln.full == expected
