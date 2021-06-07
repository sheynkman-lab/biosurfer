gene = gd['PAX5']

import random

for i in range(100):

    grp = isogroup.GencodeGroup(random.sample(gene.orfs, 2))
    # grp = isogroup.GencodeGroup([gene['PAX5-218'], gene['PAX5-214']])

    grp.enter()

    def get_overlapping_pos_in_second_orf_if_exists(pos, orf):
        abs_coord = pos.coord
        for pos in orf.chain:
            if pos.coord == abs_coord:
                return pos
        return None

    def determine_frame_of_trans_ordinal(frm):
        """Based on shift in positions (0, 1, 2), determine 1-based frame shift."""
        if frm == 0:
            return 1
        elif frm == 2:
            return 2
        elif frm == 1:
            return 3

    for res1 in grp.repr_orf.res_chain:
        frmr1 = isofeature.FrameResidue(1)
        map = isomap.Res_Featr_Mapping(res1, frmr1, grp=grp)
        # find the overlapping res in second orf, if exists, then create frm
        pos1 = res1.p1
        pos2 = get_overlapping_pos_in_second_orf_if_exists(pos1, grp.other_orfs[0])
        if pos2 and pos2.res:
            res2 = pos2.res
            frm_ord = determine_frame_of_trans_ordinal(pos2.frm)
            frmr2 = isofeature.FrameResidue(frm_ord)
            map = isomap.Res_Featr_Mapping(res2, frmr2, grp=grp)

    print(grp.repr_orf.name)
    print(grp.repr_orf.current_grp)
    print(''.join([res.aa for res in grp.repr_orf.res_chain]))
    print(''.join([str(res.map.featr.ord) if res.map else '-' for res in grp.repr_orf.res_chain]))
    print(grp.other_orf.name)
    print(grp.other_orf.current_grp)
    print(''.join([res.aa for res in grp.other_orf.res_chain]))
    print(''.join([str(res.map.featr.ord) if res.map else '-' for res in grp.other_orf.res_chain]))
    print('\n\n')
    grp.exit()


