import molsysmt as msm

def is_head_up(lipid):

    position_head = msm.get(lipid, element='atom', selection='atom_name=="P"', coordinates=True)
    position_com = msm.structure.get_center(lipid)

    if position_head[0,0,2]>=position_com[0,0,2]:
        return True
    else:
        return False

