import molsysmt as msm
import omembrane as omem
from tqdm import tqdm

def work_lipid(lipid):
    lipid = msm.structure.align_principal_axes(lipid, axes=[[0,0,1],[1,0,0],[0,1,0]])
    lipid = msm.structure.center(lipid, selection='all', center_of_selection='atom_name=="P"')
    if not omem.lipid.is_head_up(lipid):
        lipid = msm.structure.flip(lipid)
    return lipid

filename = './charmmgui/conf1/popc_1.crd'
lipid = msm.convert(filename)

for ii in tqdm(range(2, 1000)):
    filename = './charmmgui/conf1/popc_'+str(ii)+'.crd'
    aux_lipid = msm.convert(filename)
    aux_lipid = work_lipid(aux_lipid)
    msm.append_structures(lipid, aux_lipid)

for ii in tqdm(range(1, 1000)):
    filename = './charmmgui/conf2/popc_'+str(ii)+'.crd'
    aux_lipid = msm.convert(filename)
    aux_lipid = work_lipid(aux_lipid)
    msm.append_structures(lipid, aux_lipid)


msm.convert(lipid, 'popc.msmpk')
coordinates = msm.get(lipid, coordinates=True)
msm.convert(coordinates, 'popc_conformations.xyznpy')

