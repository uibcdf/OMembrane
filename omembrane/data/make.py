import molsysmt as msm
import omembrane as omem
from parmed.charmm import CharmmParameterSet
import requests
import tarfile
from tqdm import tqdm
import os
import shutil

# Auxiliary functions

def work_lipid(lipid):
    lipid = msm.structure.align_principal_axes(lipid, axes=[[0,0,1],[1,0,0],[0,1,0]])
    lipid = msm.structure.center(lipid, selection='all', center_of_selection='atom_name=="P"')
    if not omem.lipid.is_head_up(lipid):
        lipid = msm.structure.flip(lipid)
    return lipid

# Charmm forcefield

toppar_dir = 'charmm_topology/toppar/'
toppar_files = [toppar_dir+'top_all36_lipid.rtf',
                toppar_dir+'par_all36_lipid.prm']

params = CharmmParameterSet(*toppar_files)

# POPC

lipid_name = 'popc'
print(f'...{lipid_name}')

url_coors = f'https://www.charmm-gui.org/archive/lipid/{lipid_name}.tar.gz'

response = requests.get(url_coors)

if response.status_code == 200:
    with open(f'{lipid_name}.tar.gz', 'wb') as file:
        file.write(response.content)
    with tarfile.open(f'{lipid_name}.tar.gz', 'r:gz') as tar:
        nombres_archivos = tar.getnames()
        tar.extractall()
else:
    print("Error al descargar el archivo tar.gz. Código de estado HTTP:", response.status_code)

parmed_structure = params.residues[lipid_name.upper()].to_structure()
molsys = msm.convert(parmed_structure)
molsys = msm.remove(molsys, structure_indices='all')

for ii in tqdm(range(1,1001)):
    filename = f'{lipid_name}/conf1/{lipid_name}_{ii}.crd'
    aux_lipid = msm.convert(filename)
    aux_lipid = work_lipid(aux_lipid)
    msm.append_structures(molsys, aux_lipid)

for ii in tqdm(range(1,1001)):
    filename = f'{lipid_name}/conf2/{lipid_name}_{ii}.crd'
    aux_lipid = msm.convert(filename)
    aux_lipid = work_lipid(aux_lipid)
    msm.append_structures(molsys, aux_lipid)

msm.convert(molsys, f'{lipid_name}.h5msm')

shutil.rmtree(f'{lipid_name}')
os.remove(f'{lipid_name}.tar.gz')

