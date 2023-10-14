import molsysmt as msm
import omembrane as omem
from parmed.charmm import CharmmParameterSet
import requests
import tarfile
from tqdm import tqdm
import os
import shutil

toppar_dir = '../omembrane/data/charmm_topology/toppar/'
toppar_files = [toppar_dir+'top_all36_lipid.rtf',
                toppar_dir+'par_all36_lipid.prm']

params = CharmmParameterSet(*toppar_files)

def work_lipid(lipid):
    lipid = msm.structure.align_principal_axes(lipid, axes=[[0,0,1],[1,0,0],[0,1,0]])
    lipid = msm.structure.center(lipid, selection='all', center_of_selection='atom_name=="P"')
    if not omem.lipid.is_head_up(lipid):
        lipid = msm.structure.flip(lipid)
    return lipid

for lipid_name in params.residues:

    if lipid_name in omem.lipid.name:
        lipid_code = lipid_name.lower()
    else:
        print(lipid_name, 'not in name lists')
        continue

    print(lipid_name, lipid_code)

    parmed_structure = params.residues[lipid_name].to_structure()
    molsys = msm.convert(parmed_structure)
    molsys = msm.remove(molsys, structure_indices='all')

    if os.path.exists(lipid_code):
        continue

    os.mkdir(lipid_code)

    msm.convert(parmed_structure, lipid_code+'/'+lipid_code+'.psf')

    src_web = 'https://www.charmm-gui.org/archive/csml/'+lipid_code+'.pdb'
    response = requests.get(src_web, stream=True, auth=('user', 'pass'))


    if response.status_code == 200:
        with open(lipid_code+'/'+lipid_code+'.pdb', 'w') as fff:
            fff.write(response.text)
    else:
        print(src_web)
        print('     ...without pdb')

    src_web = 'https://www.charmm-gui.org/archive/lipid/'+lipid_code+'.tar.gz'
    response = requests.get(src_web, stream=True, auth=('user', 'pass'))

    if response.status_code == 200:
        with open('foo.tar.gz', 'wb') as fff:
            fff.write(response.raw.read())
        file = tarfile.open('foo.tar.gz')
        file.extractall(path="./foo")

        for ii in tqdm(range(1,1001)):
            filename = './foo/'+lipid_code+'/conf1/'+lipid_code+'_'+str(ii)+'.crd'
            aux_lipid = msm.convert(filename)
            aux_lipid = work_lipid(aux_lipid)
            msm.append_structures(molsys, aux_lipid)

        for ii in tqdm(range(1, 1001)):
            filename = './foo/'+lipid_code+'/conf2/'+lipid_code+'_'+str(ii)+'.crd'
            aux_lipid = msm.convert(filename)
            aux_lipid = work_lipid(aux_lipid)
            msm.append_structures(molsys, aux_lipid)

        os.remove('foo.tar.gz')
        shutil.rmtree('./foo')

        msm.convert(molsys, lipid_code+'/'+lipid_code+'.msmpk')
        coordinates = msm.get(molsys, coordinates=True)
        msm.convert(coordinates, lipid_code+'/'+lipid_code+'_conformations.xyznpy')

        del(coordinates, aux_lipid)

    else:

        print(src_web)
        print('     ...without coordinates')

    del(parmed_structure, molsys)

    print(' ')

