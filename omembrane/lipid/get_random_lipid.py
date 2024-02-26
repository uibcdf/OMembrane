import molsysmt as msm
import numpy as np
import sys
import gc

if sys.version_info[1]==10:
    from importlib.resources import files
    def path(package, file):
        return files(package).joinpath(file)
elif sys.version_info[1] in (8,9):
    from pathlib import PurePath
    parent = PurePath(__file__).parent
    def path(package, file):
        data_dir = package.split('.')[-1]
        return parent.joinpath('../data/'+data_dir+'/'+file).__str__()

rng = np.random.default_rng()

def get_random_lipid(name, n_lipids=1):

    output=[]

    if name=='POPC':

        lipid = msm.convert(path('omembrane.data','popc.h5msm'))

        n_structures = msm.get(lipid, n_structures=True)

        structure_indices = rng.integers(0, n_structures, n_lipids)

        for ii in structure_indices:

            output.append(msm.extract(lipid, structure_indices=ii))

        del(lipid)
        gc.collect()

    if len(output)==1:
        output = output[0]

    return output

