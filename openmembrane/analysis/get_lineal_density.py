import molsysmt as msm
from molsysmt import pyunitwizard as puw
import numpy as np

def get_lineal_density(molecular_system, selection, axis=[0,0,1], bins=100, range=None, normalized=False, syntax='MolSysMT'):


    atom_indices = msm.select(molecular_system, selection=selection, syntax=syntax)
    
    n_structures = msm.get(molecular_system, n_structures=True)

    if range is None:
        
        iterator = msm.Iterator(molecular_system, selection = atom_indices, coordinates = True)

        min_total = np.inf
        max_total = -np.inf
        
        for aux_coordinates in iterator:
            proj = np.dot(aux_coordinates, axis)
            min_coordinate = puw.get_value(proj.min())
            max_coordinate = puw.get_value(proj.max())

            if min_coordinate < min_total:
                min_total = min_coordinate

            if max_coordinate > max_total:
                max_total = max_coordinate

        range = [min_total, max_total]


    bin_length = (range[1]-range[0])/bins

    frequency = np.zeros(bins)

    iterator = msm.Iterator(molecular_system, selection = atom_indices, coordinates = True)

    for aux_coordinates in iterator:
        proj = np.dot(aux_coordinates, axis)
        proj = proj.flatten()
        proj = puw.get_value(proj)

        for ii in proj:
            bin_index = int(np.floor((ii - range[0])/bin_length))
            if bin_index==bins:
                bin_index=bins-1
            frequency[bin_index]+=1

    if not normalized:
        density = frequency/(n_structures*bin_length)
    else:
        density = frequency/(n_structures*bin_length)
        density = density/density.sum()


    return np.linspace(range[0], range[1], bins+1), density

