import molsysmt as msm
from molsysmt import pyunitwizard as puw
import numpy as np

def get_lineal_density(molecular_system, selection, axis=[0,0,1], bins=100, bin_length=None, range=None,
        normalized=False, frequency=False, weight=None, syntax='MolSysMT'):

    atom_indices = msm.select(molecular_system, selection=selection, syntax=syntax)

    n_structures = msm.get(molecular_system, n_structures=True)

    length_unit = puw.get_standard_units(dimensionality={'[L]':1})

    if range is None:

        iterator = msm.Iterator(molecular_system, selection = atom_indices, coordinates = True)

        min_total = puw.quantity(np.inf,length_unit)
        max_total = puw.quantity(-np.inf,length_unit)

        for aux_coordinates in iterator:
            proj = np.dot(aux_coordinates, axis)
            min_coordinate = proj.min()
            max_coordinate = proj.max()

            if min_coordinate < min_total:
                min_total = min_coordinate

            if max_coordinate > max_total:
                max_total = max_coordinate

        range = [min_total, max_total]

    else:

        range=puw.quantity(range)


    if bin_length is None:
        bin_length = (range[1]-range[0])/bins
    else:
        bin_length=puw.quantity(bin_length)
        bins = (range[1]-range[0])/bin_length
        bins = int(np.ceil(bins))
        range[1] = range[0]+bins*bin_length

    freq = np.zeros(bins)

    iterator = msm.Iterator(molecular_system, selection = atom_indices, coordinates = True)

    bin_length_value = puw.get_value(bin_length)
    range_0_value = puw.get_value(range[0])

    for aux_coordinates in iterator:
        aux_coordinates = puw.get_value(aux_coordinates)
        proj = np.dot(aux_coordinates, axis)
        proj = proj.flatten()

        for ii in proj:
            bin_index = int(np.floor((ii - range_0_value)/bin_length_value))
            if bin_index==bins:
                bin_index=bins-1
            freq[bin_index]+=1

    if frequency:
        density = freq
    elif not normalized:
        density = freq/(n_structures*bin_length)
    else:
        density = freq/(n_structures*bin_length)
        density = density/density.sum()

    return np.linspace(range[0], range[1], bins+1), density

