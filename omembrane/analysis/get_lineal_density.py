import molsysmt as msm
from molsysmt import pyunitwizard as puw
import numpy as np

def get_lineal_density(molecular_system, selection, structure_index=0,
        axis=[0,0,1], bins=100, normalized=False, frequency=False,
        weight=None, syntax='MolSysMT'):

    length_unit = puw.get_standard_units(dimensionality={'[L]':1})

    atom_indices = msm.select(molecular_system, selection=selection, syntax=syntax)
    coordinates = msm.get(molecular_system, selection=atom_indices, structure_indices=structure_index, coordinates=True)
    coordinates = puw.convert(coordinates, to_unit=length_unit)
    box = puw.convert(box, to_unit=length_unit)


    # Hay que incluir aquí el wrap

    box = puw.get_value(box)
    coordinates = puw.get_value(coordinates)

    min_q = np.inf
    max_q = -np.inf

    for ii in [0,1]:
        for jj in [0,1]:
            for kk in [0,1]:

                vertex = ii*box[structure_index][0]+jj*box[structure_index][1]+kk*box[structure_index][2]
                proj = np.dot(vertex, axis)

                if proj<min_q:
                    min_q = proj
                elif proj>max_q:
                    max_q = proj

    bin_length = (max_q-min_q)/bins

    freq = np.zeros(bins)

    coordinates = msm.get(molecular_system, selection=atom_indices, structure_indices=structure_index, coordinates=True)
    coordinates = puw.get_value(coordinates, to_unit=length_unit)

    for aux_coordinates in coordinates[0,:]:
        proj = np.dot(aux_coordinates, axis)
        bin_index = int(np.floor((proj - min_q)/bin_length))
        if bin_index==bins:
            bin_index=bins-1
        if bin_index>bins:
            print(aux_coordinates)
            print(box)
        freq[bin_index]+=1

    if frequency:
        density = freq
    elif not normalized:
        density = freq/bin_length
    else:
        density = freq/bin_length
        density = density/density.sum()

    return np.linspace(min_q, max_q, bins+1), density

