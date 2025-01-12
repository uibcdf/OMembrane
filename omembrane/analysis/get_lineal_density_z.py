import molsysmt as msm
from molsysmt import pyunitwizard as puw
import numpy as np

def get_lineal_density_z(molecular_system, selection, structure_index=0, bins=100, normalized=False, frequency=False,
        weights=None, box_z_origin='0 nanometers', box_z_center=None, center_of_selection=False,
        center_coordinates=None, syntax='MolSysMT'):

    length_unit = puw.get_standard_units(dimensionality={'[L]':1})

    box = msm.get(molecular_system, structure_indices=structure_index, box=True)
    box = puw.convert(box, to_unit=length_unit)

    box_origin = None
    if box_z_origin is not None:
        value = puw.get_value(box_z_origin, to_unit=length_unit)
        box_origin = puw.quantity([0,0,value], length_unit)

    box_center = None
    if box_z_center is not None:
        value = puw.get_value(box_z_center, to_unit=length_unit)
        box_center = puw.quantity([0,0,value], length_unit)

    atom_indices = msm.select(molecular_system, selection=selection, syntax=syntax)
    coordinates = msm.get(molecular_system, selection=atom_indices, structure_indices=structure_index, coordinates=True)
    coordinates = puw.convert(coordinates, to_unit=length_unit)

    wrapped_dict = msm.pbc.wrap_to_pbc({'coordinates':coordinates, 'box':box}, box_origin=box_origin, box_center=box_center)
    coordinates = wrapped_dict['coordinates']

    box = puw.get_value(box)
    coordinates = puw.get_value(coordinates)

    min_q = np.inf
    max_q = -np.inf

    for ii in [0,1]:
        for jj in [0,1]:
            for kk in [0,1]:

                vertex = ii*box[structure_index][0]+jj*box[structure_index][1]+kk*box[structure_index][2]
                proj = vertex[2]

                if proj<min_q:
                    min_q = proj
                elif proj>max_q:
                    max_q = proj

    bin_length = (max_q-min_q)/bins

    freq = np.zeros(bins)

    for aux_coordinates in coordinates[0,:]:
        proj = aux_coordinates[2]
        bin_index = int(np.floor((proj - min_q)/bin_length))
        freq[bin_index]+=1

    if frequency:
        density = freq
    elif not normalized:
        density = freq/bin_length
    else:
        density = freq/bin_length
        density = density/density.sum()

    return np.linspace(min_q, max_q, bins+1), density

