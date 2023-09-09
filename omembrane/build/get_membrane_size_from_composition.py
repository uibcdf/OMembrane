import numpy as np
from openmm import unit as u

surface_area = {'POPC':68.3*u.angstroms**2}
radius = {'POPC':4.66*u.angstroms}

def get_membrane_size_from_composition(lower_composition, upper_composition):

    upper_area = 0.0 * u.angstroms**2
    lower_area = 0.0 * u.angstroms**2

    for lipid in upper_composition:
        upper_area += upper_composition[lipid]*surface_area[lipid]

    for lipid in lower_composition:
        lower_area += lower_composition[lipid]*surface_area[lipid]

    filter = abs(upper_area-lower_area)<30.0*u.angstroms**2

    if filter:

        area = abs(upper_area+lower_area)/2.0

        clearance_factor = 4.0/np.pi
        clearance_factor = 1.00

        area = area*clearance_factor

        length_x = np.sqrt(area)
        length_y = length_x

        return length_x, length_y

    raise ValueError

