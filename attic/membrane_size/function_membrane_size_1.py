import numpy as np
import simtk.unit as unit

def get_membrane_size(n_POPC = 0 , n_DPPC = 0 , n_DOPC = 0):

    surface_area_POPC = 68.3 * unit.angstroms**2 
    surface_area_DPPC = 63  * unit.angstroms**2 
    surface_area_DOPC = 69.7 * unit.angstroms**2

    total_area_POPC = n_POPC * surface_area_POPC
    total_area_DPPC = n_DPPC * surface_area_DPPC
    total_area_DOPC = n_DOPC * surface_area_DOPC

    total_area_lipids = total_area_POPC + total_area_DPPC + total_area_DOPC

    x = np.sqrt(total_area_lipids)
    y = np.sqrt(total_area_lipids)

    print ("The surface area of the membrane is", total_area_lipids)

    return x, y 

