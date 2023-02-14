import numpy as np
import simtk.unit as unit

def get_n_lipids_from_size(total_area = 10000 * unit.angstroms**2, p_POPC = 0.5 , p_DPPC = 0.3, p_DOPC = 0.2):


    area_POPC = (4 / np.pi)* 68.3 * unit.angstroms**2
    area_DPPC = (4 / np.pi)* 63  * unit.angstroms**2
    area_DOPC = (4 / np.pi)* 69.7 * unit.angstroms**2

    r_POPC = p_POPC / np.min([p_POPC, p_DPPC, p_DOPC])
    r_DPPC = p_DPPC / np.min([p_POPC, p_DPPC, p_DOPC])
    r_DOPC = p_DOPC / np.min([p_POPC, p_DPPC, p_DOPC])


    area_unidad = r_POPC * area_POPC + r_DPPC * area_DPPC + r_DOPC * area_DOPC
    
    no_unidades_totales = total_area / area_unidad
    print (no_unidades_totales)

    n_POPC = int(np.round(no_unidades_totales *r_POPC))
    n_DPPC = int(np.round(no_unidades_totales *r_DPPC))
    n_DOPC = int(np.round(no_unidades_totales *r_DOPC))

    n_total_lipids = n_POPC + n_DPPC + n_DOPC 

    print ("The n_POPC in the membrane is", n_POPC, "(", n_POPC / n_total_lipids,")")
    print ("The n_DPPC in the membrane is", n_DPPC, "(", n_DPPC / n_total_lipids,")")
    print ("The n_DOPC in the membrane is", n_DOPC, "(", n_DOPC / n_total_lipids,")")


    print ("The result area is", n_POPC * area_POPC + n_DPPC * area_DPPC + n_DOPC * area_DOPC)
