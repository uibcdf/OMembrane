"""
Unit and regression test for the analysis.get_lineal_density function of the openmembrane package.
"""

# Import package, test suite, and other packages as needed
import openmembrane as omem
import molsysmt as msm
import numpy as np

def test_get_lineal_density_1():
    molsys = msm.convert(omem.demo['POPC']['POPC.msmpk'])
    bins_edges, density = omem.analysis.get_lineal_density(molsys, selection='atom_name=="P"', axis=[0,0,1], range=[0.0, 10.0],
                                                       bins=20, normalized=False)
    good_density = np.zeros([20])
    good_density[2] = 0.92
    good_density[3] = 77.88
    good_density[4] = 196.16
    good_density[5] = 19.04
    good_density[10] = 1.24
    good_density[11] = 88.88
    good_density[12] = 184.16
    good_density[13] = 19.72
    assert bins_edges.shape==(21,)
    assert np.allclose(density, good_density)

def test_get_lineal_density_2():
    molsys = msm.convert(omem.demo['POPC']['POPC.msmpk'])
    bins_edges, density = omem.analysis.get_lineal_density(molsys, selection='atom_name=="P"', axis=[0,0,1],
                                                       bins=20, normalized=True)
    good_density = np.array([0.00360544, 0.04727891, 0.19156463, 0.1947619 , 0.05965986,
       0.00312925, 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.00020408, 0.00414966,
       0.05666667, 0.18367347, 0.18986395, 0.06122449, 0.00421769])
    assert bins_edges.shape==(21,)
    assert np.isclose(bins_edges[0], 1.28126097)
    assert np.isclose(bins_edges[20], 6.982036113739014)
    assert np.allclose(density, good_density)


