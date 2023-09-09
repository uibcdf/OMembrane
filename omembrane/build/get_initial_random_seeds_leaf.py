import openmm as mm
from openmm import unit as u
import numpy as np
import math
from tqdm import tqdm

rng = np.random.default_rng()


def _get_max_force(f):
    fmax = 0.0*f.unit**2
    for aa in f:
        mod2 = aa[0]**2+aa[1]**2+aa[2]**2
        if fmax<mod2:
            fmax=mod2
    return np.sqrt(fmax)


def _distance_pbc(coors1, coors2, box):

    x = coors2[0] - coors1[0]
    y = coors2[1] - coors1[1]

    Lx = box[0]
    Ly = box[1]

    if x>= 0.5*Lx:
        x+= -Lx
    elif x< -0.5*Lx:
        x+= Lx

    if y>= 0.5*Ly:
        y+= -Ly
    elif x< -0.5*Ly:
        y+= Ly

    d = math.sqrt(x**2+y**2)

    return d


def _pbc_wrap(positions, length_x, length_y):
    n = positions.shape[0]
    for ii in range(n):
        x = positions[ii,0]
        lx = x/length_x
        if (lx>=1.0) or (lx<0.0):
            positions[ii,0]=x-np.floor(lx)*length_x
        y = positions[ii,1]
        ly = y/length_y
        if (ly>=1.0) or (ly<0.0):
            positions[ii,1]=y-np.floor(ly)*length_y
    return positions


def _initial_positions(radii, length_x, length_y):

    length_x = length_x._value
    length_y = length_y._value

    n_discs = len(radii)

    coordinates = np.zeros((n_discs,2))

    for ii in range(n_discs):

        is_good = False
        attempts = 0

        while not is_good:

            x = rng.uniform(low=0.0, high=length_x)
            y = rng.uniform(low=0.0, high=length_y)

            is_good = True

            for jj in range(ii):

                dd = _distance_pbc([x,y], coordinates[jj], box=[length_x, length_y])

                if dd <= radii[ii]+radii[jj]:

                    is_good = False
                    break

            attempts += 1

            if attempts >= 50000:
                raise ValueError


        coordinates[ii,:] = [x,y]

    return coordinates


def get_initial_random_seeds_leaf(lipid_radius, length_x, length_y):

    factor = 0.7

    n_discs = len(lipid_radius)
    initial_disc_radius = [factor*radius._value for radius in lipid_radius]

    coordinates = _initial_positions(initial_disc_radius, length_x, length_y)

    discs_charge=[]*u.elementary_charge
    discs_sigma=[]*u.angstrom
    discs_epsilon=[]*u.kilocalories_per_mole

    for radius in lipid_radius:
        sigma = 0.89*(2*radius)
        epsilon = 0.1*0.238*u.kilocalories_per_mole # argon
        charge = 0.0*u.elementary_charge
        discs_charge.append(charge)
        discs_sigma.append(sigma)
        discs_epsilon.append(epsilon)

    discs_positions = np.zeros([n_discs,3], dtype=float)*u.angstroms
    discs_positions[:,0] = coordinates[:,0]*u.angstroms
    discs_positions[:,1] = coordinates[:,1]*u.angstroms

    # System
    system = mm.System()

    # Periodic box vectors.
    a = [length_x._value, 0.0, 0.0]*u.angstrom
    b = [0.0, length_y._value, 0.0]*u.angstrom
    c = [0.0, 0.0, 100*max(lipid_radius)._value]*u.angstrom

    system.setDefaultPeriodicBoxVectors(a, b, c)

    nb = mm.NonbondedForce()
    nb.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
    nb.setCutoffDistance(3.0*max(discs_sigma))
    nb.setUseDispersionCorrection(True)

    for charge, sigma, epsilon in zip(discs_charge, discs_sigma, discs_epsilon):
        system.addParticle(40.0*u.amu)
        nb.addParticle(charge, sigma, epsilon)

    _ = system.addForce(nb)

    # Add a restrining potential to keep atoms in z=0
    energy_expression = 'k * (z^2)'
    force = mm.CustomExternalForce(energy_expression)
    force.addGlobalParameter('k', 1000)
    for disc_index in range(n_discs):
        force.addParticle(disc_index, [])
    _ = system.addForce(force)

    # Definición del estado termodinámico y el integrador.

    step_size = 0.1*u.femtoseconds
    temperature = 100*u.kelvin
    friction = 10.0/u.picosecond # Damping para la dinámica de Langevin

    integrator = mm.LangevinIntegrator(temperature, friction, step_size)

    # Creación de la plataforma.

    platform_name = 'CUDA'
    platform = mm.Platform.getPlatformByName(platform_name)

    # Creación del contexto.

    context = mm.Context(system, integrator, platform)

    # Condiciones iniciales

    initial_positions  = discs_positions

    initial_velocities = np.zeros([n_discs, 3], np.float32) * u.angstroms/u.picoseconds

    context.setPositions(initial_positions)
    context.setVelocities(initial_velocities)

    f = context.getSystem().getForce(0)

    for aux_factor in tqdm(np.linspace(factor, 1.0, num=1000, endpoint=True)):

        for index, charge, sigma, epsilon in zip(range(n_discs), discs_charge, discs_sigma, discs_epsilon):
            f.setParticleParameters(index, charge, aux_factor*sigma, epsilon)

        f.updateParametersInContext(context)
        state = context.getState(getEnergy=True, getForces=True)
        Ep = state.getPotentialEnergy()
        max_force = _get_max_force(state.getForces())

        #print('>>>', aux_factor, Ep, max_force)

        if max_force>100.0*max_force.unit:

            mm.LocalEnergyMinimizer.minimize(context)

            context.getIntegrator().step(50000)

            mm.LocalEnergyMinimizer.minimize(context)
            state = context.getState(getEnergy=True, getForces=True)
            Ep = state.getPotentialEnergy()
            max_force = _get_max_force(state.getForces())

            #print('  >', aux_factor, Ep, max_force)

    positions = context.getState(getPositions=True).getPositions(asNumpy=True)

    positions = _pbc_wrap(positions, length_x, length_y)

    return positions

