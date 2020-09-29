import numpy as np

def build_leaflet(length_x, length_y, radii):

    # Plano definido en x en [0, length_x] e y en [0, length_y]
    # z para las esferas será de momento 0.0

    n_spheres = len(radii)

    coordinates = np.zeros([n_spheres, 3])

    random_generator = np.random.default_rng()

    for sphere_index in range(n_spheres):

        x = random_generator.uniform(0.0, length_x, 1)
        y = random_generator.uniform(0.0, length_y, 1)

        coordinates[sphere_index, 0] = x
        coordinates[sphere_index, 1] = y

    return coordinates

