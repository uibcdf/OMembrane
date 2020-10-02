import numpy as np

def build_leaflet(length_x, length_y, radii):

    area_disponible = length_x * length_y
    area_spheres = 0
    for radius in radii:
        area_spheres += np.pi*radius**2

    if area_spheres>area_disponible:
        print("Caja demasiado pequeña")
        return

    # Plano definido en x en [0, length_x] e y en [0, length_y]
    # z para las esferas será de momento 0.0

    n_spheres = len(radii)

    coordinates = np.zeros([n_spheres, 3])

    random_generator = np.random.default_rng()

    for sphere_index in range(n_spheres):

        entra = False

        while entra==False:

            x = random_generator.uniform(0.0, length_x, 1)
            y = random_generator.uniform(0.0, length_y, 1)

            entra = True

            for ii in range(sphere_index):

                dist = np.sqrt( (x-coordinates[ii,0])**2 + (y-coordinates[ii,1])**2 )
                if dist < radii[ii] + radii[sphere_index]:
                    entra = False
                    break

        coordinates[sphere_index, 0] = x
        coordinates[sphere_index, 1] = y

    return coordinates

