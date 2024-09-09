"""TODO"""

from classes import Protein, Sphere
import numpy as np



PROTEIN = Protein("toy/1a5d.pdb")

counter = 0

print(len(PROTEIN.list_atoms),"****")

for atom in PROTEIN.list_atoms:
    sphere = Sphere(n=92) 
    sphere.scale_and_move(atom.coord, atom.radius)

    atom.connectivity(PROTEIN.list_atoms, 5)


    for point in sphere.lattice:
        for neighbour in atom.neighbour:
            if neighbour.id == atom.id: continue

            dist = neighbour.distance(point.coord)

            if dist < (neighbour.radius + (1.52*2)):
                point.set_accessibility(False)
            else:
                point.set_accessibility(True)

    accessible_point = [point.accessibility for point in sphere.lattice].count(True)

    area = 4 * np.pi * (atom.radius)**2

    cover_factor = area / 92

    atom.accessibility = accessible_point * cover_factor

list_access = [a.accessibility for a in PROTEIN.list_atoms]

print(sum(list_access))