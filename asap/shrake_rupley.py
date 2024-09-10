"""TODO"""

from classes import Protein, Sphere
import numpy as np

sonde = 1.52

PROTEIN = Protein("toy/1a5d.pdb")

counter = 0

print(len(PROTEIN.atoms),"****")



for atom in PROTEIN.atoms:
    sphere = Sphere(n=92) 
    sphere.scale_and_move(atom.coord, (atom.radius+sonde))

    atom.connectivity(PROTEIN.atoms, 15)

    accessible_point = 0
    for point in sphere.lattice:

        for neighbour in atom.neighbour:
            if neighbour.id == atom.id: continue

            dist = neighbour.distance(point.coord)



            if dist < (neighbour.radius + (1.52)):
                point.set_accessibility(False)
                break

        if point.accessibility == True:
            accessible_point += 1

    area = 4 * np.pi * (atom.radius+sonde)**2

    cover_factor = area / 92

    atom.accessibility = accessible_point * cover_factor

    

list_access = [a.accessibility for a in PROTEIN.atoms]



print(sum(list_access))