import Bio.PDB as pdb
import math
import numpy as np
import os


class Sphere:
    def __init__(self, center:list=[0,0,0], radius:int|float=1, method:str="golden", n:int=92):
        self.center = center
        self.radius = radius

        if method == "golden":
            self.lattice = self._golden_spirale(n)
        else:
            raise Exception("Unkown method to generate spirale.")


    def _golden_spirale(self, n:int=92):
        """
        Generate sphere forming points around the atom center.

        Uses the golden spiral method to creates evenly spaces points
        in the shape of a 3 dimensional sphere, like a sphere shaped lattice.

        Args:
            n (int): Amount of points forming the sphere.

        Returns:
            list_points (list of tuples): Coordinates (x,y,z) of the points in the sphere.
        """
        golden_ratio = (1 + 5**0.5)*0.5

        list_points = []

        # Compute the n envenly spaced points 
        for i in range(0,n):
            theta = 2 * np.pi * i / golden_ratio
            phi = np.arccos(1 - 2*(i+0.5)/n)

            x = np.cos(theta) * np.sin(phi)
            y = np.sin(theta) * np.sin(phi)
            z = np.cos(phi)

            list_points.append(Point([x,y,z]))
                
        return list_points
    

    def scale_and_move(self, center:list=[0,0,0], radius:int|float=None):
        # Use Van der Walls surface radius by default
        if radius is None:
            radius = self.radius
        # Extract amount of points generated
        n = len(self.lattice)

        list_points = []
        for i in range(0,n):
            x,y,z = self.lattice[i]

            # Scale to match radius
            x = x * radius
            y = y * radius
            z = z * radius
            # Move around atom coordinates
            x = x + center[0]
            y = y + center[1]
            z = z + center[2]

            list_points.append(Point([x,y,z]))

        return list_points

class Point:
    """Composite class of the Sphere. Represent a point of its lattice."""
    def __init__(self, coord:list):
        self.coord = coord

    def set_accessiblity(self, accesibility:bool=True):
        self.accessibility = accesibility

class Protein:
    def __init__(self, pdb_path:str, name:str="protein"):
        self.name = name
        
        if pdb_path is not None: # Handle empty path
            if not os.path.exists(pdb_path):
                raise FileNotFoundError(f"No PDB file found at: {pdb_path}")   
        
        # Extract structure from pdb file
        self.structure = pdb.PDBParser().get_structure(name, pdb_path)

        self.list_atoms = []
        # Store the protein atoms as Atom object in a list
        for model in self.structure:
            for chain in model:
                for amino_acid in chain:
                    for atom in amino_acid:
                        elem = atom.element
                        id = atom.id
                        residue = amino_acid
                        coord = list(atom.coord)
                        self.list_atoms.append(Atom(id, elem, residue, coord, 2))

        
class Atom:
    def __init__(self, id, elem, residue, coord, radius):
        '''
        Atoms are part of proteins
        '''
        self.id = id
        self.elem = elem
        self.residue = residue
        self.coord = coord
        self.radius = radius


    def distance(self, coord2:list):
        """
        Distance calculation between the atom and a given point.
        
        Args
            coord (list): A list of 3 coordinates with which distance is
                calculated.
        
        Returns:
            distance (float): Distance between the atom and 
                the given point.
        """
        # Extract the Atom coordinates
        coord1 = self.coord
        # Compute distance
        sum = (coord2[0] - coord1[0])**2 + (coord2[1] - coord2[1])**2 + (coord2[2] - coord1[2])**2
        distance = math.sqrt(sum)

        return distance
    

    def connectiviy(self, atoms_list:list, threshold:int|float=5):
        """Compute the neighbour of the Atom with a max distance.
        
        Args
            atoms_list (list): list of Atom object with which the connectivity
                is computed.
            threshold (float or int): maximum distance to characterize a
                neighbouring atom.
        
        Returns
            neighbours (list): list of Atoms considered as neighbours
        """
        # Create (or reset) attribute to store the neighbour
        self.neighbour = []
        # Check distance with atoms in the list
        for atom in atoms_list:
            distance = self.distance(atom.coord)
            # If distance is less than the threshold it's a neighbour
            if distance < threshold:
                self.neighbour.append(atom)


if __name__ == "__main__":

    PROTEIN = Protein("toy/1l2y.pdb")

    atome = PROTEIN.list_atoms[75]

    sphere = Sphere(n=50)

    counter = 0

    print(type(sphere.scale_and_move(atome.coord, atome.radius)))

""" 
sphere = Sphere(n=92) 
for atom in protein.list_atoms:
    for s


"""

