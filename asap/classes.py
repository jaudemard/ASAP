import Bio.PDB as pdb
import math
import numpy as np
import os

ATOMIC_RADII = {
        "H": 1.200,
        "HE": 1.400,
        "C": 1.700,
        "N": 1.550,
        "O": 1.520,
        "F": 1.470,
        "NA": 2.270,
        "MG": 1.730,
        "P": 1.800,
        "S": 1.800,
        "CL": 1.750,
        "K": 2.750,
        "CA": 2.310,
        "NI": 1.630,
        "CU": 1.400,
        "ZN": 1.390,
        "SE": 1.900,
        "BR": 1.850,
        "CD": 1.580,
        "I": 1.980,
        "HG": 1.550,
    }

pdb.Residue.Residue


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
            x,y,z = self.lattice[i].coord

            # Scale to match radius
            x = x * radius
            y = y * radius
            z = z * radius
            # Move around atom coordinates
            x = x + center[0]
            y = y + center[1]
            z = z + center[2]

            list_points.append(Point([x,y,z]))

        self.lattice = list_points

class Point:
    """Composite class of the Sphere. Represent a point of its lattice."""
    def __init__(self, coord:list):
        self.coord = coord
        self.accessibility = True

    def set_accessibility(self, accessibility:bool=True):
        self.accessibility = accessibility

class Protein:
    def __init__(self, pdb_path:str=None, name:str=None, model:int=0):
        """
        
        """
        if name is None:
            file_name = pdb_path.split("/")[-1].split(".")[0]

        self.name = file_name

        self.structure = self._extract_structure(pdb_path)
        self.chains = self._extract_chains()
        self.residues = self._extract_residues()
        self.atoms = self._extract_atoms()

    def _extract_structure(self, pdb_path:str):
        if pdb_path is not None: # Handle empty path
            # Handle wrong path
            if not os.path.exists(pdb_path):
                raise FileNotFoundError(f"No PDB file found at: {pdb_path}")
            # Handle wrong format (extension)
            elif pdb_path.split('.')[-1] != "pdb":
                raise FileNotFoundError(f"No PDB file found at: {pdb_path}") 
        
        # Extract structure from pdb file
        return pdb.PDBParser().get_structure(self.name, pdb_path)

    def _extract_chains(self):
        list_chains = []
        for chain in self.structure:
            list_chains.append(chain)
        return list_chains
    
    def _extract_residues(self):
        self.list_residue = []
        # Extract information about chain, amino acid ans atoms
        for chain in self.chains:
            for amino_acid in chain:
                # Do not include water molecule
                if amino_acid.resname == "HOH": continue
                for atom in amino_acid:
                    elem = atom.element
                    id = atom.id
                    residue = amino_acid.resname
                    coord = list(atom.coord)
                    radius = ATOMIC_RADII[elem]
                    self.list_atoms.append(Atom(id, elem, residue, coord, radius))
        

        
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
        self.area = 4 * np.pi * (ATOMIC_RADII[elem])**2


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
    

    def connectivity(self, atoms_list:list, threshold:int|float=5):
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
                    point.accessibility = False
                else:
                    point.accessibility = True

        accessible_point = [point.accessibility for point in sphere.lattice].count(True)

        area = 4 * np.pi * (atom.radius)**2

        cover_factor = area / 92

        atom.accessibility = accessible_point * cover_factor

    list_access = [a.accessibility for a in PROTEIN.list_atoms]
    
    print(sum(list_access))

        
