import Bio.PDB as pdb
import math
import numpy as np
import os
import sys

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
        golden_ratio = np.pi * (3 - np.sqrt(5))

        offset = 2/ n

        list_points = []

        # Compute the n envenly spaced points 
        for i in range(0,n):
            y = i * offset - 1 + (offset/2)
            r = np.sqrt(1- y * y)
            phi = i * golden_ratio
            x = np.cos(phi) * r
            z = np.sin(phi) * r

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
        TODO
        """
        if name is None:
            file_name = pdb_path.split("/")[-1].split(".")[0]
            self.name = file_name
        # Handle wrong path
        if not os.path.exists(pdb_path):
            raise FileNotFoundError(f"No PDB file found at: {pdb_path}")
        # Handle wrong format (extension)
        elif pdb_path.split('.')[-1] != "pdb":
            raise FileNotFoundError(f"No PDB file found at: {pdb_path}")

        self.structure, self.chains, self.residues, self.atoms = self._extract_structure(pdb_path, model)

    def _extract_structure(self, pdb_path:str, model:int=0): 
        """
        TODO
        """

        chains = []
        residues = []
        atoms = []

        # Use biopython parser to extract whole structure and its content
        structure = pdb.PDBParser().get_structure(self.name, pdb_path)

        for chain in structure[model]:
            # Store the chain's residues
            chain_res = []

            for amino_acid in chain:
                # Exclude water molecules
                if amino_acid.resname == "HOH": continue
                # Store the residue's atoms
                res_atoms = []

                for atom in amino_acid:
                    # Extract data about the atom
                    elem = atom.element
                    id = atom.id
                    residue = amino_acid.resname
                    coord = list(atom.coord)
                    radius = ATOMIC_RADII[elem]

                    # Create an object of the Atom class
                    atom_object = Atom(id, elem, residue, coord, radius)

                    # Add atoms to the structure list
                    atoms.append(atom_object)
                    # Add atoms to the residue's atoms list
                    res_atoms.append(atom_object)
                
                # Extract data about the residue
                id = amino_acid.id
                type = amino_acid.resname

                # Create an object of the residue class
                residue_object = Residue(id, type, res_atoms)

                # Add residue to the structure list
                residues.append(residue_object)
                # Add residue to the chain class
                chain_res.append(residue_object)
            
            # Extract data about the chain
            id = chain.id

            # Create an object of the Chain class
            chain_object = Chain(id, chain_res)

            # Add chain to the structure list
            chains.append(chain_object)
        
        return structure, chains, residues, atoms

class Chain:
    def __init__(self, id, rescontent:list=None):
        self.id = id
        self.rescontent = rescontent
        
        if rescontent is not None:
            self.area = np.sum([res.area for res in self.rescontent])
            self.accessibility = self.set_accessibility()

    def set_accessibility(self, accessibility:float|int=None):
        if accessibility is None:
            accessibility


class Residue:
    def __init__(self, id, type, atomscontent:list=None):
        self.id = id
        self.type = type
        self.atomscontent = atomscontent

        self.area = np.sum([atom.area for atom in self.atomscontent])

        
class Atom:
    def __init__(self, id, elem, residue, coord, radius):
        """TODO"""
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
            coord2 (list): A list of 3 coordinates with which distance is
                calculated.
        Returns:
            distance (float): Distance between the atom and 
                the given point.
        """
        # Extract the Atom coordinates
        coord1 = np.array(self.coord)
        coord2 = np.array(coord2)
        distance = np.linalg.norm(coord1 - coord2)

        return distance
    

    def connectivity(self, atoms_list:list, threshold:int|float=15):
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
    sys.exit()

        
