import Bio.PDB as pdb
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


class Sphere:
    def __init__(self, center:list=[0,0,0], radius:int|float=1, method:str="golden", n:int=92):
        """
        A sphere is a 3-dimensional object made of a center point and 
        a lattice. The lattice is done using the golden ratio method, to space
        evenly point around the center.
        """
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

            coord = np.array([x,y,z])

            list_points.append(Point(coord))
                
        return list_points
    

    def scale_and_move(self, center:np.array, radius:int|float):
        """
        Scale up or down the sphere and move it to a new center.

        Args:
            center (np.array): An array of 3 coordinates
            raidus (int or float): New radius to scale the sphere        
        """
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

            coord = np.array([x,y,z])

            list_points.append(Point(coord))

        self.lattice = list_points

class Point:
    def __init__(self, coord:list):
        """Composite class of the Sphere. Represent a point of its lattice."""
        self.coord = coord
        self.accessibility = True # Accessibility of the point by the probe

    def set_accessibility(self, accessibility:bool=True):
        """Update accessibility of the point."""
        self.accessibility = accessibility


class Protein:
    def __init__(self, pdb_path:str=None, model:int=0):
        """
        Protein is extracted from a pdb file, include its atoms,
        residues and chains.
        """
        # Handle wrong path
        if not os.path.exists(pdb_path):
            raise FileNotFoundError(f"No PDB file found at: {pdb_path}")
        # Handle wrong format (extension)
        elif pdb_path.split('.')[-1] != "pdb":
            raise FileNotFoundError(f"No PDB file found at: {pdb_path}")
        
        # Get protein name
        file_name = pdb_path.split("/")[-1].split(".")[0]
        self.name = file_name
        # Get structure instances
        self.atoms, self.residues, self.chains, self.structure = self._extract_structure(pdb_path, model)
        # Area is the sum of the chains area
        self.area = np.sum([chain.area for chain in self.chains])
        # Default accessibility is the area
        self.accessibility = self.area

    def _extract_structure(self, pdb_path:str, model:int=0): 
        """
        Extract a structure from a pdb file.

        Parse the pdb file with biopython PDB parser, to get chains, 
        residues and the atoms. Create instance of Chain, Residue and
        Atom class for each of them. Return list of Chain, Residue,
        Atom and the protein's structure.

        It excludes water molecules.

        Args:
            pdb_path (str): Path of the pdb file.
            model (int): number of model to use if different models are
                in the pdb file.
        Returns:
            atoms (list of Atom): A list of Atom, from the protein's atoms.
            residues (list of Residue): A list of Residue from the protein's residues.
            chains (list of Chain): A list of Chain from the protein's chains.
            structure: The Bio.PDB.PDBParser structure object
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
                    coord = np.array(atom.coord)
                    radius = ATOMIC_RADII[elem]

                    # Create an object of the Atom class
                    atom_object = Atom(id, elem, residue, coord, radius)
                    # Add atoms to the structure list
                    atoms.append(atom_object)
                    # Add atoms to the residue's atoms list
                    res_atoms.append(atom_object)
                
                # Extract data about the residue
                id = amino_acid.id[1]
                type = amino_acid.resname
                res_chain = chain.id

                # Create an object of the residue class
                residue_object = Residue(id, type, res_chain, res_atoms)

                # Add residue to the structure list
                residues.append(residue_object)
                # Add residue to the chain class
                chain_res.append(residue_object)
            
            # Extract data about the chain
            id = chain.id
            protein = self.name

            # Create an object of the Chain class
            chain_object = Chain(id, chain_res, protein)

            # Add chain to the structure list
            chains.append(chain_object)
        
        return atoms, residues, chains, structure
    
    
    def update_accessibility(self):
        """
        Update accessibility with the sum of the chain accessibility.
        """
        self.accessibility = np.sum([chain.accessibility for chain in self.chains])

class Chain:
    def __init__(self, id, rescontent:list=None, protein:str=None):
        """
        Chain is part of a protein.
        """
        self.id = id
        self.rescontent = rescontent or []
        self.protein = protein or "protein"

        # Area of a chain is the sum of the residues area
        self.area = np.sum([res.area for res in rescontent])
        # Default accessibility is the area
        self.accessibility = self.area # Default is full accessibility
        
    def update_accessibility(self):
        """
        Update accessibility with the sum of the residues accessibility.
        """
        # Chain accessibility is the sum of its residue accessibility
        self.accessibility = np.sum([res.accessibility for res in self.rescontent])
            

class Residue:
    def __init__(self, id:str, type:str, chain:Chain=None, atomscontent:list=None):
        """
        Residue is part of a protein chain, or a standalone amino acid.
        """
        self.id = id
        self.type = type
        self.atomscontent = atomscontent or []
        self.chain = chain or "Chain"
        
        # Area of a residue is the sum of the atoms area
        self.area = np.sum([atom.area for atom in self.atomscontent])
        # Default accessibility is the area
        self.accessibility = self.area # Default is full accessibility
    
    def update_accessibility(self):
        """
        Update accessibility with the sum of the atoms accessibility.
        """
        # Amino acid accessibility is the sum of its atoms accessibility
        self.accessibility = np.sum([atom.accessibility for atom in self.atomscontent])
        
class Atom:
    def __init__(self, id:str, elem:str, coord:np.array, radius:float|int, residue:str=None):
        """
        An Atom is part of a residue, or a standalone atom.
        """
        self.id = id
        self.elem = elem
        self.residue = residue # Name of the parent residue
        self.coord = coord
        self.radius = radius # Van der Walls raidus

        # Van der Walls surface
        self.area = 4 * np.pi * (ATOMIC_RADII[elem])**2
    
        self.accessibility = self.area # Default is all surface is accessible

        # Closest atoms
        self.neighbour = []


    def distance(self, coord2:np.array):
        """
        Distance calculation between the atom and a given point.
        
        Args
            coord2 (np.array): A list of 3 coordinates with which distance is
                calculated.
        Returns:
            distance (float): Distance between the atom and 
                the given point.
        """
        # Extract the Atom coordinates
        coord1 = self.coord
        coord2 = coord2
        distance = np.linalg.norm(coord1 - coord2)
        return distance
    

    def connectivity(self, atoms_list:list, threshold:int|float=None, probe_size:int|float=None):
        """Compute the neighbour of the Atom with a max distance.
        
        Args
            atoms_list (list): list of Atom object with which the connectivity
                is computed.
            threshold (float or int): maximum distance to characterize a
                neighbouring atom.
            probe_size (float or int): Size of the probe, to compute ideal threshold automatically.
        
        Returns
            neighbours (list): list of Atoms considered as neighbours
        """
        if threshold is None:
            pass # TODO
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
    print(PROTEIN.chains[0].id)
        
