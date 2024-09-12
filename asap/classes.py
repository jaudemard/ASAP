import Bio.PDB as pdb
import numpy as np
import os
from scipy.spatial import KDTree

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
    def __init__(self, center=np.array([0, 0, 0]), radius=1, method="golden", n:int=92):
        """
        A sphere is a 3-dimensional object made of a center point and 
        a lattice. The lattice is done using the golden ratio method to space
        evenly points around the center.
        """
        self.center = center
        self.radius = radius
        self.n = n

        if method == "golden":
            self.lattice = self._golden_spirale()
        elif method == "saff_kuijlaars":
            self.lattice = self._saaf_kuijlaars()
        else:
            raise ValueError("Unknown method to generate spiral.")

    def _golden_spirale(self):
        """
        Generate sphere forming points using the golden spiral or fibonacci method.
        """
        golden_ratio = np.pi * (3 - np.sqrt(5))
        offset = 2 / n
        points = []

        for i in range(self.n):
            y = i * offset - 1 + offset / 2
            r = np.sqrt(1 - y ** 2)
            phi = i * golden_ratio
            x = np.cos(phi) * r
            z = np.sin(phi) * r

            points.append(Point(np.array([x, y, z])))

        return points

    def _saff_kuijlaars(self):
        """
        Generate sphere forming points usin the Saaf and Kuiklaars method.
        """
        points = [] # Store points of the lattice
        for k in range(self.n):
            theta = np.arccos(1 - 2 * (k + 0.5) / self.n)
            phi = np.pi * (1 + np.sqrt(5)) * k
            
            # Get coordinates
            x = np.sin(theta) * np.cos(phi) 
            y = np.sin(theta) * np.sin(phi)
            z = np.cos(theta)
            
            points.append(Point(np.array([x, y, z])))
        
        return points

    def scale_and_move(self, center, radius):
        """Scale up or down the sphere and move it to a new center."""
        radius = radius or self.radius
        new_lattice = [
            Point(point.coord * radius + center) for point in self.lattice
        ]
        self.lattice = new_lattice

class Point:
    def __init__(self, coord):
        """
        A point representing a coordinate in 3D space.

        Args:
            coord (list or np.array): 3D coordinates (x, y, z).
        """
        self.coord = np.array(coord)
        self.accessibility = True  # Default is that the point is accessible
    
    def set_accessibility(self, accessibility=True):
        """
        Update the accessibility of the point.

        Args:
            accessibility (bool): True if accessible, False if blocked.
        """
        self.accessibility = accessibility

class Protein:
    def __init__(self, pdb_path, model=0):
        # Handle wrong path
        if not os.path.exists(pdb_path):
            raise FileNotFoundError(f"No PDB file found at: {pdb_path}")
        # Handle wrong extension
        elif pdb_path.split(".")[-1] != "pdb":
            raise ValueError(f"Invalid file format: {pdb_path}")
        
        self.name = os.path.basename(pdb_path).split(".")[0]
        self.atoms, self.residues, self.chains, self.structure = self._extract_structure(pdb_path, model)
        self.area = np.sum([chain.area for chain in self.chains])
        self.accessibility = self.area

    def _extract_structure(self, pdb_path, model=0):
        """
        Extract structure from a PDB file using BioPython parser.

        Args:
            pdb_path (str): Path to the pdb_file
            model (int): Model to use if there is several in the pdb.
        Returns:
            atoms (list of Atom): Atoms list in the protein
            residues (list of Residue): Residues list in the protein
            chains (list of Chain): Chains list in the protien
            structure: The BioPython structure object builded from the file.
        """
        atoms, residues, chains = [], [], []

        # Extract the structure with Biopython
        structure = pdb.PDBParser().get_structure(self.name, pdb_path)

        # Parse of the structure
        for chain in structure[model]:
            chain_res = []
            for amino_acid in chain:
                # Exclude Water
                if amino_acid.resname == "HOH": continue
                # Exclue Hetero Atoms
                if amino_acid.id[0] == "H": continue
                res_atoms = [] # to store the residue atoms
                for atom in amino_acid:
                    coord = np.array(atom.coord)
                    elem = atom.element
                    radius = ATOMIC_RADII[elem]
                    atom_obj = Atom(atom.id, elem, coord, radius, amino_acid.resname)

                    # Create an atom
                    atoms.append(atom_obj) # Add to the protein's atoms list
                    res_atoms.append(atom_obj) # Add to the residue's atoms list
                
                # Create a Residue
                residue_obj = Residue(amino_acid.id[1], amino_acid.resname, chain.id, res_atoms)
                residues.append(residue_obj) # Add to the protein's residues list
                chain_res.append(residue_obj) # Add to the chain's residues list

            # Create a Chain
            chain_obj = Chain(chain.id, chain_res, self.name)
            chains.append(chain_obj) # Add to the protein's chains list

        return atoms, residues, chains, structure

    def update_accessibility(self):
        """
        Update accessibility with the sum of the protein's chains accessibility.
        """
        self.accessibility = np.sum([chain.accessibility for chain in self.chains])


class Chain:
    def __init__(self, id, rescontent=None, protein=None):
        """A chain is part of a Protein"""
        self.id = id
        self.rescontent = rescontent or []
        self.protein = protein or "protein"

        # Surface area is the sum its residues area
        self.area = np.sum([res.area for res in self.rescontent])
        self.accessibility = self.area

    def update_accessibility(self):
        """
        Update accessibility with the sum of the chain's residues accessibility.
        """
        self.accessibility = np.sum([res.accessibility for res in self.rescontent])


class Residue:
    def __init__(self, id, type, chain=None, atomscontent=None):
        """A Residue is part of a Chain's"""
        self.id = id
        self.type = type
        self.atomscontent = atomscontent or []
        self.chain = chain or "Chain"
        
        # Surface area is the sum of its atoms Van Der Walls surface
        self.area = np.sum([atom.area for atom in self.atomscontent])
        self.accessibility = self.area

    def update_accessibility(self):
        """
        Update the accessibility with the sum of the residue's atoms accessibility.
        """
        self.accessibility = np.sum([atom.accessibility for atom in self.atomscontent])


class Atom:
    def __init__(self, id:str, elem:str, coord:np.array, radius:int|float, residue:str=None):
        """An Atom is part of a Residue."""
        self.id = id
        self.elem = elem
        self.coord = coord
        self.residue = residue
        self.radius = radius # Radius of the Van der Walls space around the atom
        self.area = 4 * np.pi * (ATOMIC_RADII[elem]) ** 2
        self.accessibility = self.area
        self.neighbour = [] # Closest atoms

    def distance(self, coord2:np.array):
        """Compute distance between an atom and a point.
        Args:
            coord2 (np.array): Coordinates of a point in a 3D space.
        Returns:
            Distance between self coordinates and the point.
        """
        return np.linalg.norm(self.coord - coord2)

    def connectivity(self, atoms_list:list, threshold:int|float):
        """Find nearest atoms using KDTree.

        Args:
            atoms_list (list of Atom): A list of Atoms within which
                neighbours are searched for.
            threshold (int or float): Maximum distance to qualify an atom as
                a neighbour.
        """
        atom_coords = np.array([atom.coord for atom in atoms_list])
        # Build a KDTree for the atom coordinates
        tree = KDTree(atom_coords)
        # Find indices of atoms within the threshold distance
        indices = tree.query_ball_point(self.coord, r=threshold)
        # Use the indices to find neighboring atoms
        self.neighbour = [atoms_list[i] for i in indices if not np.array_equal(atoms_list[i].coord, self.coord)]


if __name__ == "__main__":
    PROTEIN = Protein("toy/1l2y.pdb")
    print(PROTEIN.chains[0].id)
        
