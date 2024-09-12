import argparse
import asap.classes as classes
import numpy as np
import csv
import os
import sys


def command_help():
    """
    Print the command help.
    """
    print("asap help")
    print(__doc__)


def shrake_rupley_cli(command_args):
    """shrake_rupley command line interface"""
    # Create the parser
    parser = argparse.ArgumentParser(
        usage="""usage:
            asap shrake_rupley --pdb=FILE --probe=FLOAT [--output=<output>] [--model=<model>] [--point=<point>]

        Options:
            -h --help  Show help.
            --pdb=<pdb>  Path to the PDB structure
            --probe=<probe>  Size of the probe
            --output=<output>  Output directory for the output files
            --model=<model>  Model used if the pdb includes several of them
            --point=<point>  Number of points in the sphere lattice
        """
    )
    # Define the arguments
    parser.add_argument('--pdb', required=True, help='Path to the PDB structure')
    parser.add_argument('--probe', type=float, required=True, help='Size of the probe')
    parser.add_argument('--output', default='.', help='Output directory for the output files')
    parser.add_argument('--model', type=int, default=0, help='Model used if the pdb includes several of them')
    parser.add_argument('--point', type=int, default=100, help='Number of points in the sphere lattice')

    # Parse the arguments
    args = parser.parse_args(command_args)

    # Ensure output directory exists
    if not os.path.exists(args.output):
        raise FileNotFoundError(f"No directory at: {args.output}")

    # Run function
    shrake_rupley(args.pdb, args.model, args.probe, args.point, args.output)


def shrake_rupley(pdb_file, model, probe, sphere_point, output):
    """From a pdb file, get the Solvent Accessible Surface Area.

    Using a Shrake and Rupley based method of a rolling probe,
    the function assess the protein's atoms accessible surface.
    Create a log with total accessibility and per chain accessibility,
    as well as an atomic level accessibility description.

    Args:
        pdb_file (str): path to the PDB file
        model (int): model to use
        probe (float|int): size of the probe
        sphere_point (int): amount of point to create golden ratio sphere
        output (str): output directory
    """
    # Gets the protein structure
    protein = classes.Protein(pdb_path=pdb_file, model=model)

    for atom in protein.atoms:
        # Create a sphere with golden ratio lattice
        sphere = classes.Sphere(n=sphere_point, method="golden")
        sphere.scale_and_move(center=atom.coord, radius=(atom.radius+probe))

        # Compute closest atoms to limit calculation
        atom.connectivity(protein.atoms, threshold=15)

        accessible_point = 0
        for point in sphere.lattice:

            for neighbour in atom.neighbour:
                # Get distance with neighbour
                dist = neighbour.distance(point.coord)
                # Check if the point accessibility by the probe is blocked
                if dist < (neighbour.radius + (probe)):
                    point.set_accessibility(False)
                    break
            # Point is accessible by the probe if no atoms blocks it   
            if point.accessibility == True:
                accessible_point += 1

        area = 4 * np.pi * (atom.radius+probe)**2
        cover_factor = area / sphere_point
        # Update atom's accessibility
        atom.accessibility = accessible_point * cover_factor

    # Propagate accessibility to residues
    for residue in protein.residues:
        residue.update_accessibility()
    # Propagate accessibility to chains
    for chain in protein.chains:
        chain.update_accessibility()
    # Update the protein accessibility
    protein.update_accessibility()

    # Create log
    log_output(output=output, protein=protein, model=model, probe=probe)
    # Create csv with detailed atomic accessibility
    atomic_accessibility(output=output, protein=protein)
    # Create csv with detailed residue accessibility and relative accessibility
    residue_accessibility(output=output, protein=protein)

    
def log_output(output:str, protein:classes.Protein, model:int, probe:float|int):
    """Create Log of the shrake_rupley command.
    
    Args:
        protein (asap.classes.Protein): A protein object
        model (int): Model of the protein
        probe (int or float): Size of the probe used as solvent
    """
    # 
    if not os.path.exists(output):
        raise FileNotFoundError(f"No directory at: {output}")

    # Create log file
    with open(f"{output}/{protein.name}.log", "w+") as log:
        log.write(f"Solvent accessible surface area for {protein.name}, model {model}.\n"
                    f"Size of the probe: {probe}.\n"
                    f"Method: Shrake and Rupley, with golden ratio spiral.\n"
                    f"Exclude water.\n"
                    f"Exclude hetero-atoms.\n"
                    f"{len(protein.chains)} Chain\n"
                    f"{len(protein.residues)} Residues\n"
                    f"{len(protein.atoms)} Atoms\n"
                    f"Accessible surface area by the solvent (Squared Angstrom):"
                    f"Protein: {round(protein.accessibility,4)}\n"
                    f"Protein relative SASA: {round(((protein.accessibility)*100/protein.max_asa),4)}%\n")
        # Get chain data
        for chain in protein.chains:
            log.write(f"Chain {chain.id}: {round(chain.accessibility,4)}\n")
        
        # Redirect to atomic's level file
        log.write(f"See {output}/{protein.name}.csv for atomic accessibility.")


def atomic_accessibility(output, protein):
    """
    Save protein atom accessibility data to a CSV file.

    Columns: Atom ID, Residue Name, Residue Number, Chain Name, Accessibility,
    and Percentage Accessibility relative to total surface area.

    Args:
        protein (Protein): Protein instance with atom accessibility data.
        output (str): Path to save the CSV file (default "protein_accessibility.csv").
    """
    # Handle wrong output path
    if not os.path.exists(output):
        raise FileNotFoundError(f"No directory at: {output}")

    # Creat csv file
    with open(f"{output}/{protein.name}.csv", mode='w+', newline='') as file:
        writer = csv.writer(file)
        
        # Write header row
        writer.writerow(["atom_id", "residue_name", "residue_id", "chain_name", "accessibility", "percent_accessibility"])
        
        # Fetch atomic's level data
        for atom in protein.atoms:
            atom_id = atom.id
            residue_name = atom.residue  # Name of the residue this atom belongs to
            residue_id = next(res.id for res in protein.residues if atom in res.atomscontent)
            chain_name = next(chain.id for chain in protein.chains if any(res.id == residue_id for res in chain.rescontent))
            accessibility = atom.accessibility
            total_area = atom.area  # The total area is just the surface area of the atom itself
            percent_accessibility = (accessibility / total_area) * 100  # Calculate percentage accessibility
            
            # Write data row
            writer.writerow([atom_id, residue_name, residue_id, chain_name, round(accessibility,4), round(percent_accessibility,4)])

def residue_accessibility(output, protein):
    """
    Save protein atom accessibility data to a CSV file.

    Columns: Atom ID, Residue Name, Residue Number, Chain Name, Accessibility,
    and Percentage Accessibility relative to total surface area.

    Args:
        protein (Protein): Protein instance with atom accessibility data.
        output (str): Path to save the CSV file (default "protein_accessibility.csv").
    """
    # Handle wrong output path
    if not os.path.exists(output):
        raise FileNotFoundError(f"No directory at: {output}")

    # Creat csv file
    with open(f"{output}/{protein.name}_relative.csv", mode='w+', newline='') as file:
        writer = csv.writer(file)
        
        # Write header row
        writer.writerow(["residue_name", "residue_id", "chain_name", "accessibility", "relative_accessibility"])
        
        # Fetch atomic's level data
        for residue in protein.residues:
            residue_name = residue.type
            residue_id = residue.id  # Name of the residue this atom belongs to
            chain_name = next(chain.id for chain in protein.chains if any(res.id == residue_id for res in chain.rescontent))
            accessibility = residue.accessibility
            relative_accessibility = (residue.accessibility/residue.max_asa) * 100
            
            # Write data row
            writer.writerow([residue_name, residue_id, chain_name, round(accessibility,4), round(relative_accessibility,4)])
