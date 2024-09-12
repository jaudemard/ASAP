"""
usage:
    asap shrake_rupley --pdb=FILE --probe INT [--output=DIR] [--model INT] [--point INT]

option:
    -h  --help    Show help.
    --pdb_file    Path to the PDB structure
    --probe    Size of the probe
    --output    Output directory for the output files
    --prot_name    Given name of the protein
    --model    Model used if the pdb include several of them
    --point    Amount of paint created around each atom
"""
import classes
import numpy as np
import docopt
import csv
import os


def command_help():
    """
    Print the command help.
    """
    print(docopt.docopt(__doc__))


def shrake_rupley_cli(command_args):
    """skrake_rupley command line interface"""
    args = docopt.docopt(__doc__, argv=command_args)

    pdb_file = args["--pdb"]
    probe = float(args["--probe"])

    output = args["--output"] or "."
    # Handle unexisting output file to not waist the execution time
    if not os.path.exists(output):
        raise FileNotFoundError(f"No directory at: {output}")
    
    model = int(args["--model"] or 0)
    sphere_point = int(args["--point"] or 100)

    # Run function
    shrake_rupley(pdb_file, model, probe, sphere_point, output)


def shrake_rupley(pdb_file, model, probe, sphere_point, output):
    # Gets the protein structure
    protein = classes.Protein(pdb_path=pdb_file, model=model)

    for atom in protein.atoms:
        # Create a sphere with golden ratio lattice
        sphere = classes.Sphere(n=sphere_point)
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
                    f"Accessible surface area by the solvent (Squared Angstrum):"
                    f"Protein: {round(protein.accessibility,4)}\n")
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


shrake_rupley("toy/1a5d.pdb", model = 0, probe = 1.52, sphere_point = 92, output="res")

