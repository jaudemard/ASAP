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

"""

import classes
import numpy as np
import docopt

def command_help():
    """
    Print the command help.
    """
    print(docopt.docopt(__doc__))

def shrake_rupley_cli(command_args):
    args = docopt.docopt(__doc__, argv=command_args)

    pdb_file = args["--pdb"]
    probe = args["--probe"]

    output = args["--output"] or "."
    model = args["--model"] or 0
    sphere_point = args["--point"] or 100

    shrake_rupley(pdb_file, model, probe, sphere_point, output)

def shrake_rupley(pdb_file, model, probe, sphere_point):
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
                if neighbour.id == atom.id: continue

                # Get distance with neighbour
                dist = neighbour.distance(point.coord)
                # Check if the point accessibility by the probe is blocked
                if dist < (neighbour.radius + (probe)):
                    point.set_accessibility(False)
                    break

            if point.accessibility == True:
                accessible_point += 1

        area = 4 * np.pi * (atom.radius+probe)**2
        cover_factor = area / sphere_point
        atom.accessibility = accessible_point * cover_factor
    protein.accessibility = np.sum([atom.accessibility for atom in protein.atoms])

    for residue in protein.residues:
        residue.update_accessibility()
    
    for chain in protein.chains:
        chain.update_accessibility()

    # Update the protein accessibility
    protein.update_accessibility()

    print(protein.accessibility)
    
    def shrake_rupley_output(output, protein, model, execution_time):

        # Create Log
        with open(f"{output}/{protein.name}.log", "w+") as log:
            log.write(f"Solvent accessible surface area for {protein.name}, model {model}.\n",
                      f"Size of the probe: {probe}.\n",
                      "Method: Shrake and Rupley, with golden ratio spiral.\n",
                      f"Exclude Water Molecule (HOH).\n",
                      f"{len(protein.chains)} Chain\n",
                      f"{len(protein.residues)} Residues\n",
                      f"{len(protein.atoms)} Atoms\n"
                      f"Total accessibility: {protein.accessibility}\n.",
                      f"It took {execution_time}s.")

shrake_rupley("toy/1a5d.pdb", name = "1a5d", model = 0, probe = 1.52, sphere_point = 92)