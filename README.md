# ASAP

The **ASAP** (**Acessible Surface Area Project**) propose an easy python based tool to calculate a protein's accessible surface area to a solvant (SASA).

## Description

The `shrake_rupley` command of asap, aims to compute a protein's SASA from a pdb file.

**The method is based on the work of Shrake and Rupley (1971), as it uses a rolling probe process to capture which part of the atoms surface is accessible or not.**

For each atom, a lattice is create in a sphere shape, around the atom center. The sphere radius is the sum of the Van der Waals radius and the radius of the probe.  
For each point of the lattice, if its distance to another atom is lesser than that atom's Van der Waals and probe radius, then the point is not accessible.  
The accessible surface is then proportionnal to the amount of accessible points.

To space envely the points on the sphere, we chose to use tu Fibonacci Spirale also called the Golden Ratio Spirale.

`shrake_rupley` produces 3 output:
* **`protein.log`**: Data about the Protein, its total and relative accessibility
* **`protein.csv`**: Accessibility at an atomic level
* **`protein_relative.csv`**: Accessibility at a residue level, with relative accessibility per residue.

The relative per residue calculation uses the residues accessibility in a Gly-X-Gly conformation, with X the residue, as the references. They are provided from the work of Tien and Meyer (2013).

## Requirements

Required package can be installed through a yaml file, called `asap_env.yml`.

Using conda, run the command below:  
```bash
conda env create -f asap_env.yml
conda activate asap
```

## Installation

The installation is conducted through pip.

Go into the repository ASAP and execute the following command  
```bash
pip install .
```
You may check if it was correctly installed by running `asap --help`, or with the provided toy and basic usage commands.

/!\ If you get a ModuleNotFound when running `asap`

Use a virtual environment.

Start by running to create the virtual environment.
````bash
python -m venv asap
cd asap
````

Then, activate the virtual environment
```bash
source bin/activate
```

Install the requirements, from the requirements.txt in the git.
```bash
python -m pip install -r requirements.txt
```

Clone the git repository from there
```bash
git clone https://github.com/jaudemard/ASAP
```

Now you can install the module and run it
```bash
python -m pip install .
```


## Basic Usage

We'll use the pdb structure [1L2Y](https://www.rcsb.org/structure/1L2Y) as a toy.

Run a simple accessible surface calculation, without a probe:
```bash
asap shrake_rupley --pdb toy/1l2y.pdb --probe 0
```

Use a probe the size of an oxygen Van der Waals area, to mimic a water molecule, on the third model of the file.
```bash
asap shrake_rupley --pdb toy/1l2y.pdb --probe 1.52 --model 2 --output results/
```

## Citations

**For the method of the rolling probe as a way to determine accessibility**  
A. Shrake et J. A. Rupley, « Environment and exposure to solvent of protein atoms. Lysozyme and insulin », J Mol Biol, vol. 79, nᵒ 2, p. 351‑371, sept. 1973, doi: 10.1016/0022-2836(73)90011-9.


**For the algorithm generating envely spaced point on the surface of Atoms**  
E. B. Saff et A. B. J. Kuijlaars, « Distributing many points on a sphere », The Mathematical Intelligencer, vol. 19, nᵒ 1, p. 5‑11, déc. 1997, doi: 10.1007/BF03024331.


**For the max Accessibility Surface Area references**  
M. Z. Tien, A. G. Meyer, D. K. Sydykova, S. J. Spielman, et C. O. Wilke, « Maximum Allowed Solvent Accessibilites of Residues in Proteins », PLoS One, vol. 8, nᵒ 11, p. e80635, nov. 2013, doi: 10.1371/journal.pone.0080635.


