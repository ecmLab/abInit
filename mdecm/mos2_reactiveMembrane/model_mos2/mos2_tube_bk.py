import numpy as np
import os
from ase.build import mx2
from ase import Atoms
from ase.io import write
from Bio.PDB import PDBParser, PDBIO

# length means the length of tube=length*a
def build_mos2_nanotube(nx, ny, vacuum):
    # Approximate tube radius
    a = 3.1828  # lattice constant of MoS2 in-plane (approx.)
    mo_circumference = np.sqrt(3)* a * ny
    mo_radius = mo_circumference / (2 * np.pi)
    theta_o   = np.arctan(0.91880/(mo_radius + 1.56339))
    so_radius = 0.91880/np.sin(theta_o)
    theta_i   = np.arctan(0.91880/(mo_radius - 1.56339))
    si_radius = 0.91880/np.sin(theta_i)
    # Create flat monolayer 1T-MoS2 unit
    slab = newmx2('Mo2S4', kind='1T', a=a, size=(nx, ny, 1), vacuum=0)
    positions = slab.get_positions()
    positions = positions - positions[0][:]
    slab.set_positions(positions)
    # Roll into a cylinder
    symbols = slab.get_chemical_symbols()
    for i, pos in enumerate(positions):
        x, y, z = pos
        sym = symbols[i]

        if sym == "Mo":
           r = mo_radius
           theta = y/r
        elif sym == "S":
            if z > 0:
                r = so_radius-0.04
            elif z < 0:
                r = si_radius-0.04
            theta = y/mo_radius
        # if sym == "Mo":
        #     print(theta, np.sin(theta))
        new_y =  r * np.sin(theta)
        new_z =  r * np.cos(theta)
        positions[i] = [x, new_y, new_z]
    slab.set_positions(positions)
    return slab

def newmx2(formula='Mo2S4', kind='2H', a=3.1828, size=(1, 1, 1), vacuum=None):
#   kind : {'2H', '1T'} '2H': mirror-plane symmetry; '1T': inversion symmetry

    if kind == '2H':
        basis = [(0, 0, 0),
                 (2 / 3, 1 / 3, 0.0),
                 (2 / 3, 1 / 3, -0.0)]
    elif kind == '1T':
        basis = [(0.0, 0.0,      0.0),
                 (0.5, 0.5,      0.0),
                 (0.5, 0.5/3.0,  0.4912),
                 (0.5, 0.5/3.0, -0.4912),
                 (0.0, 2.0/3.0,  0.4912),
                 (0.0, 2.0/3.0, -0.4912)]
    else:
        raise ValueError('Structure not recognized:', kind)
#    cell = [[a, 0, 0], [-a / 2, a * 3**0.5 / 2, 0], [0, 0, a]]
    cell = [a, a * np.sqrt(3), a, 90, 90, 90]
    atoms = Atoms(formula, cell=cell, pbc=(1, 1, 0))
    atoms.set_scaled_positions(basis)
    if vacuum is not None:
        atoms.center(vacuum, axis=2)
    atoms = atoms.repeat(size)
    return atoms

def rename_resids(pdb_in:str, pdb_out:str=None):
    """
    For every 6 ATOMs, assign the same resid; resid increases by 1 each block.
    Example: atoms 1–6 resid=1, atoms 7–12 resid=2, etc.
    """
    if pdb_out is None:
        root, ext = os.path.splitext(pdb_in)
        pdb_out = root + "_residfix" + ext

    atom_idx = 0
    with open(pdb_in, "r") as fin, open(pdb_out, "w") as fout:
        for line in fin:
            if line.startswith(("ATOM")):
                atom_idx += 1
                resid = 1 + (atom_idx - 1) // 6
                parts = line.rstrip("\n").split()
                parts[4] = str(resid)
                if (atom_idx % 6) == 1:
                    parts[2] = "Mo1"
                elif (atom_idx % 6) == 2:
                    parts[2] = "Mo2"
                elif (atom_idx % 6) == 3:
                    parts[2] = "S1"
                elif (atom_idx % 6) == 4:
                    parts[2] = "S2"
                elif (atom_idx % 6) == 5:
                    parts[2] = "S3"
                elif (atom_idx % 6) == 0:
                    parts[2] = "S4"
                else:
                    print('Name_error')
# Convert to standard PDB format:
# www.biostat.jhsph.edu/~iruczins/teaching/260.655/links/pdbformat.pdf
                parts[1] = parts[1].rjust(6,' ')
                parts[2] = parts[2].rjust(4,' ')
                parts[3] = parts[3].rjust(3,' ')
                parts[4] = parts[4].rjust(5,' ')
                parts[5] = parts[5].rjust(11,' ')
                parts[6] = parts[6].rjust(7,' ')
                parts[7] = parts[7].rjust(7,' ')
                parts[8] = parts[8].rjust(5,' ')
                parts[9] = parts[9].rjust(5,' ')
                parts[10] = parts[10].rjust(11,' ')

                line = " ".join(parts) + "\n"
#                if atom_idx == 9:
#                   print(line)

            fout.write(line)
    return pdb_out

# Build the nanotube
# ny=10  for r=8.77380A in radius and 10*sqrt(3)*3.18  in circumference
# ny=20  for r=17.5477A in radius and 20*sqrt(3)*3.18  in circumference
# ny=50  for r=43.8692A in radius and 50*sqrt(3)*3.18  in circumference
# ny=100 for r=87.7385A in radius and 100*sqrt(3)*3.18 in circumference
iny = 100
mos2_tube = build_mos2_nanotube(nx=7, ny=iny, vacuum=10)

# Save as PDB data file
write(f"tmp.pdb", mos2_tube, format="proteindatabank")
# fix several issues in the pdb file
# 1. increase resid of each unit cell (each unit cell has two MoS2)
# 2. rename atom name in each unit cell with Mo1 Mo2 S1 S2 S3 S4
rename_resids(f"tmp.pdb", pdb_out=f"mos2_N{iny}.pdb")

# Save as LAMMPS data file and XYZ for visualization
write(f"mos2_N{iny}.xyz", mos2_tube)

print("MoS₂ nanotube structure generated and saved")
