import numpy as np
import os
from ase.build import mx2
from ase import Atoms
from ase.io import write
from Bio.PDB import PDBParser, PDBIO

# length means the length of tube=length*a
def build_mos2_nanotube(nx, ny, vacuum):

    a  = 3.1828  # nearest distance of Mo-Mo
    a1 = 2.4127  # bond length of Mo-S
    theta_inc = np.pi/ny
    # Approximate Mo tube radius
    radius_mo = np.sqrt(3)/2 * a / theta_inc
    rsin      = radius_mo * np.sin(theta_inc)
    rcos      = radius_mo * np.cos(theta_inc)
    # Compute radius and initial theta of innter layer S
    zS        = np.sqrt(a1*a1 - np.square(rsin/2 + a*a/8/rsin))
    theta_i0  = np.arctan((rsin/2 - a*a/8/rsin) / (rcos - zS))
    radius_si = (rcos - zS) / np.cos(theta_i0)
    # Compute radius and initial theto of outer lyaer S
    theta_o0  = np.arctan((rsin/2 - a*a/8/rsin) / (rcos + zS))
    radius_so = (rcos + zS) / np.cos(theta_o0)

    # Create flat monolayer 1T-MoS2 unit
    slab = newmx2('Mo2S4', kind='1T', a=a, size=(nx, ny, 1), vacuum=0)
    positions = slab.get_positions()
    positions = positions - positions[0][:]
    slab.set_positions(positions)

#   Not sure why but need to do some minor correction to theta and raduius of S atoms
    if ny == 10:
        si_correction = 0.2912
        so_correction = 0.2792
    elif ny == 20:
        si_correction = 0.1456
        so_correction = 0.1422
    elif ny == 50:
        si_correction = 0.0580
        so_correction = 0.0575
    elif ny == 100:
        si_correction = 0.0289
        so_correction = 0.028459

    mo1 = np.array([0.0, 0.0, radius_mo])
    mo2 = np.array([positions[1][0], radius_mo * np.sin(theta_inc),radius_mo * np.cos(theta_inc)])
    Ntmp = 10000
#  correction of innter layer of S
    radius_si = radius_si + si_correction    # radius correction
    for itmp in range(0, Ntmp):
        dtheta  = -0.01 + itmp * 0.02/Ntmp
        if ny == 10:
            dtheta  = dtheta * 4
        theta_i  = theta_i0   + dtheta
        s4       = np.array([positions[3][0], radius_si * np.sin(theta_i),radius_si * np.cos(theta_i)])
        dis4_1   = np.linalg.norm(s4 - mo1)
        dis4_2   = np.linalg.norm(s4 - mo2)
        if np.abs(dis4_1 - dis4_2) < 0.0001:
            break
#  correct outer layer of S
    radius_so = radius_so + so_correction    # radius correction
    for itmp in range(0, Ntmp):
        dtheta  = -0.01 + itmp * 0.02/Ntmp
        if ny == 10:
            dtheta  = dtheta * 4
        theta_o  = theta_o0   + dtheta
        s3       = np.array([positions[2][0], radius_so * np.sin(theta_o),radius_so * np.cos(theta_o)])
        dis3_1   = np.linalg.norm(s3 - mo1)
        dis3_2   = np.linalg.norm(s3 - mo2)
        if np.abs(dis3_1 - dis3_2) < 0.0001:
            break

    print("Radius of Innter S layer: ", radius_si)
    print("Radius of Mo layer: ", radius_mo)
    print("Radius of Outer  S layer: ", radius_so)
    print(theta_i-theta_i0, theta_o-theta_o0, a1, dis4_1, dis4_2, dis3_1, dis3_2)

    # Roll into a cylinder
    symbols = slab.get_chemical_symbols()
    for i, pos in enumerate(positions):
        x, y, z = pos
        sym = symbols[i]
        if sym == "Mo":
           r = radius_mo
           theta = y/r
        elif sym == "S":
            if z > 0:
                r = radius_so
                theta = theta_o
            elif z < 0:
                r = radius_si
                theta = theta_i
            theta = theta + (y/(np.sqrt(3)*a) - 1/6.0) * 2*theta_inc

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
