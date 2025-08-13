import numpy as np
from ase.build import mx2
from ase import Atoms
from ase.io import write

# length means the length of tube=length*a
def build_mos2_nanotube(nx, ny, vacuum):
    # Approximate tube diameter
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

    # coords = slab.get_positions()
    # center = np.mean(coords, axis=0)
    # centered_coords = coords - center
    # slab.set_positions(centered_coords)

    return slab

def newmx2(formula='Mo2S4', kind='2H', a=3.1828, size=(1, 1, 1), vacuum=None):
#   kind : {'2H', '1T'} '2H': mirror-plane symmetry; '1T': inversion symmetry

    if kind == '2H':
        basis = [(0, 0, 0),
                 (2 / 3, 1 / 3, 0.0),
                 (2 / 3, 1 / 3, -0.0)]
    elif kind == '1T':
        # basis = [(2/3, 1/3,  0.0),
        #          (1/3, 2/3,  0.4912),
        #          (1/3, 2/3, -0.4912)]
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

# Build the nanotube
# ny=10  for r=8.77380A in radius and 10*sqrt(3)*3.18  in circumference
# ny=20  for r=17.5477A in radius and 20*sqrt(3)*3.18  in circumference
# ny=50  for r=43.8692A in radius and 50*sqrt(3)*3.18  in circumference
# ny=100 for r=87.7385A in radius and 100*sqrt(3)*3.18 in circumference
mos2_tube = build_mos2_nanotube(nx=7, ny=20, vacuum=10)

# Save as LAMMPS data file and XYZ for visualization
#write('mos2_r50.data', mos2_tube, format='lammps-data', atom_style='atomic')
write('mos2_r18.xyz', mos2_tube)
write("mos2_r18.pdb", mos2_tube, format="proteindatabank")

print("MoS₂ nanotube structure generated and saved")
