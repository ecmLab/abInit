import numpy as np
from ase.build import mx2
from ase import Atoms
from ase.io import write

# length means the length of tube=length*a
def build_mos2_nanotube(m=20, length=3, vacuum=10):
    # Approximate tube diameter and chiral vector (armchair)
    a = 3.18  # lattice constant of MoS2 in-plane (approx.)
    tube_circumference = a * m
    tube_radius = tube_circumference / (2 * np.pi)

    # Create flat monolayer MoS2 unit
    slab = mx2('MoS2', kind='1T', a=a, thickness=3.19, size=(1, m, length), vacuum=0)
    slab.center(vacuum=vacuum, axis=2)

    # Roll into a cylinder
    anglCon = 2 * np.pi / tube_circumference
    positions = slab.get_positions()
    symbols = slab.get_chemical_symbols()  # <-- This line is missing
    for i, pos in enumerate(positions):
        x, y, z = pos
        sym = symbols[i]

        if sym == "Mo":
           r = tube_radius
           theta = 2*(-x)*anglCon
        elif sym == "S":
            if x > 1:
                dl    = 2*y
            else:
                dl    = (y - np.sqrt(3)*np.abs(x))/2

            r     = tube_radius + dl
# classify S atoms into inner/outer shell based on z or proximity to Mo
            if dl > 1 or x >1:    # Outer S layer
                theta = 2*(1.59-x)*anglCon
            elif dl < 1:                 # Inner S layer
                theta = 2*(-x)*anglCon + 1.59*anglCon

        new_x =  r * np.cos(theta)
        new_y =  r * np.sin(theta)
        positions[i] = [new_x, new_y, z]

    slab.set_positions(positions)

    coords = slab.get_positions()
    center = np.mean(coords, axis=0)
    centered_coords = coords - center
    slab.set_positions(centered_coords)

    return slab

# Build the nanotube
# m=20  for r=1.0120nm in radius corresponds to 20*3.18 in circumference
# m=40  for r=2.0240nm in radius corresponds to 40*3.18 in circumference
# m=100 for r=5.0611nm in radius corresponds to 100*3.18 in circumference
# m=200 for r=10.122nm in radius corresponds to 200*3.18 in circumference
mos2_tube = build_mos2_nanotube(m=40, length=7, vacuum=10)

# Save as LAMMPS data file and XYZ for visualization
write('mos2_r20.data', mos2_tube, format='lammps-data', atom_style='atomic')
write('mos2_r20.xyz', mos2_tube)
write("mos2_r20.pdb", mos2_tube, format="proteindatabank")

print("MoS₂ nanotube structure generated and saved")
