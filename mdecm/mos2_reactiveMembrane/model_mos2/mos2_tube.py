import numpy as np
from ase.build import mx2
from ase import Atoms
from ase.io import write

# Only generate armchair (n=m) type nanotube.
# length means the length of tube=length*a
def build_mos2_nanotube(m=10, length=3, vacuum=10):
    # Approximate tube diameter and chiral vector (armchair)
    a = 3.18  # lattice constant of MoS2 in-plane (approx.)
    tube_radius = (np.sqrt(3) * a * m) / (2 * np.pi)
    tube_circumference = 2 * np.pi * tube_radius

    # Create flat monolayer MoS2 unit
    slab = mx2('MoS2', kind='2H', a=a, thickness=3.13, size=(1, 2*m, length), vacuum=0)
    slab.center(vacuum=vacuum, axis=2)

    # Roll into a cylinder
    positions = slab.get_positions()
    for i, pos in enumerate(positions):
        x, y, z = pos
        theta = 2 * np.pi * y / tube_circumference
        new_x = tube_radius * np.cos(theta)
        new_y = tube_radius * np.sin(theta)
        positions[i] = [new_x, new_y, z]

    slab.set_positions(positions)

    # Center the structure at the origin
#    center_of_mass = np.mean(positions, axis=0)
#    positions -= center_of_mass  # Shift all positions so the center of mass is at the origin
#    slab.set_positions(positions)

    coords = slab.get_positions()
    center = np.mean(coords, axis=0)
    centered_coords = coords - center
    slab.set_positions(centered_coords)

    return slab

# Build the nanotube
mos2_tube = build_mos2_nanotube(m=11, length=3, vacuum=10)

# Save as LAMMPS data file and XYZ for visualization
write('mos2_r10.data', mos2_tube, format='lammps-data', atom_style='atomic')
write('mos2_r10.xyz', mos2_tube)
write("mos2_r10.pdb", mos2_tube, format="proteindatabank")

print("MoS₂ nanotube structure generated and saved")
