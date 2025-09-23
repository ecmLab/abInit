from pymatgen.core import Lattice, Structure
from pymatgen.io.vasp import Poscar

# Approximate lattice parameter for spinel (Fd-3m)
a = 8.25
lattice = Lattice.cubic(a)

# Define species (24 total)
# 8a: Li × 8
# 16d: Ti × 6, Cr × 1, Li × 1 (Cr/Ti/Li mixed site)
# 32e: O × 8
species = ['Li']*8 + ['Ti']*6 + ['Cr'] + ['Li'] + ['O']*8

# Define corresponding fractional coordinates (24 total)
coords = [
    # Li at 8a
    [0.125, 0.125, 0.125],
    [0.875, 0.875, 0.125],
    [0.875, 0.125, 0.875],
    [0.125, 0.875, 0.875],
    [0.375, 0.375, 0.375],
    [0.625, 0.625, 0.375],
    [0.625, 0.375, 0.625],
    [0.375, 0.625, 0.625],

    # Ti at 16d (6 atoms)
    [0.5, 0.5, 0.5],
    [0.0, 0.0, 0.5],
    [0.0, 0.5, 0.0],
    [0.5, 0.0, 0.0],
    [0.25, 0.75, 0.75],
    [0.75, 0.25, 0.25],

    # Cr at 16d (1 atom)
    [0.25, 0.25, 0.75],

    # Li stuffed into 16d (1 atom)
    [0.75, 0.75, 0.25],

    # O at 32e (8 atoms, approximate)
    [0.261, 0.261, 0.261],
    [0.739, 0.739, 0.261],
    [0.739, 0.261, 0.739],
    [0.261, 0.739, 0.739],
    [0.113, 0.113, 0.613],
    [0.887, 0.887, 0.613],
    [0.387, 0.113, 0.113],
    [0.613, 0.887, 0.887],
]

# Create structure and write POSCAR
structure = Structure(lattice, species, coords)
poscar = Poscar(structure)
poscar.write_file("POSCAR_Li1.25Cr0.25Ti1.5O4.vasp")
print("✅ POSCAR written to POSCAR_Li1.25Cr0.25Ti1.5O4.vasp")
