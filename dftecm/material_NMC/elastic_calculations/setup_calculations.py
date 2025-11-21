#!/usr/bin/env python3
"""
Setup VASP calculations for elastic properties of NMC materials
as a function of Li content.

This script:
1. Reads CIF files from crystal_strucuture folder
2. Converts to POSCAR format
3. Creates structures with different Li contents (100%, 75%, 50%, 25%)
4. Generates INCAR, KPOINTS files for relaxation and elastic calculations
"""

import os
import re
import numpy as np
from pathlib import Path
import shutil

# Base paths
BASE_DIR = Path(__file__).parent.parent
CIF_DIR = BASE_DIR / "crystal_strucuture"
CALC_DIR = BASE_DIR / "elastic_calculations"

# Li content levels to calculate
LI_CONTENTS = [1.00, 0.75, 0.50, 0.25]

def parse_cif(cif_path):
    """Parse CIF file and extract structure information."""
    with open(cif_path, 'r') as f:
        content = f.read()

    # Extract lattice parameters
    a = float(re.search(r'_cell_length_a\s+([\d.]+)', content).group(1))
    b = float(re.search(r'_cell_length_b\s+([\d.]+)', content).group(1))
    c = float(re.search(r'_cell_length_c\s+([\d.]+)', content).group(1))
    alpha = float(re.search(r'_cell_angle_alpha\s+([\d.]+)', content).group(1))
    beta = float(re.search(r'_cell_angle_beta\s+([\d.]+)', content).group(1))
    gamma = float(re.search(r'_cell_angle_gamma\s+([\d.]+)', content).group(1))

    # Convert angles to radians
    alpha_rad = np.radians(alpha)
    beta_rad = np.radians(beta)
    gamma_rad = np.radians(gamma)

    # Calculate lattice vectors
    # Standard crystallographic convention
    ax = a
    ay = 0.0
    az = 0.0

    bx = b * np.cos(gamma_rad)
    by = b * np.sin(gamma_rad)
    bz = 0.0

    cx = c * np.cos(beta_rad)
    cy = c * (np.cos(alpha_rad) - np.cos(beta_rad) * np.cos(gamma_rad)) / np.sin(gamma_rad)
    cz = np.sqrt(c**2 - cx**2 - cy**2)

    lattice = np.array([
        [ax, ay, az],
        [bx, by, bz],
        [cx, cy, cz]
    ])

    # Extract atomic positions
    atoms = []

    # Find the atom site loop
    # Try different CIF formats
    patterns = [
        # Format 1: label, occupancy, fract_x, fract_y, fract_z, ...
        r'(\w+)\s+(\d+\.?\d*)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+\w+\s+[\d.?]+\s+(\w+)',
        # Format 2: label, type, fract_x, fract_y, fract_z, occupancy
        r'(\w+)\s+(\w+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([\d.]+)',
    ]

    # Try to find atom positions after loop_ _atom_site
    atom_section = re.search(r'loop_\s*\n\s*_atom_site', content)
    if atom_section:
        # Get everything after the loop header
        after_loop = content[atom_section.start():]

        # Find all atom lines (lines that start with element symbol)
        lines = after_loop.split('\n')
        in_data = False
        for line in lines:
            line = line.strip()
            if not line or line.startswith('#'):
                if in_data:
                    break
                continue
            if line.startswith('_'):
                continue
            if line.startswith('loop_') and in_data:
                break

            # Try to parse atom line
            parts = line.split()
            if len(parts) >= 5:
                label = parts[0]
                # Extract element from label (e.g., Li1 -> Li, Ni10 -> Ni)
                element = re.match(r'([A-Za-z]+)', label).group(1)

                # Find coordinates - they should be decimal numbers
                coords = []
                for p in parts[1:]:
                    try:
                        val = float(p)
                        if -2 < val < 2:  # Likely a fractional coordinate
                            coords.append(val)
                        if len(coords) == 3:
                            break
                    except ValueError:
                        if element in ['Li', 'Ni', 'Mn', 'Co', 'O'] and len(coords) == 0:
                            # This might be the element symbol column
                            element = p if p in ['Li', 'Ni', 'Mn', 'Co', 'O'] else element

                if len(coords) == 3:
                    atoms.append({
                        'element': element,
                        'position': coords,
                        'label': label
                    })
                    in_data = True

    return {
        'lattice': lattice,
        'a': a, 'b': b, 'c': c,
        'alpha': alpha, 'beta': beta, 'gamma': gamma,
        'atoms': atoms
    }


def write_poscar(structure, filename, title="NMC structure"):
    """Write structure to POSCAR format."""
    lattice = structure['lattice']
    atoms = structure['atoms']

    # Group atoms by element in VASP order: Li, Ni, Mn, Co, O
    element_order = ['Li', 'Ni', 'Mn', 'Co', 'O']
    grouped = {el: [] for el in element_order}

    for atom in atoms:
        el = atom['element']
        if el in grouped:
            grouped[el].append(atom['position'])

    # Remove empty groups
    elements = [el for el in element_order if grouped[el]]
    counts = [len(grouped[el]) for el in elements]

    with open(filename, 'w') as f:
        f.write(f"{title}\n")
        f.write("1.0\n")

        # Lattice vectors
        for vec in lattice:
            f.write(f"  {vec[0]:16.10f}  {vec[1]:16.10f}  {vec[2]:16.10f}\n")

        # Element symbols and counts
        f.write("  " + "  ".join(elements) + "\n")
        f.write("  " + "  ".join(map(str, counts)) + "\n")

        # Direct coordinates
        f.write("Direct\n")
        for el in elements:
            for pos in grouped[el]:
                # Wrap coordinates to [0, 1)
                pos = [p - np.floor(p) for p in pos]
                f.write(f"  {pos[0]:16.10f}  {pos[1]:16.10f}  {pos[2]:16.10f}\n")

    return elements, counts


def remove_li_atoms(structure, fraction_to_keep):
    """Remove Li atoms to achieve target Li content."""
    atoms = structure['atoms'].copy()
    li_atoms = [a for a in atoms if a['element'] == 'Li']
    other_atoms = [a for a in atoms if a['element'] != 'Li']

    n_li_total = len(li_atoms)
    n_li_keep = int(round(n_li_total * fraction_to_keep))

    if n_li_keep == 0:
        kept_li = []
    elif n_li_keep >= n_li_total:
        kept_li = li_atoms
    else:
        # Remove Li atoms uniformly distributed
        # Use deterministic selection based on position
        li_atoms_sorted = sorted(li_atoms, key=lambda a: (a['position'][2], a['position'][0], a['position'][1]))
        step = n_li_total / n_li_keep
        indices = [int(i * step) for i in range(n_li_keep)]
        kept_li = [li_atoms_sorted[i] for i in indices]

    new_structure = structure.copy()
    new_structure['atoms'] = other_atoms + kept_li
    return new_structure, n_li_keep, n_li_total


def write_incar_relax(filename, system_name):
    """Write INCAR for structural relaxation."""
    content = f"""# INCAR for structural relaxation of {system_name}
# Elastic property calculation - Step 1: Relaxation

SYSTEM = {system_name}

# Electronic relaxation
ENCUT = 600        # Plane-wave cutoff (eV)
PREC = Accurate
EDIFF = 1E-6       # Electronic convergence
NELMIN = 4
NELM = 200
LREAL = Auto       # Real space projection

# Ionic relaxation
IBRION = 2         # Conjugate gradient
ISIF = 3           # Relax ions, cell shape, and volume
NSW = 200          # Max ionic steps
EDIFFG = -0.01     # Force convergence (eV/A)

# DFT+U settings (Dudarev approach)
LDAU = .TRUE.
LDAUTYPE = 2
LDAUL = -1 2 2 2 -1    # Li Ni Mn Co O (d-orbitals for TM)
LDAUU = 0.0 6.0 4.5 3.3 0.0
LDAUJ = 0.0 0.0 0.0 0.0 0.0
LMAXMIX = 4

# Spin polarization
ISPIN = 2          # Spin-polarized
MAGMOM = 96*0.0    # Will be adjusted based on system size

# Smearing
ISMEAR = 0         # Gaussian smearing
SIGMA = 0.05

# Output
LWAVE = .FALSE.
LCHARG = .FALSE.
LORBIT = 11        # DOSCAR and lm-decomposed PROCAR

# Parallelization
NCORE = 4
KPAR = 2

# vdW correction (important for layered materials)
IVDW = 12          # DFT-D3 with BJ damping
"""
    with open(filename, 'w') as f:
        f.write(content)


def write_incar_elastic(filename, system_name):
    """Write INCAR for elastic constant calculation."""
    content = f"""# INCAR for elastic constant calculation of {system_name}
# Elastic property calculation - Step 2: Elastic constants

SYSTEM = {system_name}

# Electronic relaxation
ENCUT = 600        # Plane-wave cutoff (eV)
PREC = Accurate
EDIFF = 1E-7       # Tight electronic convergence
NELMIN = 4
NELM = 200
LREAL = Auto

# Elastic constants calculation
IBRION = 6         # Finite differences for elastic tensor
ISIF = 3           # Calculate stress tensor
NFREE = 4          # Number of displacements (2 or 4)
POTIM = 0.015      # Step size for finite differences

# DFT+U settings (same as relaxation)
LDAU = .TRUE.
LDAUTYPE = 2
LDAUL = -1 2 2 2 -1
LDAUU = 0.0 6.0 4.5 3.3 0.0
LDAUJ = 0.0 0.0 0.0 0.0 0.0
LMAXMIX = 4

# Spin polarization
ISPIN = 2
MAGMOM = 96*0.0    # Will be adjusted

# Smearing
ISMEAR = 0
SIGMA = 0.05

# Output
LWAVE = .FALSE.
LCHARG = .FALSE.
LORBIT = 11

# Parallelization
NCORE = 4
KPAR = 2

# vdW correction
IVDW = 12
"""
    with open(filename, 'w') as f:
        f.write(content)


def write_kpoints(filename, lattice, target_density=40):
    """
    Write KPOINTS file based on lattice vectors.
    target_density: k-points per reciprocal Angstrom
    """
    # Calculate reciprocal lattice vector lengths
    a_len = np.linalg.norm(lattice[0])
    b_len = np.linalg.norm(lattice[1])
    c_len = np.linalg.norm(lattice[2])

    # Calculate k-point grid
    ka = max(1, int(round(target_density / a_len)))
    kb = max(1, int(round(target_density / b_len)))
    kc = max(1, int(round(target_density / c_len)))

    content = f"""Automatic mesh
0
Gamma
{ka} {kb} {kc}
0 0 0
"""
    with open(filename, 'w') as f:
        f.write(content)

    return ka, kb, kc


def write_job_script(filename, system_name, calc_type):
    """Write a sample job submission script."""
    content = f"""#!/bin/bash
#SBATCH --job-name={system_name}_{calc_type}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=24:00:00
#SBATCH --partition=normal

module load vasp/6.3.0

cd $SLURM_SUBMIT_DIR

mpirun -np $SLURM_NTASKS vasp_std > vasp.out 2>&1

echo "Job completed at $(date)"
"""
    with open(filename, 'w') as f:
        f.write(content)


def main():
    """Main function to set up all calculations."""

    # Find all CIF files
    cif_files = list(CIF_DIR.glob("*.cif"))
    print(f"Found {len(cif_files)} CIF files")

    # Create summary file
    summary_file = CALC_DIR / "calculation_summary.txt"

    with open(summary_file, 'w') as summary:
        summary.write("NMC Elastic Property Calculations Summary\n")
        summary.write("=" * 60 + "\n\n")

        for cif_file in sorted(cif_files):
            nmc_name = cif_file.stem  # e.g., "NMC-111"
            print(f"\nProcessing {nmc_name}...")
            summary.write(f"\n{nmc_name}\n")
            summary.write("-" * 40 + "\n")

            # Parse CIF file
            try:
                structure = parse_cif(cif_file)
            except Exception as e:
                print(f"  Error parsing {cif_file}: {e}")
                summary.write(f"  ERROR: Could not parse CIF file: {e}\n")
                continue

            n_li_original = len([a for a in structure['atoms'] if a['element'] == 'Li'])
            n_atoms_total = len(structure['atoms'])

            summary.write(f"  Original: {n_atoms_total} atoms, {n_li_original} Li\n")
            summary.write(f"  Lattice: a={structure['a']:.3f}, b={structure['b']:.3f}, c={structure['c']:.3f}\n")

            # Create folder for this NMC type
            nmc_dir = CALC_DIR / nmc_name
            nmc_dir.mkdir(exist_ok=True)

            # Create calculations for each Li content
            for li_frac in LI_CONTENTS:
                li_pct = int(li_frac * 100)
                li_dir = nmc_dir / f"Li_{li_pct:03d}"

                # Create subdirectories
                relax_dir = li_dir / "relax"
                elastic_dir = li_dir / "elastic"
                relax_dir.mkdir(parents=True, exist_ok=True)
                elastic_dir.mkdir(parents=True, exist_ok=True)

                # Create structure with modified Li content
                mod_structure, n_li_kept, n_li_total = remove_li_atoms(structure, li_frac)

                if n_li_total == 0 and li_frac > 0:
                    # Skip if no Li in original structure (e.g., NMC-RS)
                    print(f"  Skipping Li_{li_pct:03d} - no Li in structure")
                    continue

                system_name = f"{nmc_name}_Li{li_pct}"

                # Write POSCAR files
                poscar_title = f"{nmc_name} Li={li_frac:.2f} ({n_li_kept}/{n_li_total} Li)"
                elements, counts = write_poscar(mod_structure, relax_dir / "POSCAR", poscar_title)
                shutil.copy(relax_dir / "POSCAR", elastic_dir / "POSCAR")

                # Calculate total atoms for MAGMOM
                total_atoms = sum(counts)

                # Write INCAR files
                write_incar_relax(relax_dir / "INCAR", system_name)
                write_incar_elastic(elastic_dir / "INCAR", system_name)

                # Update MAGMOM in INCAR files based on actual atom count
                for incar_path in [relax_dir / "INCAR", elastic_dir / "INCAR"]:
                    with open(incar_path, 'r') as f:
                        content = f.read()
                    # Create proper MAGMOM string
                    magmom_parts = []
                    for el, cnt in zip(elements, counts):
                        if el == 'Ni':
                            magmom_parts.append(f"{cnt}*1.0")  # Ni: ~1 muB
                        elif el == 'Mn':
                            magmom_parts.append(f"{cnt}*4.0")  # Mn: ~4 muB (high spin)
                        elif el == 'Co':
                            magmom_parts.append(f"{cnt}*0.5")  # Co: ~0.5 muB (low spin)
                        else:
                            magmom_parts.append(f"{cnt}*0.0")
                    magmom_str = " ".join(magmom_parts)
                    content = re.sub(r'MAGMOM = .*', f'MAGMOM = {magmom_str}', content)

                    # Update LDAUL based on elements present
                    ldaul_parts = []
                    ldauu_parts = []
                    for el in elements:
                        if el == 'Li':
                            ldaul_parts.append("-1")
                            ldauu_parts.append("0.0")
                        elif el == 'Ni':
                            ldaul_parts.append("2")
                            ldauu_parts.append("6.0")
                        elif el == 'Mn':
                            ldaul_parts.append("2")
                            ldauu_parts.append("4.5")
                        elif el == 'Co':
                            ldaul_parts.append("2")
                            ldauu_parts.append("3.3")
                        elif el == 'O':
                            ldaul_parts.append("-1")
                            ldauu_parts.append("0.0")

                    content = re.sub(r'LDAUL = .*', f'LDAUL = {" ".join(ldaul_parts)}', content)
                    content = re.sub(r'LDAUU = .*', f'LDAUU = {" ".join(ldauu_parts)}', content)
                    content = re.sub(r'LDAUJ = .*', f'LDAUJ = {" ".join(["0.0"]*len(elements))}', content)

                    with open(incar_path, 'w') as f:
                        f.write(content)

                # Write KPOINTS files
                ka, kb, kc = write_kpoints(relax_dir / "KPOINTS", structure['lattice'])
                shutil.copy(relax_dir / "KPOINTS", elastic_dir / "KPOINTS")

                # Write job scripts
                write_job_script(relax_dir / "submit.sh", system_name, "relax")
                write_job_script(elastic_dir / "submit.sh", system_name, "elastic")

                print(f"  Created Li_{li_pct:03d}: {n_li_kept}/{n_li_total} Li, "
                      f"k-mesh: {ka}x{kb}x{kc}")
                summary.write(f"  Li_{li_pct:03d}: {n_li_kept}/{n_li_total} Li, "
                             f"{sum(counts)} atoms, k-mesh: {ka}x{kb}x{kc}\n")

    # Create master run script
    master_script = CALC_DIR / "run_all.sh"
    with open(master_script, 'w') as f:
        f.write("""#!/bin/bash
# Master script to submit all calculations
# First run relaxations, then elastic calculations

echo "Submitting relaxation calculations..."
for dir in */Li_*/relax; do
    if [ -d "$dir" ]; then
        echo "Submitting $dir"
        cd "$dir"
        # Uncomment the line below for your scheduler:
        # sbatch submit.sh
        cd - > /dev/null
    fi
done

echo ""
echo "After relaxations complete, copy CONTCAR to elastic/POSCAR:"
echo "for dir in */Li_*/relax; do"
echo '    elastic_dir="${dir/relax/elastic}"'
echo '    cp "$dir/CONTCAR" "$elastic_dir/POSCAR"'
echo "done"
echo ""
echo "Then submit elastic calculations..."
""")

    # Create POTCAR generation helper
    potcar_helper = CALC_DIR / "generate_potcar.sh"
    with open(potcar_helper, 'w') as f:
        f.write("""#!/bin/bash
# Helper script to generate POTCAR files
# Modify VASP_PP_PATH to your pseudopotential directory

VASP_PP_PATH="/path/to/vasp/potentials"

# Recommended PAW potentials for NMC:
# Li_sv (1s2s2p), Ni_pv (3p4s3d), Mn_pv (3p4s3d), Co (4s3d), O (2s2p)

for dir in */Li_*/relax */Li_*/elastic; do
    if [ -d "$dir" ]; then
        echo "Generating POTCAR for $dir"
        cd "$dir"

        # Read elements from POSCAR line 6
        elements=$(sed -n '6p' POSCAR)

        # Generate POTCAR
        rm -f POTCAR
        for el in $elements; do
            case $el in
                Li) cat "$VASP_PP_PATH/PAW_PBE/Li_sv/POTCAR" >> POTCAR ;;
                Ni) cat "$VASP_PP_PATH/PAW_PBE/Ni_pv/POTCAR" >> POTCAR ;;
                Mn) cat "$VASP_PP_PATH/PAW_PBE/Mn_pv/POTCAR" >> POTCAR ;;
                Co) cat "$VASP_PP_PATH/PAW_PBE/Co/POTCAR" >> POTCAR ;;
                O)  cat "$VASP_PP_PATH/PAW_PBE/O/POTCAR" >> POTCAR ;;
            esac
        done

        cd - > /dev/null
    fi
done

echo "POTCAR generation complete!"
""")

    print(f"\n{'='*60}")
    print("Setup complete!")
    print(f"Calculation directory: {CALC_DIR}")
    print(f"Summary file: {summary_file}")
    print("\nNext steps:")
    print("1. Edit generate_potcar.sh with your VASP_PP_PATH")
    print("2. Run: bash generate_potcar.sh")
    print("3. Submit relaxation jobs first")
    print("4. After relaxations, copy CONTCAR to elastic/POSCAR")
    print("5. Submit elastic calculations")


if __name__ == "__main__":
    main()
