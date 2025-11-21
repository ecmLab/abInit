#!/usr/bin/env python3
"""
Update all INCAR files to use SCAN+rVV10 functional instead of PBE+U+D3.
"""

import os
from pathlib import Path

CALC_DIR = Path(__file__).parent

def update_incar(incar_path):
    """Update a single INCAR file to use SCAN+rVV10."""
    with open(incar_path, 'r') as f:
        lines = f.readlines()

    new_lines = []
    skip_lines = {'LDAU', 'LDAUTYPE', 'LDAUL', 'LDAUU', 'LDAUJ', 'LMAXMIX', 'IVDW'}

    for line in lines:
        # Check if this line should be skipped
        stripped = line.strip()
        if stripped.startswith('#'):
            # Update comments about DFT+U
            if 'DFT+U' in stripped or 'Dudarev' in stripped:
                new_lines.append('# SCAN+rVV10 functional (no U correction needed)\n')
                continue
            elif 'vdW correction' in stripped.lower():
                new_lines.append('# van der Waals: rVV10 (built into SCAN+rVV10)\n')
                continue

        # Skip DFT+U and D3 related lines
        skip = False
        for key in skip_lines:
            if stripped.startswith(key):
                skip = True
                break

        if skip:
            continue

        new_lines.append(line)

    # Find where to insert SCAN+rVV10 settings (after LREAL line)
    insert_idx = None
    for i, line in enumerate(new_lines):
        if 'LREAL' in line:
            insert_idx = i + 1
            break

    # SCAN+rVV10 settings to add
    scan_settings = """
# SCAN+rVV10 functional (meta-GGA with nonlocal vdW)
METAGGA = SCAN
LUSE_VDW = .TRUE.
BPARAM = 15.7      # rVV10 parameter for SCAN
LASPH = .TRUE.     # Non-spherical contributions (required for meta-GGA)
ADDGRID = .TRUE.   # Additional support grid for accuracy

"""

    if insert_idx:
        new_lines.insert(insert_idx, scan_settings)
    else:
        # Insert after electronic relaxation section
        for i, line in enumerate(new_lines):
            if 'Electronic relaxation' in line:
                # Find the end of this section
                for j in range(i+1, len(new_lines)):
                    if new_lines[j].strip() == '' or new_lines[j].startswith('#'):
                        if new_lines[j].strip() == '':
                            new_lines.insert(j+1, scan_settings)
                            break
                break

    # Write updated INCAR
    with open(incar_path, 'w') as f:
        f.writelines(new_lines)

    return True


def main():
    # Find all INCAR files
    incar_files = list(CALC_DIR.glob("*/Li_*/*/INCAR"))

    print(f"Found {len(incar_files)} INCAR files to update")

    updated = 0
    for incar_path in sorted(incar_files):
        try:
            update_incar(incar_path)
            print(f"  Updated: {incar_path.relative_to(CALC_DIR)}")
            updated += 1
        except Exception as e:
            print(f"  ERROR: {incar_path}: {e}")

    print(f"\nUpdated {updated}/{len(incar_files)} INCAR files to SCAN+rVV10")
    print("\nKey changes:")
    print("  - Removed: LDAU, LDAUTYPE, LDAUL, LDAUU, LDAUJ, LMAXMIX, IVDW")
    print("  - Added: METAGGA=SCAN, LUSE_VDW=.TRUE., BPARAM=15.7, LASPH=.TRUE.")


if __name__ == "__main__":
    main()
