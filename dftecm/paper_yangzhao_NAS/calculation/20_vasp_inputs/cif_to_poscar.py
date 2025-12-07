#!/usr/bin/env python3
import math
import re
import sys
from pathlib import Path

def _clean_num(s: str) -> float:
    s = s.strip().strip("'")
    s = re.sub(r"\([^)]*\)", "", s)  # remove parentheses like 0.123(4)
    return float(s)

def parse_cif(path: Path):
    a = b = c = alpha = beta = gamma = None
    atoms = []  # list of (label, elem, x,y,z, occ)
    lines = path.read_text().splitlines()
    i = 0
    # simple key-value first
    for ln in lines:
        if ln.strip().startswith('_cell_length_a'):
            a = _clean_num(ln.split()[-1])
        elif ln.strip().startswith('_cell_length_b'):
            b = _clean_num(ln.split()[-1])
        elif ln.strip().startswith('_cell_length_c'):
            c = _clean_num(ln.split()[-1])
        elif ln.strip().startswith('_cell_angle_alpha'):
            alpha = _clean_num(ln.split()[-1])
        elif ln.strip().startswith('_cell_angle_beta'):
            beta = _clean_num(ln.split()[-1])
        elif ln.strip().startswith('_cell_angle_gamma'):
            gamma = _clean_num(ln.split()[-1])
    # parse atom site loop
    # find loop_ that contains _atom_site_fract_x/y/z
    while i < len(lines):
        if lines[i].strip().startswith('loop_'):
            headers = []
            j = i + 1
            while j < len(lines) and lines[j].strip().startswith('_'):
                headers.append(lines[j].strip())
                j += 1
            if any(h.endswith('_fract_x') for h in headers):
                # determine column indices
                def idx(name):
                    for k, h in enumerate(headers):
                        if h.endswith(name):
                            return k
                    return None
                idx_label = idx('_atom_site_label')
                idx_type = idx('_atom_site_type_symbol')
                idx_x = idx('_atom_site_fract_x')
                idx_y = idx('_atom_site_fract_y')
                idx_z = idx('_atom_site_fract_z')
                idx_occ = idx('_atom_site_occupancy')
                # read data rows until blank or next loop_
                k = j
                while k < len(lines) and lines[k].strip() and not lines[k].strip().startswith('loop_') and not lines[k].strip().startswith('_'):
                    parts = re.split(r"\s+", lines[k].strip())
                    try:
                        label = parts[idx_label] if idx_label is not None else parts[idx_type]
                        elem = parts[idx_type] if idx_type is not None else re.sub(r"[^A-Za-z]", "", label)
                        # strip charges like Na+, S2-
                        elem = re.match(r"[A-Za-z]+", elem).group(0)
                        x = _clean_num(parts[idx_x])
                        y = _clean_num(parts[idx_y])
                        z = _clean_num(parts[idx_z])
                        occ = _clean_num(parts[idx_occ]) if idx_occ is not None else 1.0
                        atoms.append((label, elem, x, y, z, occ))
                    except Exception:
                        pass
                    k += 1
                break
            i = j
        i += 1
    if None in (a, b, c, alpha, beta, gamma):
        raise ValueError('Failed to parse cell parameters from CIF.')
    return (a, b, c, alpha, beta, gamma), atoms

def lattice_vectors(a, b, c, alpha_deg, beta_deg, gamma_deg):
    alpha = math.radians(alpha_deg)
    beta = math.radians(beta_deg)
    gamma = math.radians(gamma_deg)
    va = (a, 0.0, 0.0)
    vb = (b * math.cos(gamma), b * math.sin(gamma), 0.0)
    cx = c * math.cos(beta)
    cy = c * (math.cos(alpha) - math.cos(beta) * math.cos(gamma)) / (math.sin(gamma) if abs(math.sin(gamma)) > 1e-8 else 1.0)
    cz_sq = c**2 - cx**2 - cy**2
    cz = math.sqrt(max(cz_sq, 0.0))
    vc = (cx, cy, cz)
    return va, vb, vc

def write_poscar(path_out: Path, title: str, lv, atom_list):
    # group by element
    from collections import OrderedDict
    groups = OrderedDict()
    for label, elem, x, y, z, occ in atom_list:
        if occ == 0:
            continue
        groups.setdefault(elem, []).append((x, y, z))
    elems = list(groups.keys())
    counts = [len(groups[e]) for e in elems]
    with open(path_out, 'w') as f:
        f.write(f"{title}\n")
        f.write("1.0\n")
        for v in lv:
            f.write(f"  {v[0]:.10f}  {v[1]:.10f}  {v[2]:.10f}\n")
        f.write("  " + "  ".join(elems) + "\n")
        f.write("  " + "  ".join(str(c) for c in counts) + "\n")
        f.write("Direct\n")
        for e in elems:
            for (x, y, z) in groups[e]:
                f.write(f"  {x:.10f}  {y:.10f}  {z:.10f}\n")

def main():
    if len(sys.argv) < 3:
        print("Usage: cif_to_poscar.py in.cif out.POSCAR")
        sys.exit(1)
    cif = Path(sys.argv[1])
    out = Path(sys.argv[2])
    (a, b, c, alpha, beta, gamma), atoms = parse_cif(cif)
    lv = lattice_vectors(a, b, c, alpha, beta, gamma)
    write_poscar(out, cif.stem, lv, atoms)
    print(f"Wrote POSCAR to {out}")

if __name__ == '__main__':
    main()

