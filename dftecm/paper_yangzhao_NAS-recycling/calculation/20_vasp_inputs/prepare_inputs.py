#!/usr/bin/env python3
import argparse
import shutil
from math import ceil
from pathlib import Path

from cif_to_poscar import lattice_vectors, parse_cif, write_poscar


def _dot(u, v):
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]


def _cross(u, v):
    return (
        u[1]*v[2] - u[2]*v[1],
        u[2]*v[0] - u[0]*v[2],
        u[0]*v[1] - u[1]*v[0],
    )


def _norm(v):
    return (v[0]**2 + v[1]**2 + v[2]**2) ** 0.5


def kmesh_from_length(lv, target_ksp=0.2):
    a_vec, b_vec, c_vec = lv
    vol = _dot(a_vec, _cross(b_vec, c_vec))
    b1 = tuple(2 * 3.141592653589793 * x / vol for x in _cross(b_vec, c_vec))
    b2 = tuple(2 * 3.141592653589793 * x / vol for x in _cross(c_vec, a_vec))
    b3 = tuple(2 * 3.141592653589793 * x / vol for x in _cross(a_vec, b_vec))
    lens = [_norm(b1), _norm(b2), _norm(b3)]
    return [max(1, int(ceil(L / target_ksp))) for L in lens]


def write_kpoints(path_out: Path, mesh):
    with open(path_out, 'w') as f:
        f.write("Automatic mesh\n0\nGamma\n")
        f.write(f"{mesh[0]} {mesh[1]} {mesh[2]}\n")
        f.write("0 0 0\n")


def main():
    ap = argparse.ArgumentParser(description="Prepare VASP inputs from CIF")
    ap.add_argument('--cif', required=True)
    ap.add_argument('--out', required=True, help='Output directory base (created)')
    ap.add_argument('--ksp', type=float, default=0.2, help='Target k-point spacing (1/Ang)')
    args = ap.parse_args()

    cif = Path(args.cif)
    out = Path(args.out)
    (a, b, c, alpha, beta, gamma), atoms = parse_cif(cif)
    lv = lattice_vectors(a, b, c, alpha, beta, gamma)
    mesh = kmesh_from_length(lv, args.ksp)

    template_dir = Path(__file__).parent / 'templates'
    for mode in ['relax', 'static']:
        d = out / mode
        d.mkdir(parents=True, exist_ok=True)
        write_poscar(d / 'POSCAR', cif.stem, lv, atoms)
        write_kpoints(d / 'KPOINTS', mesh)
        shutil.copy(template_dir / f'INCAR.{mode}', d / 'INCAR')
    print(f"Prepared VASP inputs in {out}")


if __name__ == '__main__':
    main()
