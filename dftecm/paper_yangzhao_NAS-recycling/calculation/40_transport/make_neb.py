#!/usr/bin/env python3
import argparse, os
from pathlib import Path

def main():
    ap = argparse.ArgumentParser(description='Scaffold NEB directory tree')
    ap.add_argument('--out', required=True, help='NEB output directory')
    ap.add_argument('--images', type=int, default=5, help='Number of intermediate images')
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    nimg = args.images
    # create 00..(nimg+1) directories
    for i in range(nimg+2):
        d = out / f"{i:02d}"
        d.mkdir(exist_ok=True)
    # copy INCAR.NEB and KPOINTS template
    here = Path(__file__).parent / 'neb_template'
    os.system(f"cp '{here}/INCAR.NEB' '{out}/INCAR'")
    os.system(f"cp '{here}/KPOINTS' '{out}/KPOINTS'")
    print(f"NEB scaffold created in {out}. Place initial POSCAR in 00 and final POSCAR in {nimg+1:02d}.")

if __name__ == '__main__':
    main()

