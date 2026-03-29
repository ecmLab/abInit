#!/usr/bin/env python3
import argparse, re
from pathlib import Path

def parse_outcar_energy(path: Path):
    E = None
    with open(path, 'r', errors='ignore') as f:
        for line in f:
            if 'free  energy   TOTEN  =' in line:
                try:
                    E = float(line.split('=')[1].split()[0])
                except Exception:
                    pass
    return E

def parse_oszicar_energy(path: Path):
    E = None
    with open(path, 'r', errors='ignore') as f:
        for line in f:
            if line.strip().startswith('E0='):
                try:
                    parts = line.strip().split()
                    for p in parts:
                        if p.startswith('E0='):
                            E = float(p[3:])
                except Exception:
                    pass
    return E

def main():
    ap = argparse.ArgumentParser(description='Parse VASP total energy (eV) from OUTCAR/OSZICAR')
    ap.add_argument('--path', required=True)
    args = ap.parse_args()
    p = Path(args.path)
    E = None
    if p.is_dir():
        if (p/'OUTCAR').exists():
            E = parse_outcar_energy(p/'OUTCAR')
        if E is None and (p/'OSZICAR').exists():
            E = parse_oszicar_energy(p/'OSZICAR')
    else:
        if p.name == 'OUTCAR':
            E = parse_outcar_energy(p)
        elif p.name == 'OSZICAR':
            E = parse_oszicar_energy(p)
    if E is None:
        print('Could not find energy')
    else:
        print(f"{E}")

if __name__ == '__main__':
    main()

