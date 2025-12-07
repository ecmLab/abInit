#!/usr/bin/env python3
import argparse, yaml, math
from utils_thermo import mu_ideal_gas, H2O_SHOMATE_298_1000, H2S_SHOMATE_298_1000, pH2O_from_RH

def load_yaml(p):
    with open(p, 'r') as f:
        return yaml.safe_load(f)

def reaction_energy(reac, energies, T=298.15, RH=None, p_H2O_bar=None, p_H2S_bar=None):
    # Sum solids from energies (eV) and gases from chemical potentials (convert J/mol -> eV per molecule)
    eV_per_J = 1.0/1.602176634e-19
    NA = 6.02214076e23
    J_per_eV_per_molecule = 1.602176634e-19
    total_eV = 0.0
    # handle gas partial pressures
    if p_H2O_bar is None and RH is not None:
        # assume 1 bar total
        # RH here as fraction (0-1), convert at given T (C)
        p_H2O_bar = pH2O_from_RH(T-273.15, RH)
    if p_H2S_bar is None:
        p_H2S_bar = 1e-6  # default 1 ppm
    mu_H2O_Jmol = mu_ideal_gas(T, max(p_H2O_bar or 1e-9, 1e-12), H2O_SHOMATE_298_1000)
    mu_H2S_Jmol = mu_ideal_gas(T, p_H2S_bar, H2S_SHOMATE_298_1000)
    mu_H2O_eV_per = mu_H2O_Jmol/NA * eV_per_J
    mu_H2S_eV_per = mu_H2S_Jmol/NA * eV_per_J

    for sp, coeff in reac['stoichiometry'].items():
        if sp == 'H2O_g':
            total_eV += coeff * mu_H2O_eV_per
        elif sp == 'H2S_g':
            total_eV += coeff * mu_H2S_eV_per
        else:
            if sp not in energies['solids']:
                raise KeyError(f"Missing energy for solid: {sp}")
            total_eV += coeff * energies['solids'][sp]
    return total_eV

def main():
    ap = argparse.ArgumentParser(description='Compute reaction energies ΔE/ΔG')
    ap.add_argument('--energies', required=True, help='YAML with solids energies (eV/f.u.)')
    ap.add_argument('--reactions', required=True, help='reactions.yaml')
    ap.add_argument('--T', type=float, default=298.15, help='Temperature (K)')
    ap.add_argument('--RH', type=float, default=None, help='Relative humidity (0-1)')
    ap.add_argument('--pH2O', type=float, default=None, help='p_H2O in bar (overrides RH)')
    ap.add_argument('--pH2S', type=float, default=1e-6, help='p_H2S in bar')
    args = ap.parse_args()

    energies = load_yaml(args.energies)
    rxns = load_yaml(args.reactions)['reactions']
    for r in rxns:
        dE = reaction_energy(r, energies, T=args.T, RH=args.RH, p_H2O_bar=args.pH2O, p_H2S_bar=args.pH2S)
        print(f"{r['name']}: Δ(energy) = {dE:.6f} eV (per reaction as written)")

if __name__ == '__main__':
    main()

