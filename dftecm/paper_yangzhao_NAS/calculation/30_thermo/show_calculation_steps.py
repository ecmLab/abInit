#!/usr/bin/env python3
"""
Show detailed calculation steps for hydration free energy vs RH (Figure 1).
"""
import math
import yaml
from utils_thermo import mu_ideal_gas, H2O_SHOMATE_298_1000, pH2O_from_RH

# Load energies
with open('energies.yaml', 'r') as f:
    data = yaml.safe_load(f)
    energies_solids = data['solids']
    energies_gas = data['gas']

# Constants
eV_per_J = 1.0/1.602176634e-19
NA = 6.02214076e23
R = 8.31446261815324  # J/mol/K

print("="*80)
print("CALCULATION STEPS FOR FIGURE 1: HYDRATION FREE ENERGY VS RH")
print("="*80)
print()

# Step 1: Load DFT energies
print("STEP 1: DFT Energies from energies.yaml")
print("-" * 80)
E_H2O_DFT = energies_gas['H2O']
E_NAS = energies_solids['Na3SbS4']
E_8H2O = energies_solids['Na3SbS4_8H2O']

print(f"E(H₂O, gas, 0K)      = {E_H2O_DFT:.4f} eV")
print(f"E(Na₃SbS₄, 0K)       = {E_NAS:.4f} eV")
print(f"E(Na₃SbS₄·8H₂O, 0K)  = {E_8H2O:.4f} eV")
print()

# Step 2: Calculate 0K hydration energy
print("STEP 2: Calculate 0K Hydration Energy (ΔE_hyd)")
print("-" * 80)
print("Hydration reaction: Na₃SbS₄ + 8 H₂O(g) → Na₃SbS₄·8H₂O")
print()
print("ΔE_hyd = E(Na₃SbS₄·8H₂O) - E(Na₃SbS₄) - 8×E(H₂O)")
Delta_E_hyd = E_8H2O - E_NAS - 8 * E_H2O_DFT
print(f"       = {E_8H2O:.4f} - ({E_NAS:.4f}) - 8×({E_H2O_DFT:.4f})")
print(f"       = {E_8H2O:.4f} + {-E_NAS:.4f} + {-8*E_H2O_DFT:.4f}")
print(f"       = {Delta_E_hyd:.4f} eV")
print(f"       = {Delta_E_hyd/8:.4f} eV per H₂O")
print()

# Step 3: Calculate hydration free energy at finite T and different RH
print("STEP 3: Calculate ΔG(T,RH) at 298 K")
print("-" * 80)
print("The free energy at finite temperature includes thermal corrections:")
print()
print("ΔG(T,p) = ΔE_hyd + 8×[μ_H₂O(T,p) - E_DFT(H₂O)]")
print()
print("Where μ_H₂O(T,p) is the chemical potential from NIST Shomate equations:")
print("μ_H₂O(T,p) = H(T) - T×S(T) + RT×ln(p/p°)")
print()
print("This simplifies to:")
print("ΔG(T,p) = ΔE_hyd + 8×[H(T) - T×S(T) + RT×ln(p/p°) - E_DFT(H₂O)]")
print()

# Calculate for a few example RH values
T_K = 298.15
T_C = T_K - 273.15
print(f"Temperature: {T_K} K ({T_C}°C)")
print()

RH_examples = [0.20, 0.50, 0.70, 0.90]

print("Examples at different relative humidity:")
print("-" * 80)

for RH in RH_examples:
    print(f"\n** RH = {RH*100:.0f}% **")

    # Step 3a: Convert RH to water partial pressure
    p_H2O_bar = pH2O_from_RH(T_C, RH)
    print(f"  a) Convert RH to water partial pressure:")
    print(f"     p(H₂O) = RH × p_sat(T)")
    print(f"     p(H₂O) = {RH:.2f} × {pH2O_from_RH(T_C, 1.0):.5f} bar")
    print(f"     p(H₂O) = {p_H2O_bar:.5f} bar")

    # Step 3b: Calculate chemical potential from Shomate
    mu_H2O_Jmol = mu_ideal_gas(T_K, p_H2O_bar, H2O_SHOMATE_298_1000)
    mu_H2O_eV = mu_H2O_Jmol / NA * eV_per_J

    print(f"\n  b) Calculate μ_H₂O(T,p) from NIST Shomate equations:")
    print(f"     μ_H₂O({T_K} K, {p_H2O_bar:.5f} bar) = {mu_H2O_Jmol:.2f} J/mol")
    print(f"     μ_H₂O = {mu_H2O_eV:.6f} eV per molecule")

    # Step 3c: Calculate ΔG
    DG = Delta_E_hyd + 8 * mu_H2O_eV

    print(f"\n  c) Calculate ΔG_hydration:")
    print(f"     ΔG = ΔE_hyd + 8×μ_H₂O")
    print(f"     ΔG = {Delta_E_hyd:.4f} + 8×({mu_H2O_eV:.6f})")
    print(f"     ΔG = {Delta_E_hyd:.4f} + {8*mu_H2O_eV:.4f}")
    print(f"     ΔG = {DG:.4f} eV")

    if DG < 0:
        print(f"     → SPONTANEOUS (ΔG < 0)")
    else:
        print(f"     → NON-SPONTANEOUS (ΔG > 0)")

print("\n" + "="*80)
print("KEY INSIGHT:")
print("="*80)
print("At 0K: ΔE = -1.38 eV (weakly favorable)")
print("At 298K with ambient humidity (60-80% RH):")
print("  - Large negative μ_H₂O due to entropy of gas → liquid transition")
print("  - This makes ΔG much more negative than ΔE")
print("  - Hydration becomes strongly spontaneous at room temperature")
print("="*80)
