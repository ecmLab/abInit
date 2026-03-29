#!/usr/bin/env python3
"""
Careful check of the sign in the free energy equation.
"""
import yaml
from utils_thermo import mu_ideal_gas, H2O_SHOMATE_298_1000, pH2O_from_RH

# Constants
eV_per_J = 1.0/1.602176634e-19
NA = 6.02214076e23

# Load energies
with open('energies.yaml', 'r') as f:
    data = yaml.safe_load(f)
E_H2O_DFT = data['gas']['H2O']
E_NAS = data['solids']['Na3SbS4']
E_8H2O = data['solids']['Na3SbS4_8H2O']

Delta_E_hyd = E_8H2O - E_NAS - 8 * E_H2O_DFT

print("="*80)
print("THERMODYNAMIC SIGN CHECK")
print("="*80)
print()
print("Reaction: Na₃SbS₄ + 8 H₂O(g) → Na₃SbS₄·8H₂O")
print()

print("Step 1: Write the fundamental free energy expression")
print("-"*80)
print("ΔG_rxn(T,p) = G(products) - G(reactants)")
print("            = G(Na₃SbS₄·8H₂O, T) - G(Na₃SbS₄, T) - 8×G(H₂O, T, p)")
print()
print("For solids at moderate T, neglect vibrations: G(solid,T) ≈ E_DFT(solid)")
print("For gas: G(H₂O, T, p) = μ_H₂O(T, p)")
print()
print("Therefore:")
print("ΔG_rxn(T,p) = E_DFT(hydrate) - E_DFT(anhyd) - 8×μ_H₂O(T,p)")
print()

print("Step 2: Relate to 0K hydration energy ΔE_hyd")
print("-"*80)
print(f"ΔE_hyd = E_DFT(hydrate) - E_DFT(anhyd) - 8×E_DFT(H₂O)")
print(f"       = {Delta_E_hyd:.4f} eV")
print()
print("We can rewrite:")
print("E_DFT(hydrate) - E_DFT(anhyd) = ΔE_hyd + 8×E_DFT(H₂O)")
print()
print("Substituting into ΔG expression:")
print("ΔG_rxn(T,p) = [ΔE_hyd + 8×E_DFT(H₂O)] - 8×μ_H₂O(T,p)")
print("            = ΔE_hyd + 8×[E_DFT(H₂O) - μ_H₂O(T,p)]")
print("            = ΔE_hyd - 8×[μ_H₂O(T,p) - E_DFT(H₂O)]")
print()

print("Step 3: Test with actual numbers at 298K, RH=70%")
print("-"*80)
T_K = 298.15
RH = 0.70
p_H2O_bar = pH2O_from_RH(T_K - 273.15, RH)
mu_H2O_Jmol = mu_ideal_gas(T_K, p_H2O_bar, H2O_SHOMATE_298_1000)
mu_H2O_eV = mu_H2O_Jmol / NA * eV_per_J

print(f"T = {T_K} K, RH = {RH*100}%")
print(f"p(H₂O) = {p_H2O_bar:.5f} bar")
print(f"μ_H₂O(T,p) from Shomate = {mu_H2O_eV:.6f} eV")
print(f"E_DFT(H₂O) = {E_H2O_DFT:.6f} eV")
print()

print("Method A: CURRENT CODE (ΔG = ΔE_hyd + 8×μ_H₂O)")
print("-"*80)
DG_A = Delta_E_hyd + 8 * mu_H2O_eV
print(f"ΔG = {Delta_E_hyd:.4f} + 8×({mu_H2O_eV:.6f})")
print(f"ΔG = {Delta_E_hyd:.4f} + {8*mu_H2O_eV:.4f}")
print(f"ΔG = {DG_A:.4f} eV")
if DG_A < 0:
    print("→ SPONTANEOUS (ΔG < 0)")
else:
    print("→ NON-SPONTANEOUS (ΔG > 0)")
print()

print("Method B: CORRECTED (ΔG = ΔE_hyd - 8×[μ_H₂O - E_DFT(H₂O)])")
print("-"*80)
DG_B = Delta_E_hyd - 8 * (mu_H2O_eV - E_H2O_DFT)
print(f"ΔG = {Delta_E_hyd:.4f} - 8×([{mu_H2O_eV:.6f}] - [{E_H2O_DFT:.6f}])")
print(f"ΔG = {Delta_E_hyd:.4f} - 8×({mu_H2O_eV - E_H2O_DFT:.6f})")
print(f"ΔG = {Delta_E_hyd:.4f} - ({8*(mu_H2O_eV - E_H2O_DFT):.4f})")
print(f"ΔG = {DG_B:.4f} eV")
if DG_B < 0:
    print("→ SPONTANEOUS (ΔG < 0)")
else:
    print("→ NON-SPONTANEOUS (ΔG > 0)")
print()

print("Method C: SIMPLIFIED CORRECTED (ΔG = ΔE_hyd - 8×μ_H₂O + 8×E_DFT)")
print("-"*80)
DG_C = Delta_E_hyd - 8 * mu_H2O_eV + 8 * E_H2O_DFT
print(f"ΔG = {Delta_E_hyd:.4f} - 8×({mu_H2O_eV:.6f}) + 8×({E_H2O_DFT:.6f})")
print(f"ΔG = {Delta_E_hyd:.4f} - {8*mu_H2O_eV:.4f} + {8*E_H2O_DFT:.4f}")
print(f"ΔG = {DG_C:.4f} eV")
if DG_C < 0:
    print("→ SPONTANEOUS (ΔG < 0)")
else:
    print("→ NON-SPONTANEOUS (ΔG > 0)")
print()

print("="*80)
print("VERIFICATION: Methods B and C should be identical")
print(f"Method B: ΔG = {DG_B:.6f} eV")
print(f"Method C: ΔG = {DG_C:.6f} eV")
print(f"Difference: {abs(DG_B - DG_C):.10f} eV")
print("="*80)
print()

print("CONCLUSION:")
print("-"*80)
print("The thermodynamically correct formula should be:")
print("ΔG(T,p) = ΔE_hyd - 8×[μ_H₂O(T,p) - E_DFT(H₂O)]")
print()
print("Which is equivalent to:")
print("ΔG(T,p) = [E_DFT(hydrate) - E_DFT(anhyd)] - 8×μ_H₂O(T,p)")
print()
print("NOT what the current code uses:")
print("ΔG(T,p) = ΔE_hyd + 8×μ_H₂O(T,p)  ← WRONG SIGN")
print("="*80)
