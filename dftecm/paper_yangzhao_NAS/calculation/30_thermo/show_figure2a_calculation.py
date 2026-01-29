#!/usr/bin/env python3
"""
Show detailed calculation steps for Figure 2a: Dehydration Free Energy vs Temperature
"""
import yaml
from utils_thermo import mu_ideal_gas, H2O_SHOMATE_298_1000, shomate_H, shomate_S

# Constants
eV_per_J = 1.0/1.602176634e-19
NA = 6.02214076e23
R = 8.31446261815324  # J/mol/K

# Load energies
with open('energies.yaml', 'r') as f:
    data = yaml.safe_load(f)
E_H2O_DFT = data['gas']['H2O']
E_NAS = data['solids']['Na3SbS4']
E_8H2O = data['solids']['Na3SbS4_8H2O']

# Calculate 0K hydration energy
Delta_E_hyd = E_8H2O - E_NAS - 8 * E_H2O_DFT

print("="*80)
print("CALCULATION METHOD FOR FIGURE 2a: DEHYDRATION FREE ENERGY vs TEMPERATURE")
print("="*80)
print()

print("STEP 1: Define the dehydration reaction")
print("-"*80)
print("Dehydration reaction: Na₃SbS₄·8H₂O(s) → Na₃SbS₄(s) + 8 H₂O(g)")
print()
print("The free energy of this reaction is:")
print("ΔG_dehyd(T,p) = G[Na₃SbS₄(s)] + 8×G[H₂O(g,T,p)] - G[Na₃SbS₄·8H₂O(s)]")
print()
print("Rearranging:")
print("ΔG_dehyd(T,p) = -[G(hydrate) - G(anhydrous) - 8×G(H₂O,g)]")
print("              = -ΔG_hyd(T,p)")
print()

print("STEP 2: Relate to the 0K hydration energy ΔE_hyd")
print("-"*80)
print(f"From DFT: ΔE_hyd = E(hydrate) - E(anhydrous) - 8×E_DFT(H₂O)")
print(f"          ΔE_hyd = {Delta_E_hyd:.4f} eV")
print()
print("For hydration reaction: Na₃SbS₄ + 8 H₂O(g) → Na₃SbS₄·8H₂O")
print("ΔG_hyd(T,p) = E(hydrate) - E(anhydrous) - 8×μ_H₂O(T,p)")
print()
print("Since for solids at moderate T: G(solid,T) ≈ E_DFT(solid)")
print("And μ_H₂O(T,p) = E_DFT(H₂O) + ΔG_Shomate(T,p)")
print()
print("We get:")
print("ΔG_hyd(T,p) = [E(hydrate) - E(anhydrous) - 8×E_DFT(H₂O)] - 8×ΔG_Shomate(T,p)")
print("            = ΔE_hyd - 8×ΔG_Shomate(T,p)")
print()
print("Therefore, dehydration:")
print("ΔG_dehyd(T,p) = -ΔG_hyd(T,p)")
print("              = -ΔE_hyd + 8×ΔG_Shomate(T,p)")
print("              = -ΔE_hyd + 8×μ_H₂O(T,p)")
print()
print("where μ_H₂O(T,p) is the chemical potential contribution from Shomate.")
print()

print("STEP 3: Calculate μ_H₂O(T,p) using NIST Shomate equations")
print("-"*80)
print("For H₂O(g), the chemical potential is:")
print("μ_H₂O(T,p) = G(T,p) = H(T) - T×S(T) + RT ln(p/p°)")
print()
print("Using NIST Shomate equations (298-1000 K):")
print()
print("Enthalpy: H(T) - H°(298K) = A·t + (B·t²)/2 + (C·t³)/3 + (D·t⁴)/4 - E/t + F - H")
print("          where t = T/1000, result in kJ/mol")
print()
print("Entropy:  S(T) = A·ln(t) + B·t + (C·t²)/2 + (D·t³)/3 - E/(2·t²) + G")
print("          result in J/(mol·K)")
print()
print("Shomate parameters for H₂O(g) from NIST:")
for key, val in H2O_SHOMATE_298_1000.items():
    print(f"  {key} = {val}")
print()
print("Gibbs free energy at standard pressure:")
print("G°(T) = H(T) - T×S(T)")
print()
print("Then correct for pressure:")
print("μ(T,p) = G°(T) + RT ln(p/p°)")
print("where p° = 1 bar (standard pressure)")
print()

print("STEP 4: Calculate dehydration free energy")
print("-"*80)
print("Final formula used in Figure 2a:")
print()
print("ΔG_dehyd(T,p) = -ΔE_hyd + 8 × μ_H₂O(T,p)")
print()
print(f"where ΔE_hyd = {Delta_E_hyd:.4f} eV (fixed, from DFT)")
print("      μ_H₂O(T,p) = varies with T and p (from Shomate)")
print()

# Example calculation at one point
T_example_C = 80  # °C (experimental dehydration temperature)
T_example_K = T_example_C + 273.15
p_examples = [1e-6, 1e-4, 1e-2, 0.03]  # bar

print(f"EXAMPLE CALCULATION at T = {T_example_C}°C ({T_example_K} K):")
print("-"*80)

# Show Shomate calculation details for one case
p_detailed = 1e-4
print(f"\nDetailed calculation for p = {p_detailed} bar:")
print()
t = T_example_K / 1000.0
print(f"t = T/1000 = {t:.5f}")
print()

# Calculate H(T)
H_kJmol = shomate_H(T_example_K, H2O_SHOMATE_298_1000)
print(f"H(T) - H°(298K) = {H_kJmol:.4f} kJ/mol")
print(f"                = {H_kJmol*1000:.2f} J/mol")

# Calculate S(T)
S_JmolK = shomate_S(T_example_K, H2O_SHOMATE_298_1000)
print(f"S(T) = {S_JmolK:.4f} J/(mol·K)")

# Calculate G°(T)
G_kJmol = H_kJmol - (T_example_K * S_JmolK)/1000.0
print(f"\nG°(T) = H(T) - T×S(T)")
print(f"      = {H_kJmol:.4f} - {T_example_K}×{S_JmolK:.4f}/1000")
print(f"      = {G_kJmol:.4f} kJ/mol")
print(f"      = {G_kJmol*1000:.2f} J/mol")

# Add pressure correction
import math
RT_ln_p = R * T_example_K * math.log(p_detailed / 1.0)
print(f"\nPressure correction: RT ln(p/p°)")
print(f"                   = {R:.4f} × {T_example_K} × ln({p_detailed}/1.0)")
print(f"                   = {RT_ln_p:.2f} J/mol")

mu_H2O_Jmol = G_kJmol*1000.0 + RT_ln_p
print(f"\nμ_H₂O(T,p) = G°(T) + RT ln(p/p°)")
print(f"           = {G_kJmol*1000:.2f} + {RT_ln_p:.2f}")
print(f"           = {mu_H2O_Jmol:.2f} J/mol")

mu_H2O_eV = mu_H2O_Jmol / NA * eV_per_J
print(f"           = {mu_H2O_eV:.6f} eV per molecule")

# Calculate ΔG_dehyd
DG_dehyd = -Delta_E_hyd + 8 * mu_H2O_eV
print(f"\nΔG_dehyd = -ΔE_hyd + 8×μ_H₂O")
print(f"         = -({Delta_E_hyd:.4f}) + 8×({mu_H2O_eV:.6f})")
print(f"         = {-Delta_E_hyd:.4f} + {8*mu_H2O_eV:.4f}")
print(f"         = {DG_dehyd:.4f} eV")
print()

print("\nResults at T = 80°C for all pressures in Figure 2a:")
print("-"*80)
print(f"{'Pressure (bar)':<20} {'μ_H₂O (eV)':<20} {'ΔG_dehyd (eV)':<20}")
print("-"*80)

for p in p_examples:
    mu_H2O_Jmol = mu_ideal_gas(T_example_K, p, H2O_SHOMATE_298_1000)
    mu_H2O_eV = mu_H2O_Jmol / NA * eV_per_J
    DG_dehyd = -Delta_E_hyd + 8 * mu_H2O_eV
    print(f"{p:<20.2e} {mu_H2O_eV:<20.6f} {DG_dehyd:<20.4f}")

print()
print("="*80)
print("SUMMARY OF FIGURE 2a CALCULATION METHOD")
print("="*80)
print()
print("For each (T, p) point:")
print("  1. Calculate μ_H₂O(T,p) using NIST Shomate equations for H₂O(g)")
print("  2. Convert to eV per molecule: μ_H₂O_eV = μ_H₂O_Jmol / (NA × 1.602e-19)")
print("  3. Calculate: ΔG_dehyd = -ΔE_hyd + 8×μ_H₂O_eV")
print()
print("Temperature range in Figure 2a: 20-100°C")
print("Pressure curves: 1e-6, 1e-4, 1e-2, 0.03 bar")
print()
print("Key insight:")
print("  - At ambient pressure (0.03 bar): ΔG_dehyd > 0 (endothermic)")
print("  - Under vacuum (1e-6 bar): ΔG_dehyd less positive")
print("  - Lower pressure makes μ_H₂O more negative → easier dehydration")
print("="*80)
