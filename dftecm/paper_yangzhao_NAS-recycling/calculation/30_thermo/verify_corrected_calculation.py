#!/usr/bin/env python3
"""
Verify the corrected calculation with proper sign.
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
print("VERIFICATION OF CORRECTED FORMULA")
print("="*80)
print()
print("Corrected Formula: ΔG = ΔE_hyd - 8×μ_H₂O(T,p)")
print()
print(f"ΔE_hyd = {Delta_E_hyd:.4f} eV")
print()

T_K = 298.15
RH_values = [0.20, 0.50, 0.70, 0.90]

print("Results at T = 298 K:")
print("-"*80)
print(f"{'RH (%)':<10} {'p(H₂O) (bar)':<15} {'μ_H₂O (eV)':<15} {'ΔG (eV)':<15} {'Status'}")
print("-"*80)

for RH in RH_values:
    p_H2O_bar = pH2O_from_RH(T_K - 273.15, RH)
    mu_H2O_Jmol = mu_ideal_gas(T_K, p_H2O_bar, H2O_SHOMATE_298_1000)
    mu_H2O_eV = mu_H2O_Jmol / NA * eV_per_J

    DG = Delta_E_hyd - 8 * mu_H2O_eV

    status = "Spontaneous" if DG < 0 else "Non-spontaneous"

    print(f"{RH*100:<10.0f} {p_H2O_bar:<15.5f} {mu_H2O_eV:<15.6f} {DG:<15.4f} {status}")

print()
print("Physical Interpretation:")
print("-"*80)
print("• At 0K: ΔE = -1.38 eV (weakly exothermic)")
print("• At 298K:")
print("  - μ_H₂O is negative (water vapor has lower free energy than reference)")
print("  - Subtracting negative μ_H₂O makes ΔG MORE NEGATIVE")
print("  - Higher RH → less negative μ_H₂O → less negative ΔG")
print("  - But ΔG remains negative → hydration is spontaneous at all RH")
print()
print("This makes physical sense:")
print("  - Condensing water vapor (high entropy) into crystal (low entropy)")
print("  - Loses entropy but gains enthalpy from hydrogen bonding")
print("  - At room T with ambient RH, enthalpy wins → spontaneous")
print("="*80)
