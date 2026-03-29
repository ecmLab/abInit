#!/usr/bin/env python3
"""Quick test: FeCl3 vs Li3YCl6, expected molar ratio ≈ 1.99:1"""

from mp_api.client import MPRester
from pymatgen.analysis.interface_reactions import GrandPotentialInterfacialReactivity
from pymatgen.analysis.phase_diagram import GrandPotentialPhaseDiagram, PhaseDiagram
from pymatgen.core import Composition, Element
import numpy as np

API_KEY    = "H9GD9SySZmkHxPyxI3T5ylATO9Z8XFAd"
OPEN_EL    = "Li"
V_MIN, V_MAX = 1.0, 4.2
N_VOLTAGES = 10

cathode_formula    = "FeCl3"
electrolyte_formula = "Li3YCl6"

comp1 = Composition(cathode_formula)
comp2 = Composition(electrolyte_formula)
elements = list({e.symbol for e in comp1.elements + comp2.elements} | {OPEN_EL})
print(f"Elements: {sorted(elements)}")

print("Fetching MP entries ...", flush=True)
with MPRester(API_KEY) as mpr:
    entries = mpr.get_entries_in_chemsys(elements)
print(f"{len(entries)} entries retrieved.\n")

pd_obj   = PhaseDiagram(entries)
mu_li_ref = pd_obj.get_transition_chempots(Element(OPEN_EL))[0]
print(f"μ_Li reference (most Li-rich transition): {mu_li_ref:.4f} eV\n")

n1, n2 = comp1.num_atoms, comp2.num_atoms

def print_kinks(interface, label):
    print(f"\n=== {label} ===")
    for idx, x_mix, reactivity, rxn, rxn_kJ in interface.get_kinks():
        if 0 < x_mix < 1:
            molar_ratio = (x_mix * n2) / ((1 - x_mix) * n1)
            tag = f"  ratio={molar_ratio:.3f}:1"
        else:
            tag = "  [boundary]"
        print(f"  kink {idx}: x={x_mix:.4f}  E={reactivity:.4f} eV/atom{tag}  |  {rxn}")

# Compare entries from legacy API vs new API
from pymatgen.ext.matproj import MPRester as LegacyMPRester

print("\nFetching with legacy API ...")
with LegacyMPRester(API_KEY) as mpr_old:
    entries_old = mpr_old.get_entries_in_chemsys(elements)
print(f"Legacy API: {len(entries_old)} entries")

# List stable phases for both databases
pd_old = PhaseDiagram(entries_old)
stable_old = {e.name for e in pd_old.stable_entries}
stable_new = {e.name for e in pd_obj.stable_entries}

print(f"\nStable phases in NEW API:    {sorted(stable_new)}")
print(f"Stable phases in LEGACY API: {sorted(stable_old)}")
print(f"\nNew-only phases (added since legacy): {sorted(stable_new - stable_old)}")
print(f"Legacy-only phases (removed):         {sorted(stable_old - stable_new)}")

# Run with legacy entries
mu_li_old = pd_old.get_transition_chempots(Element(OPEN_EL))[0]
for voltage in [1.0, 2.0, 3.0, 4.0, 4.2]:
    chempots = {OPEN_EL: -voltage + mu_li_old}
    gpd = GrandPotentialPhaseDiagram(entries_old, chempots)
    interface = GrandPotentialInterfacialReactivity(
        comp1, comp2, gpd,
        pd_non_grand=pd_old,
        include_no_mixing_energy=True,
        norm=True,
        use_hull_energy=False,
    )
    print(f"\n--- Legacy DB, V = {voltage:.1f} V ---")
    for idx, x_mix, reactivity, rxn, rxn_kJ in interface.get_kinks():
        if 0 < x_mix < 1:
            molar_ratio = (x_mix * n2) / ((1 - x_mix) * n1)
            tag = f"  ratio={molar_ratio:.3f}:1"
        else:
            tag = "  [boundary]"
        print(f"  kink {idx}: x={x_mix:.4f}  E={reactivity:.4f} eV/atom{tag}  |  {rxn}")
