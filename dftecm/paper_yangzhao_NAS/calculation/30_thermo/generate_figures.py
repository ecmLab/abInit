#!/usr/bin/env python3
"""
Generate publication-quality figures for Na3SbS4 recycling thermodynamics.
Produces phase diagrams and energy plots for manuscript.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Rectangle
import yaml
from utils_thermo import mu_ideal_gas, H2O_SHOMATE_298_1000, pH2O_from_RH

# Set publication-quality matplotlib parameters
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['font.size'] = 10
mpl.rcParams['axes.linewidth'] = 1.2
mpl.rcParams['xtick.major.width'] = 1.2
mpl.rcParams['ytick.major.width'] = 1.2
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['mathtext.default'] = 'regular'

# Load energies
with open('energies.yaml', 'r') as f:
    energies = yaml.safe_load(f)['solids']

# Constants
eV_per_J = 1.0/1.602176634e-19
NA = 6.02214076e23

# DFT energies (eV per formula unit)
E_NAS = energies['Na3SbS4']
E_8H2O = energies['Na3SbS4_8H2O']
E_9H2O = energies['Na3SbS4_9H2O']

# Solid-state energy differences
Delta_E_8 = E_8H2O - E_NAS
Delta_E_9 = E_9H2O - E_NAS

def hydration_free_energy(T_K, RH, n_H2O, Delta_E_solids):
    """Calculate ΔG for hydration reaction."""
    p_H2O_bar = pH2O_from_RH(T_K - 273.15, RH)
    mu_H2O_Jmol = mu_ideal_gas(T_K, p_H2O_bar, H2O_SHOMATE_298_1000)
    mu_H2O_eV = mu_H2O_Jmol / NA * eV_per_J
    DG = Delta_E_solids + n_H2O * mu_H2O_eV
    return DG

def dehydration_free_energy_vacuum(T_K, p_vac_bar, n_H2O, Delta_E_solids):
    """Calculate ΔG for dehydration under vacuum."""
    mu_H2O_Jmol = mu_ideal_gas(T_K, p_vac_bar, H2O_SHOMATE_298_1000)
    mu_H2O_eV = mu_H2O_Jmol / NA * eV_per_J
    # Dehydration: hydrate -> anhydrous + n H2O
    # ΔG_dehyd = -ΔG_hyd = -[Delta_E + n*mu] = -Delta_E - n*mu
    DG_dehyd = -Delta_E_solids - n_H2O * mu_H2O_eV
    return DG_dehyd


# ============================================================================
# FIGURE 1: Hydration Free Energy vs Relative Humidity (298 K)
# ============================================================================
print("Generating Figure 1: ΔG_hydration vs RH at 298 K...")

fig1, ax1 = plt.subplots(figsize=(6, 4.5))

RH_range = np.linspace(0.1, 0.99, 100)
DG_8_list = []
DG_9_list = []

for RH in RH_range:
    DG_8 = hydration_free_energy(298.15, RH, 8, Delta_E_8)
    DG_9 = hydration_free_energy(298.15, RH, 9, Delta_E_9)
    DG_8_list.append(DG_8)
    DG_9_list.append(DG_9)

ax1.plot(RH_range * 100, DG_8_list, 'o-', color='#2E86AB',
         linewidth=2, markersize=0, label='$\mathregular{Na_3SbS_4{\cdot}8H_2O}$')
ax1.plot(RH_range * 100, DG_9_list, 's-', color='#A23B72',
         linewidth=2, markersize=0, label='$\mathregular{Na_3SbS_4{\cdot}9H_2O}$')

# Highlight experimental RH region (60-80%)
ax1.axvspan(60, 80, alpha=0.15, color='gray', label='Typical ambient RH')

# Mark RH=68% (experimental)
ax1.axvline(68, color='red', linestyle='--', linewidth=1, alpha=0.6, label='Experimental condition')

ax1.set_xlabel('Relative Humidity (%)', fontsize=12, fontweight='bold')
ax1.set_ylabel('$\mathregular{\Delta G_{hydration}}$ (eV)', fontsize=12, fontweight='bold')
ax1.set_title('Hydration Thermodynamics at 298 K', fontsize=13, fontweight='bold')
ax1.legend(frameon=True, fontsize=9, loc='best')
ax1.grid(True, alpha=0.3, linestyle=':')
ax1.set_xlim(10, 100)
ax1.set_ylim(-125, -50)

plt.tight_layout()
plt.savefig('../../docs/figure1_hydration_vs_RH.png', dpi=300, bbox_inches='tight')
plt.savefig('../../docs/figure1_hydration_vs_RH.pdf', bbox_inches='tight')
print("  Saved: figure1_hydration_vs_RH.png/pdf")
plt.close()


# ============================================================================
# FIGURE 2: Dehydration Free Energy vs Temperature (Vacuum)
# ============================================================================
print("Generating Figure 2: ΔG_dehydration vs T under vacuum...")

fig2, ax2 = plt.subplots(figsize=(6, 4.5))

T_range_C = np.linspace(50, 400, 100)
T_range_K = T_range_C + 273.15

# Different vacuum levels
p_vacuum_levels = [1e-6, 1e-4, 1e-2, 0.03]  # bar
colors_vac = ['#006400', '#228B22', '#90EE90', '#DDA15E']
labels_vac = ['High vacuum (10⁻⁶ bar)', 'Medium vacuum (10⁻⁴ bar)',
              'Low vacuum (10⁻² bar)', 'Ambient (0.03 bar)']

for p_vac, color, label in zip(p_vacuum_levels, colors_vac, labels_vac):
    DG_dehyd_9_list = []
    for T_K in T_range_K:
        DG_dehyd_9 = dehydration_free_energy_vacuum(T_K, p_vac, 9, Delta_E_9)
        DG_dehyd_9_list.append(DG_dehyd_9)

    ax2.plot(T_range_C, DG_dehyd_9_list, '-', color=color,
             linewidth=2, label=label)

ax2.set_xlabel('Temperature (°C)', fontsize=12, fontweight='bold')
ax2.set_ylabel('$\mathregular{\Delta G_{dehydration}}$ (eV)', fontsize=12, fontweight='bold')
ax2.set_title('Dehydration Thermodynamics', fontsize=12, fontweight='bold')
ax2.set_ylim(50, 125)
ax2.legend(frameon=True, fontsize=8, loc='best')
ax2.grid(True, alpha=0.3, linestyle=':')

plt.tight_layout()
plt.savefig('../../docs/figure2_dehydration_vs_T.png', dpi=300, bbox_inches='tight')
plt.savefig('../../docs/figure2_dehydration_vs_T.pdf', bbox_inches='tight')
print("  Saved: figure2_dehydration_vs_T.png/pdf")
plt.close()


# ============================================================================
# FIGURE 3: Phase Stability Map (T vs RH)
# ============================================================================
print("Generating Figure 3: Phase stability diagram...")

fig3, ax3 = plt.subplots(figsize=(7, 5))

# Create 2D grid
T_grid_C = np.linspace(20, 100, 50)
RH_grid = np.linspace(0.1, 0.95, 50)
T_mesh, RH_mesh = np.meshgrid(T_grid_C, RH_grid)

# Calculate ΔG for each point
DG_8_grid = np.zeros_like(T_mesh)
DG_9_grid = np.zeros_like(T_mesh)

for i in range(len(RH_grid)):
    for j in range(len(T_grid_C)):
        T_K = T_grid_C[j] + 273.15
        RH = RH_grid[i]
        DG_8_grid[i, j] = hydration_free_energy(T_K, RH, 8, Delta_E_8)
        DG_9_grid[i, j] = hydration_free_energy(T_K, RH, 9, Delta_E_9)

# Determine most stable phase at each point
# Phase code: 0 = anhydrous, 1 = 8H2O, 2 = 9H2O
phase_map = np.zeros_like(T_mesh)

for i in range(T_mesh.shape[0]):
    for j in range(T_mesh.shape[1]):
        DG_0 = 0  # anhydrous reference
        DG_8 = DG_8_grid[i, j]
        DG_9 = DG_9_grid[i, j]

        # Most stable = most negative ΔG
        phases = [DG_0, DG_8, DG_9]
        phase_map[i, j] = np.argmin(phases)

# Plot phase diagram
from matplotlib.colors import ListedColormap
colors = ['#FFF8DC', '#87CEEB', '#4682B4']  # Wheat, Sky Blue, Steel Blue
cmap = ListedColormap(colors)

contour = ax3.contourf(T_mesh, RH_mesh * 100, phase_map,
                       levels=[-0.5, 0.5, 1.5, 2.5], cmap=cmap, alpha=0.8)

# Add contour lines
ax3.contour(T_mesh, RH_mesh * 100, phase_map, levels=[0.5, 1.5],
           colors='black', linewidths=1.5, linestyles='-')

# Mark experimental recycling path
ax3.plot(25, 70, 'o', color='red', markersize=10, markeredgecolor='darkred', markeredgewidth=2)
ax3.text(27, 72, 'Hydration\n(ambient)', fontsize=9, color='red', fontweight='bold')

ax3.plot(80, 20, 's', color='darkred', markersize=10, markeredgecolor='black', markeredgewidth=2)
ax3.text(82, 22, 'Dehydration\n(vacuum)', fontsize=9,
         color='darkred', fontweight='bold')

# Arrow showing recycling pathway
ax3.annotate('', xy=(80, 20), xytext=(25, 70),
            arrowprops=dict(arrowstyle='->', lw=2.5, color='black', alpha=0.6))

ax3.set_xlabel('Temperature (°C)', fontsize=12, fontweight='bold')
ax3.set_ylabel('Relative Humidity (%)', fontsize=12, fontweight='bold')
ax3.set_title('Phase Stability Diagram', fontsize=13, fontweight='bold')

# Add legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=colors[0], label='$\mathregular{Na_3SbS_4}$ (anhydrous)'),
                   Patch(facecolor=colors[1], label='$\mathregular{Na_3SbS_4{\cdot}8H_2O}$'),
                   Patch(facecolor=colors[2], label='$\mathregular{Na_3SbS_4{\cdot}9H_2O}$')]
ax3.legend(handles=legend_elements, loc='upper right', frameon=True, fontsize=9)

ax3.set_xlim(20, 100)
ax3.set_ylim(10, 95)
ax3.grid(True, alpha=0.2, linestyle=':')

plt.tight_layout()
plt.savefig('../../docs/figure3_phase_diagram.png', dpi=300, bbox_inches='tight')
plt.savefig('../../docs/figure3_phase_diagram.pdf', bbox_inches='tight')
print("  Saved: figure3_phase_diagram.png/pdf")
plt.close()


# ============================================================================
# FIGURE 4: Energy Comparison (8H2O vs 9H2O)
# ============================================================================
print("Generating Figure 4: Energy level diagram...")

fig4, ax4 = plt.subplots(figsize=(6, 5))

# Energy levels (relative to anhydrous)
E_anhydrous = 0
E_8H2O_rel = Delta_E_8
E_9H2O_rel = Delta_E_9

# At 298 K, RH=68%
RH_ref = 0.68
DG_8_ambient = hydration_free_energy(298.15, RH_ref, 8, Delta_E_8)
DG_9_ambient = hydration_free_energy(298.15, RH_ref, 9, Delta_E_9)

# Bar positions
positions = [0, 1, 2]
labels = ['Na₃SbS₄\n(anhydrous)', 'Na₃SbS₄·8H₂O', 'Na₃SbS₄·9H₂O']
energies_solid = [E_anhydrous, E_8H2O_rel, E_9H2O_rel]
energies_free = [0, DG_8_ambient, DG_9_ambient]

# Plot solid-state energies
bars1 = ax4.bar(positions, energies_solid, width=0.35,
               label='Solid-state energy (ΔE)', color='#2E86AB', alpha=0.7)

# Plot free energies at ambient conditions
bars2 = ax4.bar([p + 0.4 for p in positions], energies_free, width=0.35,
               label='Free energy at 298K, RH=68%', color='#A23B72', alpha=0.7)

ax4.set_ylabel('Energy (eV)', fontsize=12, fontweight='bold')
ax4.set_title('Energetics of Hydration', fontsize=13, fontweight='bold')
ax4.set_xticks([p + 0.2 for p in positions])
ax4.set_xticklabels(labels, fontsize=10)
ax4.legend(frameon=True, fontsize=9)
ax4.axhline(0, color='black', linestyle='-', linewidth=0.8)
ax4.grid(True, axis='y', alpha=0.3, linestyle=':')

# Add energy difference annotations
ax4.annotate('', xy=(1, E_8H2O_rel), xytext=(0, 0),
            arrowprops=dict(arrowstyle='<->', lw=1.5, color='gray'))
ax4.text(0.5, E_8H2O_rel/2, f'{E_8H2O_rel:.1f} eV',
         fontsize=9, ha='center', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

plt.tight_layout()
plt.savefig('../../docs/figure4_energy_levels.png', dpi=300, bbox_inches='tight')
plt.savefig('../../docs/figure4_energy_levels.pdf', bbox_inches='tight')
print("  Saved: figure4_energy_levels.png/pdf")
plt.close()


# ============================================================================
# FIGURE 5: Chemical Potential of Water vs Temperature
# ============================================================================
print("Generating Figure 5: μ(H₂O) vs T at different pressures...")

fig5, ax5 = plt.subplots(figsize=(6, 4.5))

T_range_K = np.linspace(273, 673, 100)
T_range_C = T_range_K - 273.15

# Different water pressures
p_H2O_levels = [1e-6, 1e-4, 1e-2, 0.031]  # bar (last is ~RH=100% at 25°C)
colors_p = ['#006400', '#228B22', '#90EE90', '#4169E1']
labels_p = ['10⁻⁶ bar (high vacuum)', '10⁻⁴ bar (medium vacuum)',
            '10⁻² bar (low vacuum)', '0.031 bar (saturated at 25°C)']

for p, color, label in zip(p_H2O_levels, colors_p, labels_p):
    mu_list = []
    for T_K in T_range_K:
        mu_Jmol = mu_ideal_gas(T_K, p, H2O_SHOMATE_298_1000)
        mu_eV = mu_Jmol / NA * eV_per_J
        mu_list.append(mu_eV)

    ax5.plot(T_range_C, mu_list, '-', color=color, linewidth=2, label=label)

ax5.set_xlabel('Temperature (°C)', fontsize=12, fontweight='bold')
ax5.set_ylabel('μ(H₂O) (eV per molecule)', fontsize=12, fontweight='bold')
ax5.set_title('Chemical Potential of Water Vapor', fontsize=13, fontweight='bold')
ax5.legend(frameon=True, fontsize=8, loc='lower left')
ax5.grid(True, alpha=0.3, linestyle=':')

plt.tight_layout()
plt.savefig('../../docs/figure5_water_chemical_potential.png', dpi=300, bbox_inches='tight')
plt.savefig('../../docs/figure5_water_chemical_potential.pdf', bbox_inches='tight')
print("  Saved: figure5_water_chemical_potential.png/pdf")
plt.close()


# ============================================================================
# Summary table
# ============================================================================
print("\n" + "="*70)
print("FIGURE GENERATION COMPLETE")
print("="*70)
print("\nGenerated figures:")
print("  1. figure1_hydration_vs_RH.png/pdf")
print("  2. figure2_dehydration_vs_T.png/pdf")
print("  3. figure3_phase_diagram.png/pdf")
print("  4. figure4_energy_levels.png/pdf")
print("  5. figure5_water_chemical_potential.png/pdf")
print("\nAll figures saved to: calculati../../docs/")
print("="*70)
