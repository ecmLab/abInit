"""
Two-figure NEB analysis:
  (1) NaXS_neb_profiles.png  — all three energy profiles in one panel
  (2) NaXS_neb_structures.png — ASE structure views of initial and saddle images

Run from calculation/:  python plot_neb.py
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
from scipy.interpolate import CubicSpline
from pathlib import Path

import ase.io
from ase.visualize.plot import plot_atoms

CALC_DIR = Path(__file__).parent
DOCS_DIR = CALC_DIR.parent / "docs"

PHASES = [
    {
        "dir":      "mp-560538_Na3AlS3",
        "label":    r"Na$_3$AlS$_3$",
        "hop":      (0, 17),          # (idx_i, idx_j) from setup_neb output
        "hop_dist": 3.225,
        "color":    "#E8A000",
        "energies": [-222.09223164, -221.93623434, -221.85916256,
                     -221.98028192, -222.11056638],
    },
    {
        "dir":      "user_Na5AlS4_Na5AlS4",
        "label":    r"Na$_5$AlS$_4$",
        "hop":      (16, 21),
        "hop_dist": 3.295,
        "color":    "#C07000",
        "energies": [-303.84804865, -303.50731989, -303.29618751,
                     -303.50731989, -303.84804865],
    },
    {
        "dir":      "JVASP-12818_Na5InS4",
        "label":    r"Na$_5$InS$_4$",
        "hop":      (0, 7),
        "hop_dist": 3.328,
        "color":    "#00BCD4",
        "energies": [-139.19245161, -138.77453993, -138.54302647,
                     -138.80038807, -139.18313359],
    },
]

# ── Figure 1: Energy profiles (single panel, all three overlaid) ──────────
fig, ax = plt.subplots(figsize=(6, 4.5))

for p in PHASES:
    E     = np.array(p["energies"])
    E_ref = min(E[0], E[-1])
    dE    = E - E_ref

    # 7-point path (normalised x: 0→1)
    x7  = np.linspace(0, 1, 7)
    dE7 = np.array([E[0]-E_ref, dE[0], dE[1], dE[2], dE[3], dE[4], E[-1]-E_ref])

    # Clamped spline: zero first derivative at endpoints → no spurious oscillation
    cs     = CubicSpline(x7, dE7, bc_type="clamped")
    x_fine = np.linspace(0, 1, 500)
    y_fine = cs(x_fine)

    # Ea from discrete images: saddle (img03) minus lower endpoint — physically correct
    Ea     = dE[2]          # E_img03 - min(E_img01, E_img05)
    y_min  = min(dE7[0], dE7[-1])   # lower endpoint energy in relative coords
    x_peak = x7[3]          # img03 position

    clr = p["color"]
    ax.plot(x_fine, y_fine, color=clr, lw=2.5, zorder=3, label=p["label"])
    ax.scatter(x7[1:6], dE7[1:6], color=clr, s=55, zorder=5,
               edgecolors="white", linewidths=0.8)
    ax.scatter([x7[0], x7[-1]], [dE7[0], dE7[-1]], color=clr, s=55,
               marker="s", zorder=5, edgecolors="white", linewidths=0.8)

    # Arrow from saddle (max) down to the lower endpoint (min)
    ax.annotate("", xy=(x_peak, y_min), xytext=(x_peak, Ea),
                arrowprops=dict(arrowstyle="<->", color=clr, lw=1.2))
    ax.text(x_peak + 0.025, (Ea + y_min) / 2,
            f"{Ea:.3f} eV", fontsize=8.5, color=clr, fontweight="bold", va="center")

ax.axhline(0, color="gray", lw=0.6, ls=":", alpha=0.6)
ax.set_xlabel("Reaction coordinate (normalised)", fontsize=11)
ax.set_ylabel("Relative energy (eV)", fontsize=11)
ax.set_xlim(-0.04, 1.04)
ax.set_ylim(bottom=-0.06)
ax.grid(True, linestyle="--", alpha=0.3)
ax.legend(fontsize=10, framealpha=0.9, loc="upper right")
ax.text(0.01, 0.98, "preliminary (NSW = 300)", transform=ax.transAxes,
        fontsize=7.5, color="gray", va="top", style="italic")

plt.tight_layout()
out1 = DOCS_DIR / "NaXS_neb_profiles.png"
plt.savefig(out1, dpi=200, bbox_inches="tight")
print(f"Saved → {out1.name}")
plt.close()

# ── Figure 2: Structure snapshots (initial + saddle per phase) ────────────
# Element colours matching the paper scheme
ase_colors = {"Na": "#FFD700", "Al": "#F4A900", "In": "#00BCD4", "S": "#FFFF00"}

fig2, axes2 = plt.subplots(2, 3, figsize=(12, 7))

for col, p in enumerate(PHASES):
    neb_dir = CALC_DIR / p["dir"] / "neb" / "hop1"
    idx_i, idx_j = p["hop"]

    for row, (img_id, title_sfx) in enumerate([("00", "Initial state"),
                                                ("03", "Saddle point")]):
        poscar = neb_dir / img_id / "POSCAR"
        ax = axes2[row][col]

        if not poscar.exists():
            ax.text(0.5, 0.5, f"{img_id}/POSCAR\nnot found",
                    ha="center", va="center", transform=ax.transAxes)
            ax.axis("off")
            continue

        atoms = ase.io.read(str(poscar), format="vasp")

        # Colour atoms by species, highlight hopping Na in red
        natoms = len(atoms)
        colors = []
        radii  = []
        for i, a in enumerate(atoms):
            if a.symbol == "Na":
                if row == 0 and i == idx_i:   # hopping Na in initial image
                    colors.append("#FF3333"); radii.append(0.55)
                else:
                    colors.append("#FFD700"); radii.append(0.45)
            elif a.symbol == "Al":
                colors.append("#6699FF"); radii.append(0.45)
            elif a.symbol == "In":
                colors.append("#00BCD4"); radii.append(0.50)
            elif a.symbol == "S":
                colors.append("#CCCC00"); radii.append(0.40)
            else:
                colors.append("gray"); radii.append(0.40)

        # Saddle: highlight the hopping Na (it has moved relative to initial)
        if row == 1:
            # hopping atom index in vacancy structure = hop_in_base
            # after removing idx_j, the new index of idx_i:
            base_idx = [i for i in range(natoms + 1) if i != idx_j]
            if idx_i in base_idx:
                new_idx = base_idx.index(idx_i)
                if new_idx < natoms:
                    colors[new_idx] = "#FF3333"
                    radii[new_idx]  = 0.55

        plot_atoms(atoms, ax, radii=radii, colors=colors,
                   rotation=("10x,10y,0z"))

        ax.set_title(f"{p['label']}  —  {title_sfx}", fontsize=9.5)
        ax.axis("off")

# Row labels
axes2[0][0].set_ylabel("Initial", fontsize=10)
axes2[1][0].set_ylabel("Saddle", fontsize=10)

# Legend
legend_items = [
    mpatches.Patch(color="#FF3333", label="Hopping Na"),
    mpatches.Patch(color="#FFD700", label="Na"),
    mpatches.Patch(color="#6699FF", label="Al"),
    mpatches.Patch(color="#00BCD4", label="In"),
    mpatches.Patch(color="#CCCC00", label="S"),
]
fig2.legend(handles=legend_items, loc="lower center", ncol=5,
            fontsize=9, framealpha=0.9, bbox_to_anchor=(0.5, 0.01))

plt.suptitle("Na$^+$ migration path — initial and saddle-point structures",
             fontsize=11, y=1.00)
plt.tight_layout(rect=[0, 0.06, 1, 1])
out2 = DOCS_DIR / "NaXS_neb_structures.png"
plt.savefig(out2, dpi=200, bbox_inches="tight")
print(f"Saved → {out2.name}")
plt.close()
