"""
Crystal structure plots for the three NEB candidate phases.
Output: docs/NaXS_crystal_structures.png

Run from calculation/:  python plot_crystal_structures.py
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from pathlib import Path

import ase.io
from ase.visualize.plot import plot_atoms

CALC_DIR = Path(__file__).parent
DOCS_DIR = CALC_DIR.parent / "docs"

PHASES = [
    {
        "contcar":   "mp-560538_Na3AlS3/relax/CONTCAR",
        "label":     r"Na$_3$AlS$_3$",
        "sg":        r"$P2_1/c$  (No. 14)",
        "rotation":  "10x,5y,0z",
    },
    {
        "contcar":   "user_Na5AlS4_Na5AlS4/relax/CONTCAR",
        "label":     r"Na$_5$AlS$_4$",
        "sg":        r"$Pbca$  (No. 61)",
        "rotation":  "10x,5y,0z",
    },
    {
        "contcar":   "JVASP-12818_Na5InS4/relax/CONTCAR",
        "label":     r"Na$_5$InS$_4$",
        "sg":        r"$P2_1/m$  (No. 11)",
        "rotation":  "10x,5y,0z",
    },
]

# Atom colours and radii
ATOM_STYLE = {
    "Na": {"color": "#FFD700", "radius": 0.50},
    "Al": {"color": "#6699FF", "radius": 0.42},
    "In": {"color": "#00BCD4", "radius": 0.50},
    "S":  {"color": "#CCCC00", "radius": 0.38},
}
DEFAULT_STYLE = {"color": "gray", "radius": 0.38}


def atom_style_lists(atoms):
    colors, radii = [], []
    for a in atoms:
        s = ATOM_STYLE.get(a.symbol, DEFAULT_STYLE)
        colors.append(s["color"])
        radii.append(s["radius"])
    return colors, radii


fig, axes = plt.subplots(1, 3, figsize=(14, 5))

for ax, phase in zip(axes, PHASES):
    contcar = CALC_DIR / phase["contcar"]
    atoms = ase.io.read(str(contcar), format="vasp")

    colors, radii = atom_style_lists(atoms)
    plot_atoms(atoms, ax, radii=radii, colors=colors,
               rotation=phase["rotation"], show_unit_cell=2)

    ax.set_title(f"{phase['label']}\n{phase['sg']}", fontsize=11, pad=8)
    ax.axis("off")

# Shared legend
legend_items = [
    mpatches.Patch(facecolor="#FFD700", edgecolor="gray", linewidth=0.5, label="Na"),
    mpatches.Patch(facecolor="#6699FF", edgecolor="gray", linewidth=0.5, label="Al"),
    mpatches.Patch(facecolor="#00BCD4", edgecolor="gray", linewidth=0.5, label="In"),
    mpatches.Patch(facecolor="#CCCC00", edgecolor="gray", linewidth=0.5, label="S"),
]
fig.legend(handles=legend_items, loc="lower center", ncol=4,
           fontsize=11, framealpha=0.9, bbox_to_anchor=(0.5, 0.01))

plt.suptitle("Crystal structures of NEB candidate phases", fontsize=12, y=1.01)
plt.tight_layout(rect=[0, 0.08, 1, 1])

out = DOCS_DIR / "NaXS_crystal_structures.png"
plt.savefig(out, dpi=200, bbox_inches="tight")
print(f"Saved → {out.name}")
plt.close()
