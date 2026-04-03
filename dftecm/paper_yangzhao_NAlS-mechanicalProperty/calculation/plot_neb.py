"""
NEB analysis for hop1_v2 (9 intermediate images, pre-relaxed endpoints):
  (1) NaXS_neb_profiles.png  — all three energy profiles in one panel
  (2) NaXS_neb_structures.png — ASE structure views of initial and saddle images

Run from calculation/:  python plot_neb.py
"""

import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
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
        "hop":      (0, 17),
        "hop_dist": 3.225,
        "color":    "#E8A000",
    },
    {
        "dir":      "user_Na5AlS4_Na5AlS4",
        "label":    r"Na$_5$AlS$_4$",
        "hop":      (16, 21),
        "hop_dist": 3.295,
        "color":    "#C07000",
    },
    {
        "dir":      "JVASP-12818_Na5InS4",
        "label":    r"Na$_5$InS$_4$",
        "hop":      (0, 7),
        "hop_dist": 3.328,
        "color":    "#00BCD4",
    },
]

N_IMAGES = 9   # intermediate images in hop1_v2


# ── Energy & force reader ─────────────────────────────────────────────────
def read_outcar(path):
    """Return (last_energy_eV, max_force_eV_per_A) from an OUTCAR."""
    text = Path(path).read_text()
    energies = re.findall(r"energy\s+without\s+entropy=\s+([\S]+)", text)
    E = float(energies[-1]) if energies else None

    sections = text.split(
        "POSITION                                       TOTAL-FORCE (eV/Angst)"
    )
    max_f = None
    if len(sections) > 1:
        forces = []
        for line in sections[-1].splitlines()[2:]:
            parts = line.split()
            if len(parts) == 6:
                try:
                    fx, fy, fz = float(parts[3]), float(parts[4]), float(parts[5])
                    forces.append((fx**2 + fy**2 + fz**2) ** 0.5)
                except Exception:
                    break
            else:
                break
        max_f = max(forces) if forces else None
    return E, max_f


def load_phase(phase):
    """Load all energies and max forces for a phase from hop1_v2 + endpoint relaxes."""
    base = CALC_DIR / phase["dir"]
    ep0  = base / "neb" / "hop1" / "00_relax" / "OUTCAR"
    ep10 = base / "neb" / "hop1" / "06_relax" / "OUTCAR"
    neb  = base / "neb" / "hop1_v2"

    E0,  _ = read_outcar(ep0)
    E10, _ = read_outcar(ep10)

    img_E, img_maxF = [], []
    for i in range(1, N_IMAGES + 1):
        E, mf = read_outcar(neb / f"{i:02d}" / "OUTCAR")
        img_E.append(E)
        img_maxF.append(mf)

    all_E   = [E0] + img_E + [E10]
    E_ref   = min(E0, E10)
    dE      = [e - E_ref for e in all_E]
    max_f   = max(img_maxF)
    return dE, max_f


# ── Figure 1: Energy profiles ─────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6, 4.5))

for p in PHASES:
    dE, max_f = load_phase(p)
    converged = max_f < 0.05

    n_total = N_IMAGES + 2          # 11 points including endpoints
    x = np.linspace(0, 1, n_total)

    cs     = CubicSpline(x, dE, bc_type="clamped")
    x_fine = np.linspace(0, 1, 500)
    y_fine = cs(x_fine)

    saddle_idx = int(np.argmax(dE))
    Ea         = dE[saddle_idx]
    y_min      = min(dE[0], dE[-1])   # lower endpoint

    clr = p["color"]
    lbl = p["label"]

    # Line style: solid if converged, dashed if not
    ls = "-" if converged else "--"
    ax.plot(x_fine, y_fine, color=clr, lw=2.5, ls=ls, zorder=3, label=lbl)

    # Image dots (exclude endpoints)
    ax.scatter(x[1:-1], dE[1:-1], color=clr, s=45, zorder=5,
               edgecolors="white", linewidths=0.8)
    ax.scatter([x[0], x[-1]], [dE[0], dE[-1]], color=clr, s=45,
               marker="s", zorder=5, edgecolors="white", linewidths=0.8)

    # Ea annotation arrow
    ax.annotate("", xy=(x[saddle_idx], y_min), xytext=(x[saddle_idx], Ea),
                arrowprops=dict(arrowstyle="<->", color=clr, lw=1.2))
    suffix = "" if converged else "*"
    ax.text(x[saddle_idx] + 0.025, (Ea + y_min) / 2,
            f"{Ea:.3f} eV{suffix}", fontsize=8.5, color=clr,
            fontweight="bold", va="center")

    print(f"{lbl}  Ea={Ea:.4f} eV  max_F={max_f:.4f} eV/Å  converged={converged}")

ax.axhline(0, color="gray", lw=0.6, ls=":", alpha=0.6)
ax.set_xlabel("Reaction coordinate (normalised)", fontsize=11)
ax.set_ylabel("Relative energy (eV)", fontsize=11)
ax.set_xlim(-0.04, 1.04)
ax.set_ylim(bottom=-0.02)
ax.grid(True, linestyle="--", alpha=0.3)
ax.legend(fontsize=10, framealpha=0.9, loc="upper right")
ax.text(0.01, 0.98,
        "Dashed: forces not fully converged (*)",
        transform=ax.transAxes, fontsize=7.5, color="gray",
        va="top", style="italic")

plt.tight_layout()
out1 = DOCS_DIR / "NaXS_neb_profiles.png"
plt.savefig(out1, dpi=200, bbox_inches="tight")
print(f"Saved → {out1.name}")
plt.close()


# ── Figure 2: Structure snapshots (initial + saddle per phase) ────────────
fig2, axes2 = plt.subplots(2, 3, figsize=(12, 7))

for col, p in enumerate(PHASES):
    neb_dir = CALC_DIR / p["dir"] / "neb" / "hop1_v2"
    idx_i, idx_j = p["hop"]

    # Saddle: middle image (img05 for 9 images, or peak image)
    dE, _ = load_phase(p)
    saddle_img = int(np.argmax(dE))   # 0-indexed over all 11; endpoints are 0 and 10

    for row, (img_id, title_sfx) in enumerate([
        (f"{1:02d}", "Initial state"),
        (f"{saddle_img:02d}", "Saddle point"),
    ]):
        # img_id 00 and 10 live in endpoint relax dirs; 01-09 in hop1_v2
        if img_id in ("00", "10"):
            ep_label = "00_relax" if img_id == "00" else "06_relax"
            poscar = CALC_DIR / p["dir"] / "neb" / "hop1" / ep_label / "CONTCAR"
        else:
            poscar = neb_dir / img_id / "POSCAR"

        ax = axes2[row][col]
        if not poscar.exists():
            ax.text(0.5, 0.5, f"{img_id}/POSCAR\nnot found",
                    ha="center", va="center", transform=ax.transAxes)
            ax.axis("off")
            continue

        atoms = ase.io.read(str(poscar), format="vasp")
        natoms = len(atoms)
        colors, radii = [], []
        for i, a in enumerate(atoms):
            if a.symbol == "Na":
                if row == 0 and i == idx_i:
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

        if row == 1:
            base_idx = [i for i in range(natoms + 1) if i != idx_j]
            if idx_i in base_idx:
                new_idx = base_idx.index(idx_i)
                if new_idx < natoms:
                    colors[new_idx] = "#FF3333"
                    radii[new_idx]  = 0.55

        plot_atoms(atoms, ax, radii=radii, colors=colors, rotation="10x,10y,0z")
        ax.set_title(f"{p['label']}  —  {title_sfx}", fontsize=9.5)
        ax.axis("off")

axes2[0][0].set_ylabel("Initial", fontsize=10)
axes2[1][0].set_ylabel("Saddle", fontsize=10)

legend_items = [
    mpatches.Patch(color="#FF3333", label="Hopping Na"),
    mpatches.Patch(color="#FFD700", label="Na"),
    mpatches.Patch(color="#6699FF", label="Al"),
    mpatches.Patch(color="#00BCD4", label="In"),
    mpatches.Patch(color="#CCCC00", label="S"),
]
fig2.legend(handles=legend_items, loc="lower center", ncol=5,
            fontsize=9, framealpha=0.9, bbox_to_anchor=(0.5, 0.01))
plt.suptitle(r"Na$^+$ migration path — initial and saddle-point structures",
             fontsize=11, y=1.00)
plt.tight_layout(rect=[0, 0.06, 1, 1])
out2 = DOCS_DIR / "NaXS_neb_structures.png"
plt.savefig(out2, dpi=200, bbox_inches="tight")
print(f"Saved → {out2.name}")
plt.close()
