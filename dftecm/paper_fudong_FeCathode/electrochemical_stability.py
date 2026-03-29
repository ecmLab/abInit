#!/usr/bin/env python3
"""
Chemical and electrochemical stability analysis for Fe-halide cathodes
with halide solid electrolytes.

Generates two types of plots (per cathode):
  1. Chemical reaction energy vs molar fraction (closed system, no Li reservoir)
  2. Electrochemical reaction energy vs molar fraction, curves colored by voltage
     (open Li system, voltage range V_MIN–V_MAX vs Li/Li+)

Rows computed here:
  Row 3: FeF3  + 5 electrolytes
  Row 4: FeCl3 + 5 electrolytes
"""

from mp_api.client import MPRester
from pymatgen.analysis.interface_reactions import (
    GrandPotentialInterfacialReactivity,
    InterfacialReactivity,
)
from pymatgen.analysis.phase_diagram import GrandPotentialPhaseDiagram, PhaseDiagram
from pymatgen.core import Composition, Element
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import pandas as pd

# ── Settings ──────────────────────────────────────────────────────────────────
API_KEY   = "H9GD9SySZmkHxPyxI3T5ylATO9Z8XFAd"   # update if expired
OPEN_EL   = "Li"
V_MIN     = 1.0    # V vs Li/Li+
V_MAX     = 4.2    # V vs Li/Li+
N_VOLTAGES = 321   # 0.01 V step over 1.0–4.2 V

# Cathodes (rows 3 & 4; row 2 = 1:1 molar FeF3+FeCl3 → Fe2F3Cl3)
CATHODES = {
    "FeF3":       "FeF3",
    "FeCl3":      "FeCl3",
    "FeF3-FeCl3": "Fe2F3Cl3",
}

# Electrolytes (columns)
ELECTROLYTES = {
    r"Li$_3$YCl$_6$":                       "Li3YCl6",
    r"Li$_{2.4}$Y$_{0.4}$Zr$_{0.6}$Cl$_6$": "Li2.4Y0.4Zr0.6Cl6",
    r"Li$_3$ScCl$_6$":                      "Li3ScCl6",
    r"Li$_{2.375}$Sc$_{0.375}$Zr$_{0.625}$Cl$_6$": "Li2.375Sc0.375Zr0.625Cl6",
    r"Li$_3$InCl$_6$":                      "Li3InCl6",
}

CMAP     = "coolwarm"     # blue = low V, red = high V
FIGSIZE  = (4.5, 4)       # per-subplot size

# ── Helpers ───────────────────────────────────────────────────────────────────

def atom_to_molar(x_atom, n1, n2):
    """Convert atom fraction of comp1 → molar fraction of comp1."""
    if x_atom <= 0.0:
        return 0.0
    if x_atom >= 1.0:
        return 1.0
    m1 = x_atom / n1
    m2 = (1.0 - x_atom) / n2
    return m1 / (m1 + m2)


def build_iface(comp1, comp2, pd_obj, entries, chempots=None):
    if chempots is None:
        return InterfacialReactivity(
            comp1, comp2, pd_obj,
            norm=True, use_hull_energy=False,
        )
    gpd = GrandPotentialPhaseDiagram(entries, chempots)
    return GrandPotentialInterfacialReactivity(
        comp1, comp2, gpd,
        pd_non_grand=pd_obj,
        include_no_mixing_energy=True,
        norm=True, use_hull_energy=False,
    )


def reaction_curve(comp1, comp2, pd_obj, entries, chempots=None):
    """Return (x_molar_list, energy_list, kink_detail_list)."""
    n1, n2 = comp1.num_atoms, comp2.num_atoms
    iface = build_iface(comp1, comp2, pd_obj, entries, chempots)
    xs, es, details = [], [], []
    for _, x_atom, energy, rxn, _kJ in iface.get_kinks():
        xm = atom_to_molar(x_atom, n1, n2)
        xs.append(xm)
        es.append(energy)
        details.append((xm, energy, str(rxn)))
    return xs, es, details


def parse_products(rxn_str):
    """Extract product phase names (no coefficients) from a reaction string."""
    if "->" in rxn_str:
        rhs = rxn_str.split("->")[1].strip()
    else:
        return rxn_str
    phases = []
    for term in rhs.split(" + "):
        tokens = term.strip().split()
        if not tokens:
            continue
        try:
            float(tokens[0])          # coefficient → skip it
            phases.append(" ".join(tokens[1:]))
        except ValueError:
            phases.append(" ".join(tokens))
    return "\n".join(phases)


def cluster_positions(xs, tol=0.02):
    """Cluster close x values and return their means."""
    if not xs:
        return []
    xs = sorted(xs)
    clusters, cur = [], [xs[0]]
    for x in xs[1:]:
        if x - cur[-1] < tol:
            cur.append(x)
        else:
            clusters.append(float(np.mean(cur)))
            cur = [x]
    clusters.append(float(np.mean(cur)))
    return clusters


# ── Plotting ──────────────────────────────────────────────────────────────────

ELEC_COLORS  = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]
ELEC_MARKERS = ["o", "s", "^", "D", "v"]

def make_plots(cathode_name, cathode_formula):
    comp1    = Composition(cathode_formula)
    n_elec   = len(ELECTROLYTES)
    voltages = np.linspace(V_MIN, V_MAX, N_VOLTAGES)
    norm     = Normalize(vmin=V_MIN, vmax=V_MAX)
    cmap     = plt.colormaps[CMAP]
    xlabel   = f"Molar fraction, {cathode_name}/(cathode+SE)"

    # Chemical: single combined axis; Electrochemical: one subplot per electrolyte
    fig_chem, ax_chem = plt.subplots(figsize=(6, 4.5))
    fig_elec, axes_elec = plt.subplots(
        1, n_elec,
        figsize=(FIGSIZE[0] * n_elec + 1, FIGSIZE[1]),
        sharey=True,
    )

    for i, (elec_label, elec_formula) in enumerate(ELECTROLYTES.items()):
        comp2    = Composition(elec_formula)
        elements = list({e.symbol for e in comp1.elements + comp2.elements} | {OPEN_EL})

        print(f"  [{i+1}/{n_elec}] {cathode_name} | {elec_formula}")
        print(f"    Elements: {sorted(elements)}", flush=True)

        # Fetch MP entries
        with MPRester(API_KEY) as mpr:
            entries = mpr.get_entries_in_chemsys(elements)
        print(f"    {len(entries)} entries retrieved.")

        pd_obj    = PhaseDiagram(entries)
        mu_li_ref = pd_obj.get_transition_chempots(Element(OPEN_EL))[0]

        color  = ELEC_COLORS[i]
        marker = ELEC_MARKERS[i]
        ax_e   = axes_elec[i]

        # ── Chemical stability (closed, no voltage) — all on one axis ────────
        xs, es, _ = reaction_curve(comp1, comp2, pd_obj, entries, chempots=None)
        ax_chem.plot(xs, es, "--", color=color, lw=1.5, zorder=2, label=elec_label)
        ax_chem.plot(xs, es, marker, color=color, ms=8, zorder=3,
                     markeredgecolor="white", markeredgewidth=0.8)

        # ── Electrochemical stability (open Li, voltage sweep) ───────────────
        for voltage in voltages:
            chempots = {OPEN_EL: -voltage + mu_li_ref}
            xs_e, es_e, details = reaction_curve(
                comp1, comp2, pd_obj, entries, chempots=chempots
            )
            color_v = cmap(norm(voltage))
            ax_e.plot(xs_e, es_e, "-", color=color_v, lw=0.8, alpha=0.85)
            # Mark kink positions with small green hollow circles
            ax_e.plot(xs_e, es_e, "o", ms=2, zorder=4,
                      markerfacecolor="none", markeredgecolor="green",
                      markeredgewidth=0.5)

        ax_e.axhline(0, color="k", lw=0.6, ls="--", alpha=0.5)
        ax_e.set_title(elec_label, fontsize=10)
        ax_e.set_xlabel(xlabel, fontsize=8)
        ax_e.set_xlim(0, 1)
        ax_e.tick_params(labelsize=8)

    # ── Chemical figure formatting ────────────────────────────────────────────
    ax_chem.axhline(0, color="k", lw=0.6, ls="--", alpha=0.5)
    ax_chem.set_xlabel(xlabel, fontsize=10)
    ax_chem.set_ylabel("Chemical reaction energy (eV/atom)", fontsize=10)
    ax_chem.set_xlim(0, 1)
    ax_chem.legend(fontsize=8, loc="upper right", framealpha=0.9,
                   bbox_to_anchor=(1.0, 1.0))
    ax_chem.tick_params(labelsize=9)
    fig_chem.suptitle(
        f"Chemical reaction energy — {cathode_name} vs halide electrolytes",
        fontsize=11,
    )
    fig_chem.tight_layout()
    out_chem = f"chemical_stability_{cathode_name}.png"
    fig_chem.savefig(out_chem, dpi=150, bbox_inches="tight")
    print(f"\n  → Saved: {out_chem}")

    # ── Electrochemical figure formatting ─────────────────────────────────────
    axes_elec[0].set_ylabel("Electrochemical reaction energy (eV/atom)", fontsize=9)
    fig_elec.suptitle(
        f"Electrochemical reaction energy — {cathode_name} vs halide electrolytes\n"
        f"(voltage {V_MIN}–{V_MAX} V vs Li/Li⁺, color: low→high V)",
        fontsize=11, y=1.01,
    )

    sm = ScalarMappable(cmap=CMAP, norm=norm)
    sm.set_array([])
    cbar = fig_elec.colorbar(sm, ax=axes_elec[-1], shrink=0.85, pad=0.02)
    cbar.set_label("Voltage (V vs Li/Li⁺)", fontsize=9)

    fig_elec.subplots_adjust(right=0.88, wspace=0.05)
    out_elec = f"electrochemical_stability_{cathode_name}.png"
    fig_elec.savefig(out_elec, dpi=150, bbox_inches="tight")
    print(f"  → Saved: {out_elec}")

    plt.close("all")


# ── Table output ──────────────────────────────────────────────────────────────

# Coarser voltage steps for the table (0.1 V)
V_TABLE_STEP = 0.1

def make_table(cathode_name, cathode_formula):
    """
    Output a CSV table:
      Voltage (V vs Li/Li+) | Electrolyte |
      Molar ratio cathode/(cathode+SE) | Reaction | Erxn (eV/atom)
    Rows with Erxn = 0 are excluded. One file per cathode, all 5 electrolytes combined.
    """
    comp1 = Composition(cathode_formula)
    voltages_table = np.arange(V_MIN, V_MAX + 1e-9, V_TABLE_STEP)
    rows  = []

    for elec_label, elec_formula in ELECTROLYTES.items():
        comp2    = Composition(elec_formula)
        elements = list({e.symbol for e in comp1.elements + comp2.elements} | {OPEN_EL})

        print(f"  Table: {cathode_name} | {elec_formula}", flush=True)
        with MPRester(API_KEY) as mpr:
            entries = mpr.get_entries_in_chemsys(elements)

        pd_obj    = PhaseDiagram(entries)
        mu_li_ref = pd_obj.get_transition_chempots(Element(OPEN_EL))[0]

        for voltage in voltages_table:
            chempots = {OPEN_EL: -voltage + mu_li_ref}
            _, _, details = reaction_curve(comp1, comp2, pd_obj, entries,
                                           chempots=chempots)
            for xm, energy, rxn_str in details:
                rows.append({
                    "Voltage (V vs Li/Li+)":           round(float(voltage), 2),
                    "Electrolyte":                      elec_formula,
                    "Molar ratio cathode/(cathode+SE)": f"{xm:.4f}",
                    "Reaction":                         rxn_str,
                    "Erxn (eV/atom)":                   round(energy, 4),
                })

    df = pd.DataFrame(rows, columns=[
        "Voltage (V vs Li/Li+)", "Electrolyte",
        "Molar ratio cathode/(cathode+SE)", "Reaction",
        "Erxn (eV/atom)",
    ])
    df = df[df["Erxn (eV/atom)"] != 0.0].reset_index(drop=True)

    out_csv = f"stability_table_{cathode_name}.csv"
    df.to_csv(out_csv, index=False)
    print(f"  → Saved: {out_csv}  ({len(df)} rows)")
    return df


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    for cathode_name, cathode_formula in CATHODES.items():
        print(f"\n{'='*65}")
        print(f"Cathode: {cathode_name}  ({cathode_formula})")
        make_plots(cathode_name, cathode_formula)
        make_table(cathode_name, cathode_formula)

    print("\nAll done.")


if __name__ == "__main__":
    main()
