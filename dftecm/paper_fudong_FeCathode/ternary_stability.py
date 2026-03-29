#!/usr/bin/env python3
"""
Option C — Ternary interface stability: FeF3 + FeCl3 (1:1 molar) vs halide electrolytes.

Treats FeF3 and FeCl3 as SEPARATE phases (not a pseudo-compound).
Scans the composition line:
    c(x) = x * (0.5*FeF3 + 0.5*FeCl3)  +  (1-x) * SE,   x in [0,1]
and computes the grand-potential reaction energy at each point:

    E_rxn(x) = E_hull_grand[c(x)]
             - [a_FeF3 * E_hull_grand(FeF3)
              + a_FeCl3 * E_hull_grand(FeCl3)
              + a_SE   * E_hull_grand(SE)]

where a_i = atom fraction of component i in the mixture.

This differs from Option A (pseudo-compound Fe2F3Cl3) in that the reference
energy explicitly uses the separate FeF3 and FeCl3 hull energies, so any
FeF3↔FeCl3 cross-reactions are captured.

Outputs per electrolyte:
  - electrochemical_ternary_{elec}.png  (voltage-colored curves, 1–4.2 V)
  - chemical_ternary_all.png            (all 5 electrolytes on one axis, no voltage)
  - ternary_table.csv                   (0.1 V step, kink compositions)
"""

from mp_api.client import MPRester
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
API_KEY    = "H9GD9SySZmkHxPyxI3T5ylATO9Z8XFAd"
OPEN_EL    = "Li"
V_MIN      = 1.0
V_MAX      = 4.2
N_VOLTAGES = 321      # 0.01 V step for plots
N_SCAN     = 300      # composition scan resolution (per voltage)
V_TABLE_STEP = 0.1
SLOPE_TOL  = 1e-3     # |ΔE_rxn/Δx| change threshold to declare a kink

COMP_A  = Composition("FeF3")
COMP_B  = Composition("FeCl3")
CATHODE_LABEL = "FeF3+FeCl3 (1:1)"

ELECTROLYTES = {
    r"Li$_3$YCl$_6$":                                    "Li3YCl6",
    r"Li$_{2.4}$Y$_{0.4}$Zr$_{0.6}$Cl$_6$":             "Li2.4Y0.4Zr0.6Cl6",
    r"Li$_3$ScCl$_6$":                                   "Li3ScCl6",
    r"Li$_{2.375}$Sc$_{0.375}$Zr$_{0.625}$Cl$_6$":      "Li2.375Sc0.375Zr0.625Cl6",
    r"Li$_3$InCl$_6$":                                   "Li3InCl6",
}

CMAP         = "coolwarm"
FIGSIZE      = (4.5, 4)
ELEC_COLORS  = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]
ELEC_MARKERS = ["o", "s", "^", "D", "v"]

# ── Core computation ───────────────────────────────────────────────────────────

def build_mixed_comp(x_cath, comp_A, comp_B, comp_SE):
    """
    Return Composition of the mixture at cathode molar fraction x_cath.
    1 unit = x_cath * [0.5 f.u. FeF3 + 0.5 f.u. FeCl3]  +  (1-x_cath) * f.u. SE
    """
    n_A  = comp_A.num_atoms
    n_B  = comp_B.num_atoms
    n_SE = comp_SE.num_atoms
    d = {}
    # contribution from A (FeF3)
    a_A = x_cath * 0.5 * n_A
    for el, frac in comp_A.fractional_composition.items():
        d[str(el)] = d.get(str(el), 0.0) + a_A * frac
    # contribution from B (FeCl3)
    a_B = x_cath * 0.5 * n_B
    for el, frac in comp_B.fractional_composition.items():
        d[str(el)] = d.get(str(el), 0.0) + a_B * frac
    # contribution from SE
    a_SE = (1.0 - x_cath) * n_SE
    for el, frac in comp_SE.fractional_composition.items():
        d[str(el)] = d.get(str(el), 0.0) + a_SE * frac
    return Composition(d)


def strip_open_el(comp, open_el=OPEN_EL):
    """Return composition with the open element (Li) removed."""
    d = {str(el): amt for el, amt in comp.items() if str(el) != open_el}
    return Composition(d) if d else None


def gpd_hull_per_total_atom(comp, gpd, open_el=OPEN_EL):
    """
    Grand-potential hull energy per TOTAL atom of comp (including Li).

    The GPCD works in the Li-free subspace: its hull energy is per non-Li atom.
    We rescale: e_per_total = e_per_nonLi * (n_nonLi / n_total).

    This gives the grand-potential energy of comp normalized to the same
    per-atom basis as the original (Li-inclusive) mixture, so that the
    reference and product energies are on a consistent scale.
    """
    n_total = comp.num_atoms
    n_li    = float(comp[Element(open_el)]) if Element(open_el) in comp else 0.0
    n_nonli = n_total - n_li

    comp_no_li = strip_open_el(comp, open_el)
    if comp_no_li is None or comp_no_li.num_atoms == 0:
        return 0.0  # pure Li: reference is the reservoir

    e_per_nonli = gpd.get_hull_energy_per_atom(comp_no_li)
    return e_per_nonli * (n_nonli / n_total)


def reaction_energy(x_cath, comp_A, comp_B, comp_SE, gpd,
                    e_A_hull, e_B_hull, e_SE_hull):
    """
    E_rxn at cathode molar fraction x_cath (eV/atom of mixture, including Li).
    e_{A,B,SE}_hull are cached grand-potential hull energies per total atom.
    """
    if x_cath <= 0.0 or x_cath >= 1.0:
        return 0.0

    n_A  = comp_A.num_atoms
    n_B  = comp_B.num_atoms
    n_SE = comp_SE.num_atoms

    a_A  = x_cath * 0.5 * n_A
    a_B  = x_cath * 0.5 * n_B
    a_SE = (1.0 - x_cath) * n_SE
    total = a_A + a_B + a_SE

    comp_mix   = build_mixed_comp(x_cath, comp_A, comp_B, comp_SE)
    e_hull_mix = gpd_hull_per_total_atom(comp_mix, gpd)

    e_ref = (a_A * e_A_hull + a_B * e_B_hull + a_SE * e_SE_hull) / total
    return e_hull_mix - e_ref


def scan_curve(comp_A, comp_B, comp_SE, gpd, n_scan=N_SCAN):
    """
    Dense scan of E_rxn(x) along the composition line.
    Returns (xs, energies) numpy arrays.
    """
    # Cache reference hull energies once (per total atom, including Li)
    e_A_hull  = gpd_hull_per_total_atom(comp_A,  gpd)
    e_B_hull  = gpd_hull_per_total_atom(comp_B,  gpd)
    e_SE_hull = gpd_hull_per_total_atom(comp_SE, gpd)

    xs = np.linspace(0.0, 1.0, n_scan)
    es = np.array([
        reaction_energy(x, comp_A, comp_B, comp_SE, gpd,
                        e_A_hull, e_B_hull, e_SE_hull)
        for x in xs
    ])
    return xs, es


def find_kinks(xs, es, tol=SLOPE_TOL):
    """
    Detect kink positions (slope discontinuities) in a piecewise-linear curve.
    Returns list of (x, e) at kink positions.
    """
    if len(xs) < 3:
        return list(zip(xs, es))
    slopes = np.diff(es) / np.diff(xs)
    kinks = [0]  # always include x=0
    for i in range(1, len(slopes)):
        if abs(slopes[i] - slopes[i-1]) > tol:
            kinks.append(i)
    kinks.append(len(xs) - 1)  # always include x=1
    # deduplicate nearby kinks
    merged, last = [], -1
    for k in kinks:
        if k - last > max(1, len(xs) // 50):
            merged.append(k)
            last = k
    return [(xs[k], es[k]) for k in merged]


def gpd_decomp_str(comp, gpd):
    """Get decomposition string from GPCD using the Li-free composition."""
    comp_no_li = strip_open_el(comp)
    if comp_no_li is None:
        return "(pure Li)"
    try:
        decomp = gpd.get_decomposition(comp_no_li)
        return " + ".join(
            f"{amt:.4f} {e.name}"
            for e, amt in sorted(decomp.items(), key=lambda t: t[1], reverse=True)
        )
    except Exception:
        return "?"


def get_decomp_str(x_cath, comp_A, comp_B, comp_SE, gpd):
    """Return human-readable reaction string at composition x_cath."""
    if x_cath <= 1e-6:
        rhs = gpd_decomp_str(comp_SE, gpd)
        return f"{comp_SE.formula} -> {rhs}"

    if x_cath >= 1.0 - 1e-6:
        comp_mix = build_mixed_comp(0.9999, comp_A, comp_B, comp_SE)
        rhs = gpd_decomp_str(comp_mix, gpd)
        return f"0.5 {comp_A.formula} + 0.5 {comp_B.formula} -> {rhs}"

    comp_mix = build_mixed_comp(x_cath, comp_A, comp_B, comp_SE)
    lhs = (f"{x_cath*0.5:.4f} {comp_A.formula} + "
           f"{x_cath*0.5:.4f} {comp_B.formula} + "
           f"{1.0-x_cath:.4f} {comp_SE.formula}")
    rhs = gpd_decomp_str(comp_mix, gpd)
    return f"{lhs} -> {rhs}"


# ── Plotting ───────────────────────────────────────────────────────────────────

def make_plots():
    n_elec   = len(ELECTROLYTES)
    voltages = np.linspace(V_MIN, V_MAX, N_VOLTAGES)
    norm_c   = Normalize(vmin=V_MIN, vmax=V_MAX)
    cmap_obj = plt.colormaps[CMAP]
    xlabel   = f"Molar fraction, cathode/(cathode+SE)"

    fig_chem, ax_chem = plt.subplots(figsize=(6, 4.5))
    fig_elec, axes_elec = plt.subplots(
        1, n_elec,
        figsize=(FIGSIZE[0] * n_elec + 1, FIGSIZE[1]),
        sharey=True,
    )

    for i, (elec_label, elec_formula) in enumerate(ELECTROLYTES.items()):
        comp_SE  = Composition(elec_formula)
        elements = list(
            {e.symbol for e in COMP_A.elements + COMP_B.elements + comp_SE.elements}
            | {OPEN_EL}
        )
        print(f"  [{i+1}/{n_elec}] {elec_formula}  elements: {sorted(elements)}", flush=True)

        with MPRester(API_KEY) as mpr:
            entries = mpr.get_entries_in_chemsys(elements)
        print(f"    {len(entries)} entries retrieved.")

        pd_obj    = PhaseDiagram(entries)
        mu_li_ref = pd_obj.get_transition_chempots(Element(OPEN_EL))[0]

        color  = ELEC_COLORS[i]
        marker = ELEC_MARKERS[i]
        ax_e   = axes_elec[i]

        # ── Chemical (closed system, no voltage) ──────────────────────────────
        gpd_chem = GrandPotentialPhaseDiagram(entries, {OPEN_EL: mu_li_ref})
        xs_c, es_c = scan_curve(COMP_A, COMP_B, comp_SE, gpd_chem)
        ax_chem.plot(xs_c, es_c, "--", color=color, lw=1.5, zorder=2, label=elec_label)
        kinks_c = find_kinks(xs_c, es_c)
        xk, ek = zip(*kinks_c) if kinks_c else ([], [])
        ax_chem.plot(xk, ek, marker, color=color, ms=8, zorder=3,
                     markeredgecolor="white", markeredgewidth=0.8)

        # ── Electrochemical (open Li, voltage sweep) ──────────────────────────
        for voltage in voltages:
            chempots = {OPEN_EL: -voltage + mu_li_ref}
            gpd_v = GrandPotentialPhaseDiagram(entries, chempots)
            xs_e, es_e = scan_curve(COMP_A, COMP_B, comp_SE, gpd_v)
            color_v = cmap_obj(norm_c(voltage))
            ax_e.plot(xs_e, es_e, "-", color=color_v, lw=0.8, alpha=0.85)
            kinks_e = find_kinks(xs_e, es_e)
            xk_e, ek_e = zip(*kinks_e) if kinks_e else ([], [])
            ax_e.plot(xk_e, ek_e, "o", ms=2, zorder=4,
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
        f"Chemical reaction energy — {CATHODE_LABEL} vs halide electrolytes",
        fontsize=11,
    )
    fig_chem.tight_layout()
    fig_chem.savefig("chemical_ternary_all.png", dpi=150, bbox_inches="tight")
    print("\n  → Saved: chemical_ternary_all.png")

    # ── Electrochemical figure formatting ─────────────────────────────────────
    axes_elec[0].set_ylabel("Electrochemical reaction energy (eV/atom)", fontsize=9)
    fig_elec.suptitle(
        f"Electrochemical reaction energy — {CATHODE_LABEL} vs halide electrolytes\n"
        f"(voltage {V_MIN}–{V_MAX} V vs Li/Li⁺, color: low→high V)",
        fontsize=11, y=1.01,
    )
    sm = ScalarMappable(cmap=CMAP, norm=norm_c)
    sm.set_array([])
    cbar = fig_elec.colorbar(sm, ax=axes_elec[-1], shrink=0.85, pad=0.02)
    cbar.set_label("Voltage (V vs Li/Li⁺)", fontsize=9)
    fig_elec.subplots_adjust(right=0.88, wspace=0.05)
    fig_elec.savefig("electrochemical_ternary_all.png", dpi=150, bbox_inches="tight")
    print("  → Saved: electrochemical_ternary_all.png")

    plt.close("all")


# ── Table ─────────────────────────────────────────────────────────────────────

def make_table():
    voltages_table = np.arange(V_MIN, V_MAX + 1e-9, V_TABLE_STEP)
    rows = []

    for elec_label, elec_formula in ELECTROLYTES.items():
        comp_SE  = Composition(elec_formula)
        elements = list(
            {e.symbol for e in COMP_A.elements + COMP_B.elements + comp_SE.elements}
            | {OPEN_EL}
        )
        print(f"  Table: {elec_formula}", flush=True)
        with MPRester(API_KEY) as mpr:
            entries = mpr.get_entries_in_chemsys(elements)

        pd_obj    = PhaseDiagram(entries)
        mu_li_ref = pd_obj.get_transition_chempots(Element(OPEN_EL))[0]

        for voltage in voltages_table:
            chempots = {OPEN_EL: -voltage + mu_li_ref}
            gpd = GrandPotentialPhaseDiagram(entries, chempots)

            xs_scan, es_scan = scan_curve(COMP_A, COMP_B, comp_SE, gpd)
            kinks = find_kinks(xs_scan, es_scan)

            for xk, ek in kinks:
                if abs(ek) < 1e-6:
                    continue
                rxn_str = get_decomp_str(xk, COMP_A, COMP_B, comp_SE, gpd)
                rows.append({
                    "Voltage (V vs Li/Li+)":           round(float(voltage), 2),
                    "Electrolyte":                      elec_formula,
                    "Molar ratio cathode/(cathode+SE)": f"{xk:.4f}",
                    "Reaction":                         rxn_str,
                    "Erxn (eV/atom)":                   round(float(ek), 4),
                })

    df = pd.DataFrame(rows, columns=[
        "Voltage (V vs Li/Li+)", "Electrolyte",
        "Molar ratio cathode/(cathode+SE)", "Reaction",
        "Erxn (eV/atom)",
    ])
    df = df[df["Erxn (eV/atom)"] != 0.0].reset_index(drop=True)
    df.to_csv("ternary_table.csv", index=False)
    print(f"  → Saved: ternary_table.csv  ({len(df)} rows)")
    return df


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("=== Ternary stability: FeF3 + FeCl3 (1:1) vs halide electrolytes ===\n")
    print("--- Generating plots ---")
    make_plots()
    print("\n--- Generating table ---")
    make_table()
    print("\nAll done.")


if __name__ == "__main__":
    main()
