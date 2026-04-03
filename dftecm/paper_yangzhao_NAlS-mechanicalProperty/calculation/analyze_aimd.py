"""
AIMD diffusivity analysis: XDATCAR → MSD → D → σ → Arrhenius
Produces:
  docs/NaXS_aimd_msd.png      — MSD vs time for all jobs
  docs/NaXS_aimd_arrhenius.png — Arrhenius plot + room-T extrapolation

Run from calculation/:  python analyze_aimd.py
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from pathlib import Path

CALC_DIR = Path(__file__).parent
DOCS_DIR = CALC_DIR.parent / "docs"

# Physical constants (SI)
E_CHARGE  = 1.60218e-19   # C
K_B       = 1.38065e-23   # J/K
N_A       = 6.02214e23

POTIM     = 2.0e-15        # 2 fs in seconds
EQUIL_FRAC = 0.2           # skip first 20% as equilibration

PHASES = [
    {"dir": "mp-560538_Na3AlS3",    "label": r"Na$_3$AlS$_3$",  "color": "#E8A000"},
    {"dir": "user_Na5AlS4_Na5AlS4", "label": r"Na$_5$AlS$_4$",  "color": "#C07000"},
]
TEMPS = [600, 800, 1000]


# ── XDATCAR reader ────────────────────────────────────────────────────────
def read_xdatcar(path):
    """Return lattice (3×3 Å), element list, Na indices, positions array (N_steps, N_atoms, 3) frac."""
    lines = Path(path).read_text().splitlines()
    scale   = float(lines[1].strip())
    lattice = np.array([l.split() for l in lines[2:5]], dtype=float) * scale
    elements = lines[5].split()
    counts   = list(map(int, lines[6].split()))
    n_atoms  = sum(counts)

    # Build per-atom element list
    elem_list = []
    for el, cnt in zip(elements, counts):
        elem_list.extend([el] * cnt)
    na_idx = [i for i, e in enumerate(elem_list) if e == "Na"]

    # Read frames
    frames = []
    i = 7
    while i < len(lines):
        if lines[i].strip().startswith("Direct"):
            frame = []
            for j in range(1, n_atoms + 1):
                frame.append(list(map(float, lines[i + j].split()[:3])))
            frames.append(frame)
            i += n_atoms + 1
        else:
            i += 1

    positions = np.array(frames)   # (N_steps, N_atoms, 3) fractional
    return lattice, na_idx, positions


# ── MSD with PBC unwrapping ───────────────────────────────────────────────
def compute_msd(positions, na_idx, lattice, equil_frac=0.2):
    """
    Sliding-window MSD for Na atoms, in Å².
    Returns (lag_times_s, msd_A2).
    """
    frac = positions[:, na_idx, :]   # (N_steps, N_Na, 3)
    N_steps, N_Na, _ = frac.shape

    # Unwrap fractional coordinates (correct for PBC jumps > 0.5)
    unwrapped = frac.copy()
    for t in range(1, N_steps):
        delta = unwrapped[t] - unwrapped[t - 1]
        delta -= np.round(delta)
        unwrapped[t] = unwrapped[t - 1] + delta

    # Convert to Cartesian (Å)
    cart = unwrapped @ lattice   # (N_steps, N_Na, 3)

    # Skip equilibration
    t_eq = int(N_steps * equil_frac)
    cart = cart[t_eq:]
    N_prod = len(cart)

    # Sliding-window MSD: average over all (origin, lag) pairs
    max_lag = N_prod // 2
    msd = np.zeros(max_lag)
    counts = np.zeros(max_lag, dtype=int)
    for origin in range(N_prod - max_lag):
        for lag in range(1, max_lag):
            disp = cart[origin + lag] - cart[origin]   # (N_Na, 3)
            msd[lag] += np.mean(np.sum(disp**2, axis=1))
            counts[lag] += 1

    msd[1:] /= counts[1:]
    lag_times = np.arange(max_lag) * POTIM   # seconds
    return lag_times[1:], msd[1:]


# ── Faster block-averaged MSD ─────────────────────────────────────────────
def compute_msd_fast(positions, na_idx, lattice, equil_frac=0.2, n_origins=200):
    """Faster MSD using random subset of time origins."""
    frac = positions[:, na_idx, :]
    N_steps, N_Na, _ = frac.shape

    unwrapped = frac.copy()
    for t in range(1, N_steps):
        delta = unwrapped[t] - unwrapped[t - 1]
        delta -= np.round(delta)
        unwrapped[t] = unwrapped[t - 1] + delta

    cart = unwrapped @ lattice

    t_eq    = int(N_steps * equil_frac)
    cart    = cart[t_eq:]
    N_prod  = len(cart)
    max_lag = N_prod // 2

    # Use evenly spaced origins
    origins = np.linspace(0, N_prod - max_lag - 1, n_origins, dtype=int)
    msd_sum = np.zeros(max_lag)
    for o in origins:
        disp = cart[o:o + max_lag] - cart[o]   # (max_lag, N_Na, 3)
        msd_sum += np.mean(np.sum(disp**2, axis=2), axis=1)
    msd = msd_sum / len(origins)

    lag_times = np.arange(max_lag) * POTIM
    return lag_times[1:], msd[1:]


# ── Diffusivity from MSD linear fit ──────────────────────────────────────
def fit_diffusivity(lag_times, msd, fit_frac=(0.2, 0.8)):
    """Linear fit MSD = 6D*t in the diffusive regime. Returns D in m²/s."""
    n = len(lag_times)
    i0 = int(n * fit_frac[0])
    i1 = int(n * fit_frac[1])
    slope, intercept, r, p, se = linregress(lag_times[i0:i1], msd[i0:i1])
    D_m2s = slope / 6.0 * 1e-20   # Å²/s → m²/s
    return D_m2s, r**2


# ── Nernst-Einstein conductivity ─────────────────────────────────────────
def nernst_einstein(D_m2s, n_Na, volume_A3, T):
    """σ in S/cm from D (m²/s), Na count, cell volume (Å³), T (K)."""
    vol_m3 = volume_A3 * 1e-30
    n_density = n_Na / vol_m3        # Na per m³
    sigma_SI  = (n_density * E_CHARGE**2 * D_m2s) / (K_B * T)   # S/m
    return sigma_SI * 1e-2           # S/m → S/cm


# ── Main loop ─────────────────────────────────────────────────────────────
results = {}   # {phase_label: {T: {"D": ..., "sigma": ..., "r2": ...}}}

fig_msd, axes_msd = plt.subplots(2, 3, figsize=(14, 7), sharey=False)
fig_msd.suptitle("Na MSD vs time (AIMD)", fontsize=12)

for row, phase in enumerate(PHASES):
    results[phase["label"]] = {}
    for col, T in enumerate(TEMPS):
        ax = axes_msd[row][col]
        xdat = CALC_DIR / phase["dir"] / "aimd" / f"{T}K" / "XDATCAR"
        print(f"Reading {phase['dir']} {T}K ...", flush=True)

        lattice, na_idx, positions = read_xdatcar(xdat)
        vol_A3  = abs(np.linalg.det(lattice))
        n_Na    = len(na_idx)
        N_steps = len(positions)

        lag_times, msd = compute_msd_fast(positions, na_idx, lattice,
                                          equil_frac=EQUIL_FRAC)
        D, r2 = fit_diffusivity(lag_times, msd)
        sigma = nernst_einstein(D, n_Na, vol_A3, T)

        # Mark non-diffusive runs: D < 0 or R² < 0.5 → unreliable
        reliable = (D > 0) and (r2 >= 0.5)
        results[phase["label"]][T] = {"D": D, "sigma": sigma, "r2": r2,
                                      "n_steps": N_steps, "reliable": reliable}

        # MSD plot
        t_ps = lag_times * 1e12
        ax.plot(t_ps, msd, color=phase["color"], lw=1.5)
        # linear fit overlay
        fit_i0 = int(len(lag_times) * 0.2)
        fit_i1 = int(len(lag_times) * 0.8)
        slope_fit = (D * 6e20)   # Å²/s
        t_fit = lag_times[fit_i0:fit_i1]
        msd_fit = slope_fit * t_fit + (msd[fit_i0] - slope_fit * lag_times[fit_i0])
        ax.plot(t_fit * 1e12, msd_fit, "k--", lw=1.2, alpha=0.7)

        ax.set_title(f"{phase['label']}  {T} K", fontsize=9)
        ax.set_xlabel("Time (ps)", fontsize=8)
        ax.set_ylabel("MSD (Å²)", fontsize=8)
        if reliable:
            info = f"D = {D:.2e} m²/s\nσ = {sigma:.2e} S/cm\nR² = {r2:.3f}"
            clr_txt = phase["color"]
        else:
            info = f"non-diffusive\n(R² = {r2:.3f})"
            clr_txt = "gray"
        ax.text(0.97, 0.05, info, transform=ax.transAxes, fontsize=7.5,
                ha="right", color=clr_txt, bbox=dict(fc="white", alpha=0.7, pad=2))
        ax.grid(True, linestyle="--", alpha=0.3)
        if not reliable:
            ax.set_facecolor("#F8F8F8")

        print(f"  steps={N_steps}  n_Na={n_Na}  D={D:.3e} m²/s  "
              f"σ={sigma:.3e} S/cm  R²={r2:.3f}  reliable={reliable}")

plt.tight_layout()
fig_msd.savefig(DOCS_DIR / "NaXS_aimd_msd.png", dpi=150, bbox_inches="tight")
print(f"\nSaved → NaXS_aimd_msd.png")
plt.close(fig_msd)


# ── Arrhenius fit & plot ───────────────────────────────────────────────────
fig_arr, ax_arr = plt.subplots(figsize=(6, 5))

T_RT   = 298.15   # room temperature K
inv_T  = 1000 / np.array(TEMPS)   # 1000/T for x-axis

for phase in PHASES:
    clr = phase["color"]
    lbl = phase["label"]

    # Only use reliable points for Arrhenius fit
    rel_T   = [T for T in TEMPS if results[lbl][T]["reliable"]]
    rel_s   = np.array([results[lbl][T]["sigma"] for T in rel_T])
    all_s   = np.array([results[lbl][T]["sigma"] for T in TEMPS])
    rel_inv = 1000 / np.array(rel_T)
    all_inv = inv_T

    # Plot all points (open = unreliable, filled = reliable)
    for T, s in zip(TEMPS, all_s):
        reliable = results[lbl][T]["reliable"]
        ax_arr.semilogy(1000/T, abs(s),
                        "o" if reliable else "^",
                        color=clr, ms=9 if reliable else 7,
                        alpha=1.0 if reliable else 0.35,
                        mec="black" if reliable else clr, mew=0.8)

    print(f"\n{lbl}  (reliable points: {rel_T})")
    for T in TEMPS:
        r = results[lbl][T]
        print(f"  {T} K:  D={r['D']:.3e} m²/s  σ={r['sigma']:.3e} S/cm  "
              f"R²={r['r2']:.3f}  reliable={r['reliable']}")

    if len(rel_T) < 2:
        print(f"  Insufficient reliable points for Arrhenius fit.")
        ax_arr.plot([], [], color=clr, label=f"{lbl}  (insufficient data)")
        continue

    y_arr = np.log(rel_s * np.array(rel_T, dtype=float))
    slope_arr, intercept_arr, _, _, _ = linregress(1 / np.array(rel_T), y_arr)
    Ea_eV = -slope_arr * K_B / E_CHARGE

    # Extrapolate to RT
    sigma_RT = np.exp(slope_arr / T_RT + intercept_arr) / T_RT

    # Fit line from max reliable T down to RT
    inv_T_ext = np.linspace(1000 / T_RT, rel_inv.max() * 1.02, 120)
    sigma_ext = np.exp(slope_arr * (inv_T_ext / 1000) + intercept_arr) / (1000 / inv_T_ext)
    ax_arr.semilogy(inv_T_ext, sigma_ext, "--", color=clr, lw=1.5, alpha=0.75)

    ax_arr.axvline(1000 / T_RT, color="gray", lw=0.8, ls=":", alpha=0.5)
    ax_arr.scatter([1000 / T_RT], [sigma_RT], color=clr, marker="*", s=200,
                   zorder=6, edgecolors="black", linewidths=0.6)

    ax_arr.plot([], [], color=clr,
                label=f"{lbl}  $E_a$={Ea_eV:.2f} eV  "
                      f"$\\sigma_{{RT}}$={sigma_RT:.1e} S/cm")
    print(f"  Ea (Arrhenius) = {Ea_eV:.3f} eV")
    print(f"  σ_RT (300 K)   = {sigma_RT:.3e} S/cm")

ax_arr.set_xlabel("1000 / T  (K$^{-1}$)", fontsize=11)
ax_arr.set_ylabel("Ionic conductivity $\\sigma$ (S/cm)", fontsize=11)
ax_arr.legend(fontsize=8.5, framealpha=0.9, loc="lower left")
ax_arr.grid(True, which="both", linestyle="--", alpha=0.3)
ax_arr.text(1000 / T_RT + 0.02, 1e-4, "300 K", fontsize=8, color="gray")
# Note on non-diffusive points
ax_arr.text(0.98, 0.98,
            "Filled circles: diffusive (R²≥0.5)\nTriangles: non-diffusive",
            transform=ax_arr.transAxes, fontsize=7.5, ha="right", va="top",
            color="gray", style="italic")
plt.tight_layout()
fig_arr.savefig(DOCS_DIR / "NaXS_aimd_arrhenius.png", dpi=200, bbox_inches="tight")
print(f"\nSaved → NaXS_aimd_arrhenius.png")
plt.close(fig_arr)
