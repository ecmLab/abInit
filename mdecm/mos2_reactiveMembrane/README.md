# HCrO₄⁻ Adsorption in MoS₂ Nanotubes

This repository hosts the analysis pipeline for molecular-dynamics simulations of chromate (HCrO₄⁻) transport through negatively charged MoS₂ nanotubes. Two pore sizes are compared – 1 nm (5 ions) and 2 nm (20 ions) – each prepared at ~45 ppb chromate with a tenfold excess of NaCl and additional Na⁺ counter-ions to neutralise the lattice charge (–7×10⁻⁴ e Å⁻²). The most recent analysis focuses on adsorption kinetics defined with a 2 Å hydration-shell cutoff and produces publication-ready figures and a LaTeX report.

## Quick Start

```bash
cd analysis

# 1. Compute adsorption metrics (Δ = 2 Å)
python adsorption_kinetics.py --delta 2 r10_1 r10_2 r10_3 r10_4 r20_1 r20_2 r20_3 r20_4

# 2. Regenerate all figures (representative trajectories, displacement, adsorption)
MPLCONFIGDIR=$(pwd)/.mplcache python generate_report_figures.py

# 3. Rebuild the PDF report
cd ../analysis/results
pdflatex -interaction=nonstopmode "MD simulation.tex"
pdflatex -interaction=nonstopmode "MD simulation.tex"
```

All derived artefacts land in `analysis/results/` and `analysis/results/report/`.

## Toolkit

| Script | Purpose |
| --- | --- |
| `analysis/trajectory_utils.py` | Reads LAMMPS dumps, unwraps PBC, detects tube centre, builds per-ion COM trajectories |
| `analysis/adsorption_kinetics.py` | Computes first-passage times, survival curves, residence fractions (configurable Δ) |
| `analysis/generate_report_figures.py` | Creates representative polar/axial/radial plots, group displacement overlays, survival and metric comparisons |
| `analysis/polar_trajectory.py`, `analysis/displacement.py` | Legacy single-simulation visualisations |
| `analysis/run_all_analysis.py` | Legacy orchestration script (kept for reference) |

Key figures and tables used in the report are:

- `representative_r10_1.png`, `representative_r20_1.png`: polar-style trajectory panels highlighting ion 2 (1 nm) and ion 3 (2 nm).
- `displacement_group_comparison.png`: averaged axial/radial displacement envelopes.
- `adsorption_survival_comparison.png`: survival probability contrast.
- `adsorption_metrics_comparison.png`: normalised success fraction, mean first passage, and wall residence fraction.
- `MD simulation.pdf`: up-to-date narrative combining the above with quantitative tables.

## Current Findings (Δ = 2 Å)

- **Survival curves:** 1 nm pores exhaust the surviving population within ~1 ns, while ~10 % of ions in 2 nm pores remain unadsorbed after 5 ns.
- **Success fraction:** 0.50 (1 nm) vs 0.10 (2 nm).
- **Mean first-passage time:** ~8.7×10² ps (1 nm) vs ~1.5×10³ ps (2 nm).
- **Wall residence fraction:** 0.29 (1 nm) vs 0.061 (2 nm).
- **Representative trajectories:** axial progress slows once ions latch onto the wall; the 1 nm ion advances >40 Å before plateauing, whereas the 2 nm ion stalls near 25 Å.

The tightened adsorption criterion sharpened the contrast between the two pore sizes without altering the qualitative conclusion: geometric confinement strongly enhances HCrO₄⁻ uptake.

## Directory Guide

```
analysis/
├── adsorption_kinetics.py        # adsorption metrics (Δ configurable)
├── generate_report_figures.py    # composite figures for the report
├── displacement.py, polar_trajectory.py  # legacy per-ion plots
├── trajectory_utils.py           # dump-reader utilities
└── results/
    ├── report/                   # generated PNGs used in the PDF
    ├── adsorption/               # CSV/PNG outputs from kinetics script
    └── MD simulation.tex/pdf     # main report sources
```

See `CONVERSATION_LOG.md` for a chronological record of design decisions and outstanding follow-ups.
