"""
Set up NVT-AIMD diffusivity calculations for the three most promising phases:
  mp-560538_Na3AlS3   (56 atoms → 1×2×1 supercell = 112 atoms)
  user_Na5AlS4_Na5AlS4 (80 atoms → 1×1×1, already large)
  JVASP-12818_Na5InS4  (20 atoms → 2×2×2 supercell = 160 atoms)

Three temperatures per phase: 600 K, 800 K, 1000 K
Output tree:
  <phase>/aimd/<T>K/{POSCAR, INCAR, KPOINTS, POTCAR, submit.sh}

Run from calculation/:  python setup_aimd.py
"""

import shutil
import numpy as np
from pathlib import Path

CALC_DIR = Path(__file__).parent
TEMPS    = [600, 800, 1000]   # K

# ── Phase registry ────────────────────────────────────────────────────────
PHASES = {
    "mp-560538_Na3AlS3":    {"supercell": (1, 2, 1)},
    "user_Na5AlS4_Na5AlS4": {"supercell": (1, 1, 1)},
    "JVASP-12818_Na5InS4":  {"supercell": (2, 2, 2)},
}


# ── Minimal POSCAR supercell builder (no external deps) ───────────────────
def read_poscar(path):
    """Return dict with keys: comment, scale, lattice, elements, counts,
    coord_type, positions (fractional)."""
    lines = Path(path).read_text().splitlines()
    comment  = lines[0]
    scale    = float(lines[1])
    lattice  = np.array([l.split() for l in lines[2:5]], dtype=float)
    elements = lines[5].split()
    counts   = list(map(int, lines[6].split()))
    ctype    = lines[7].strip()
    n_atoms  = sum(counts)
    pos = np.array([l.split()[:3] for l in lines[8:8 + n_atoms]], dtype=float)
    if ctype.lower().startswith("c"):
        # convert Cartesian → fractional
        cart = pos * scale
        pos  = np.linalg.solve(lattice.T * scale, cart.T).T
    return dict(comment=comment, scale=scale, lattice=lattice,
                elements=elements, counts=counts, positions=pos)


def make_supercell(p, sc):
    """Expand structure p by supercell tuple (na, nb, nc)."""
    na, nb, nc = sc
    new_lat = p["lattice"].copy()
    new_lat[0] *= na
    new_lat[1] *= nb
    new_lat[2] *= nc

    new_pos = []
    new_counts = []
    offset = 0
    for cnt in p["counts"]:
        atom_pos = p["positions"][offset:offset + cnt]
        expanded = []
        for ia in range(na):
            for ib in range(nb):
                for ic in range(nc):
                    shifted = (atom_pos + np.array([ia, ib, ic])) / np.array([na, nb, nc])
                    expanded.append(shifted)
        new_pos.extend(np.vstack(expanded))
        new_counts.append(cnt * na * nb * nc)
        offset += cnt

    return dict(comment=p["comment"], scale=p["scale"], lattice=new_lat,
                elements=p["elements"], counts=new_counts,
                positions=np.array(new_pos))


def write_poscar(p, path):
    lines = [p["comment"]]
    lines.append(f"  {p['scale']:.10f}")
    for row in p["lattice"]:
        lines.append("  " + "  ".join(f"{v:22.16f}" for v in row))
    lines.append("  " + "  ".join(p["elements"]))
    lines.append("  " + "  ".join(map(str, p["counts"])))
    lines.append("Direct")
    for pos in p["positions"]:
        lines.append("  " + "  ".join(f"{v:.16f}" for v in pos))
    Path(path).write_text("\n".join(lines) + "\n")


# ── File templates ────────────────────────────────────────────────────────
def incar(T, job_name):
    return f"""\
SYSTEM = NVT-AIMD {job_name} {T}K

! Electronic structure
ENCUT  = 400
EDIFF  = 1E-05
NELMIN = 4
PREC   = Normal
ALGO   = Fast
LREAL  = Auto

! Molecular dynamics (NVT, Nose-Hoover)
IBRION = 0
MDALGO = 2
NSW    = 10000
POTIM  = 2.0
TEBEG  = {T}
TEEND  = {T}
SMASS  = 0

! Symmetry & output
ISYM   = 0
LWAVE  = .FALSE.
LCHARG = .FALSE.
NBLOCK = 1
KBLOCK = 100

! Parallelisation
NPAR   = 4
"""


KPOINTS = """\
Gamma-point only
 0
Gamma
 1 1 1
 0 0 0
"""


def submit_sh(job_name, T):
    return f"""\
#!/bin/bash -l
#SBATCH --job-name={job_name}_{T}K
#SBATCH --account=membrane
#SBATCH --partition=tier3

#SBATCH --mail-user=slack:@qhteme

#SBATCH -t 3-00:00:00

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=500M

spack load vasp@6.3.2
srun vasp_std
"""


# ── Main ──────────────────────────────────────────────────────────────────
for phase_dir_name, cfg in PHASES.items():
    phase_dir = CALC_DIR / phase_dir_name
    contcar   = phase_dir / "relax" / "CONTCAR"
    potcar    = phase_dir / "relax" / "POTCAR"

    if not contcar.exists():
        print(f"WARN  missing CONTCAR: {phase_dir_name}")
        continue

    # Build supercell once
    sc = cfg["supercell"]
    p  = read_poscar(contcar)
    sp = make_supercell(p, sc)
    n_atoms = sum(sp["counts"])
    label   = "×".join(map(str, sc))
    print(f"{phase_dir_name}: {sum(p['counts'])} → {n_atoms} atoms ({label} supercell)")

    short_name = phase_dir_name.split("_")[0]   # mp-560538 / user / JVASP-12818

    for T in TEMPS:
        out_dir = phase_dir / "aimd" / f"{T}K"
        out_dir.mkdir(parents=True, exist_ok=True)

        write_poscar(sp, out_dir / "POSCAR")
        (out_dir / "INCAR").write_text(incar(T, phase_dir_name))
        (out_dir / "KPOINTS").write_text(KPOINTS)
        (out_dir / "submit.sh").write_text(submit_sh(short_name, T))

        if potcar.exists():
            shutil.copy(potcar, out_dir / "POTCAR")

        print(f"  → {out_dir.relative_to(CALC_DIR)}")

print("\nDone. Review INCAR NSW and POTIM before submitting.")
print("Submit with: python batch_submit_aimd.py")
