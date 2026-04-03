"""
Set up NEB calculations for the 3D-percolation bottleneck hop in each phase.

The bottleneck hop is the longest Na-Na distance needed to first achieve
3D percolation (diffusion in all three lattice directions). This is the
rate-limiting hop for macroscopic ionic conductivity.

Verified bottleneck hops (Na-index pairs, global atom indices):
  Na3AlS3:  Na[1] – Na[8]   (global  1–  8)  3.831 Å
  Na5AlS4:  Na[10]– Na[12]  (global 10– 12)  3.659 Å
  Na5InS4:  Na[1] – Na[8]   (global  1–  8)  3.504 Å

Creates:
  <phase>/neb/hop_bottleneck/00_relax/  — endpoint A (hopper at A, vac at B)
  <phase>/neb/hop_bottleneck/10_relax/  — endpoint B (hopper at B, vac at A)

Run from calculation/:  python setup_neb_bottleneck.py
Then submit endpoints, wait, then run:  python setup_neb_bottleneck_images.py
"""

import shutil
import numpy as np
from pathlib import Path

CALC_DIR = Path(__file__).parent
N_IMAGES = 9

PHASES = [
    {
        "dir":     "mp-560538_Na3AlS3",
        "name":    "Na3AlS3",
        "idx_hop": 1,    # Na atom index (0-based among Na atoms = global atom index here)
        "idx_vac": 8,    # Na atom that becomes vacancy (destination)
        "hop_dist": 3.831,
    },
    {
        "dir":     "user_Na5AlS4_Na5AlS4",
        "name":    "Na5AlS4",
        "idx_hop": 10,
        "idx_vac": 12,
        "hop_dist": 3.659,
    },
    {
        "dir":     "JVASP-12818_Na5InS4",
        "name":    "Na5InS4",
        "idx_hop": 1,
        "idx_vac": 8,
        "hop_dist": 3.504,
    },
]

INCAR_ENDPOINT = """\
SYSTEM = endpoint_relax_bottleneck

ENCUT  = 500
EDIFF  = 1E-06
PREC   = Accurate
ALGO   = Fast
LREAL  = .FALSE.
ADDGRID = .TRUE.
ISYM   = 0

IBRION = 2
ISIF   = 2
NSW    = 300
EDIFFG = -0.02

LWAVE  = .TRUE.
LCHARG = .FALSE.
"""

KPOINTS = """\
Gamma-point only
 0
Gamma
 1 1 1
 0 0 0
"""


def submit_endpoint_sh(job_name):
    return f"""\
#!/bin/bash -l
#SBATCH --job-name={job_name}
#SBATCH --account=membrane
#SBATCH --partition=tier3
#SBATCH --mail-user=slack:@qhteme
#SBATCH -t 5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=500M

spack load vasp@6.3.2
srun vasp_std
"""


def read_poscar(path):
    lines = Path(path).read_text().splitlines()
    scale   = float(lines[1])
    lattice = np.array([l.split() for l in lines[2:5]], dtype=float) * scale
    elements = lines[5].split()
    counts   = list(map(int, lines[6].split()))
    ctype    = lines[7].strip()
    n = sum(counts)
    pos = np.array([l.split()[:3] for l in lines[8:8+n]], dtype=float)
    if ctype.lower().startswith("c"):
        pos = np.linalg.solve(lattice.T, pos.T).T
    return dict(lattice=lattice, elements=elements, counts=counts, positions=pos)


def write_poscar(p, path, comment="NEB image"):
    lines = [comment, "  1.0"]
    for row in p["lattice"]:
        lines.append("  " + "  ".join(f"{v:22.16f}" for v in row))
    lines.append("  " + "  ".join(p["elements"]))
    lines.append("  " + "  ".join(map(str, p["counts"])))
    lines.append("Direct")
    for pos in p["positions"]:
        lines.append("  " + "  ".join(f"{v:.16f}" for v in pos))
    Path(path).write_text("\n".join(lines) + "\n")


def build_vacancy_structure(p, idx_hop, idx_vac):
    """Remove idx_vac (destination), leaving vacancy there."""
    elem_list = []
    for el, cnt in zip(p["elements"], p["counts"]):
        elem_list.extend([el] * cnt)
    base_idx = [i for i in range(len(elem_list)) if i != idx_vac]
    new_counts = []
    cursor = 0
    for el, cnt in zip(p["elements"], p["counts"]):
        kept = [i for i in range(cursor, cursor + cnt) if i in base_idx]
        new_counts.append(len(kept))
        cursor += cnt
    return dict(
        lattice=p["lattice"], elements=p["elements"],
        counts=new_counts, positions=p["positions"][base_idx].copy(),
    )


for phase in PHASES:
    phase_dir = CALC_DIR / phase["dir"]
    contcar   = phase_dir / "relax" / "CONTCAR"
    potcar    = phase_dir / "relax" / "POTCAR"
    p = read_poscar(contcar)

    # Verify hop distance
    d = p["positions"][phase["idx_vac"]] - p["positions"][phase["idx_hop"]]
    d -= np.round(d)
    dist = np.linalg.norm(d @ p["lattice"])
    print(f"{phase['name']}: hop Na[{phase['idx_hop']}]→Na[{phase['idx_vac']}]  "
          f"actual={dist:.4f} Å  expected≈{phase['hop_dist']} Å")

    hop_dir = phase_dir / "neb" / "hop_bottleneck"
    hop_dir.mkdir(parents=True, exist_ok=True)

    # Endpoint 00: hopper at A, vacancy at B
    ep0 = build_vacancy_structure(p, phase["idx_hop"], phase["idx_vac"])
    ep0_dir = hop_dir / "00_relax"
    ep0_dir.mkdir(exist_ok=True)
    write_poscar(ep0, ep0_dir / "POSCAR",
                 comment=f"{phase['name']} bottleneck ep00")
    (ep0_dir / "INCAR").write_text(INCAR_ENDPOINT)
    (ep0_dir / "KPOINTS").write_text(KPOINTS)
    (ep0_dir / "submit.sh").write_text(
        submit_endpoint_sh(f"{phase['name'][:6]}_bn_ep0"))
    if potcar.exists():
        shutil.copy(potcar, ep0_dir / "POTCAR")

    # Endpoint 10: hopper at B, vacancy at A  (swap idx_hop and idx_vac)
    ep10 = build_vacancy_structure(p, phase["idx_vac"], phase["idx_hop"])
    ep10_dir = hop_dir / "10_relax"
    ep10_dir.mkdir(exist_ok=True)
    write_poscar(ep10, ep10_dir / "POSCAR",
                 comment=f"{phase['name']} bottleneck ep10")
    (ep10_dir / "INCAR").write_text(INCAR_ENDPOINT)
    (ep10_dir / "KPOINTS").write_text(KPOINTS)
    (ep10_dir / "submit.sh").write_text(
        submit_endpoint_sh(f"{phase['name'][:6]}_bn_ep10"))
    if potcar.exists():
        shutil.copy(potcar, ep10_dir / "POTCAR")

    print(f"  → {hop_dir.relative_to(CALC_DIR)}/{{00,10}}_relax")

print("\nNext: submit all 6 endpoint relaxations, then run setup_neb_bottleneck_images.py")
