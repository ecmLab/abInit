"""
Set up hop2 NEB for Na3AlS3: the shortest INTER-cluster Na hop (3.313 Å).

Background
----------
The hop1 NEB (atoms 0–17, 3.225 Å) found Ea ≈ 0.010 eV because that hop
is *intra-cluster* — atoms 0 and 17 are within the same Na triplet
{0, 7, 17} and going nowhere useful.  Long-range Na diffusion requires
inter-cluster hops; the shortest one is atoms 11–16 at 3.313 Å.

Output
------
  mp-560538_Na3AlS3/neb/hop2/{00_relax, hop1_v2/00..10, 06_relax}/
  (endpoint relaxation dirs + interpolated images + INCAR/KPOINTS/submit)

Run from calculation/:  python setup_neb_hop2_Na3AlS3.py
Then:  python batch_submit_hop2.py
"""

import shutil
import numpy as np
from pathlib import Path

CALC_DIR  = Path(__file__).parent
PHASE_DIR = CALC_DIR / "mp-560538_Na3AlS3"
N_IMAGES  = 9

# Verified shortest inter-cluster hop
IDX_HOP = 11   # hopping Na (stays in structure)
IDX_VAC = 16   # Na that becomes vacancy (destination)
HOP_DIST = 3.3127   # Å

INCAR_ENDPOINT = """\
SYSTEM = endpoint_relax_hop2

ENCUT  = 500
EDIFF  = 1E-05
PREC   = Normal
ALGO   = Fast
LREAL  = Auto
ISYM   = 0

IBRION = 2
ISIF   = 2
NSW    = 300
EDIFFG = -0.02

LWAVE  = .TRUE.
LCHARG = .FALSE.
"""

NEB_INCAR = """\
SYSTEM = NEB_hop2

ENCUT  = 500
EDIFF  = 1E-05
PREC   = Normal
ALGO   = Fast
LREAL  = Auto
ISYM   = 0

IBRION = 3
POTIM  = 0
NSW    = 500
EDIFFG = -0.05
IMAGES = {n_images}
SPRING = -5
LCLIMB = .TRUE.

ISTART = 0
LWAVE  = .FALSE.
LCHARG = .FALSE.

NPAR   = {n_images}
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


def submit_neb_sh(job_name, n_images):
    cores_per_node = 32
    cores_per_image = cores_per_node // n_images
    ntasks = cores_per_image * n_images
    return f"""\
#!/bin/bash -l
#SBATCH --job-name={job_name}
#SBATCH --account=membrane
#SBATCH --partition=tier3
#SBATCH --mail-user=slack:@qhteme
#SBATCH -t 10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks={ntasks}
#SBATCH --ntasks-per-node={ntasks}
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
    """Remove idx_vac from the structure (creates vacancy at destination)."""
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
        lattice=p["lattice"],
        elements=p["elements"],
        counts=new_counts,
        positions=p["positions"][base_idx].copy(),
    )


def build_images(p0, p1, n_images):
    """PBC-aware linear interpolation between two endpoint structures."""
    images = []
    for k in range(n_images + 2):
        t = k / (n_images + 1)
        diff = p1["positions"] - p0["positions"]
        diff -= np.round(diff)
        pos = (p0["positions"] + t * diff) % 1.0
        images.append(dict(
            lattice=p0["lattice"],
            elements=p0["elements"],
            counts=p0["counts"],
            positions=pos,
        ))
    return images


# ── Read relaxed unit cell ─────────────────────────────────────────────────
contcar = PHASE_DIR / "relax" / "CONTCAR"
potcar  = PHASE_DIR / "relax" / "POTCAR"
p = read_poscar(contcar)

# Verify hop distance
pos_hop = p["positions"][IDX_HOP]
pos_vac = p["positions"][IDX_VAC]
d = pos_vac - pos_hop
d -= np.round(d)
dist = np.linalg.norm(d @ p["lattice"])
print(f"Hop: atom {IDX_HOP} → atom {IDX_VAC}  ({dist:.4f} Å)  [expected {HOP_DIST} Å]")

# ── Endpoint structures ────────────────────────────────────────────────────
hop_dir = PHASE_DIR / "neb" / "hop2"
hop_dir.mkdir(parents=True, exist_ok=True)

# Endpoint 00: hopping Na at site A, vacancy at site B
ep0_struct = build_vacancy_structure(p, IDX_HOP, IDX_VAC)
ep0_dir = hop_dir / "00_relax"
ep0_dir.mkdir(exist_ok=True)
write_poscar(ep0_struct, ep0_dir / "POSCAR",
             comment="Na3AlS3 hop2 endpoint 00 (Na at A, vac at B)")
(ep0_dir / "INCAR").write_text(INCAR_ENDPOINT)
(ep0_dir / "KPOINTS").write_text(KPOINTS)
(ep0_dir / "submit.sh").write_text(submit_endpoint_sh("Na3AlS3_h2_ep00"))
if potcar.exists():
    shutil.copy(potcar, ep0_dir / "POTCAR")
print(f"  Endpoint 00 → {ep0_dir.relative_to(CALC_DIR)}")

# Endpoint 10: hopping Na at site B, vacancy at site A
# Swap: IDX_VAC stays, IDX_HOP is removed, IDX_VAC is the hopper
ep10_struct = build_vacancy_structure(p, IDX_VAC, IDX_HOP)
ep10_dir = hop_dir / "10_relax"
ep10_dir.mkdir(exist_ok=True)
write_poscar(ep10_struct, ep10_dir / "POSCAR",
             comment="Na3AlS3 hop2 endpoint 10 (Na at B, vac at A)")
(ep10_dir / "INCAR").write_text(INCAR_ENDPOINT)
(ep10_dir / "KPOINTS").write_text(KPOINTS)
(ep10_dir / "submit.sh").write_text(submit_endpoint_sh("Na3AlS3_h2_ep10"))
if potcar.exists():
    shutil.copy(potcar, ep10_dir / "POTCAR")
print(f"  Endpoint 10 → {ep10_dir.relative_to(CALC_DIR)}")

print(f"\nStep 1: submit endpoint relaxations")
print(f"  cd {ep0_dir.relative_to(CALC_DIR)} && sbatch submit.sh")
print(f"  cd {ep10_dir.relative_to(CALC_DIR)} && sbatch submit.sh")
print(f"\nStep 2: after endpoints converge, run:")
print(f"  python setup_neb_hop2_images.py")
