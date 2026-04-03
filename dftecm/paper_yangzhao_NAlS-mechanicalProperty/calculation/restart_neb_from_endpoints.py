"""
Step 2 of NEB refinement: after endpoint relaxations finish, re-interpolate
images and set up a fresh CI-NEB run.

Reads:  <phase>/neb/hop1/00_relax/CONTCAR   (relaxed start)
        <phase>/neb/hop1/06_relax/CONTCAR   (relaxed end)
Writes: <phase>/neb/hop1_v2/{00..06}/POSCAR
        <phase>/neb/hop1_v2/INCAR, KPOINTS, POTCAR, submit.sh

Run from calculation/:  python restart_neb_from_endpoints.py
"""

import shutil
import numpy as np
from pathlib import Path

CALC_DIR = Path(__file__).parent
N_IMAGES = 9   # intermediate images (finer path than original hop1)

PHASES = [
    "mp-560538_Na3AlS3",
    "user_Na5AlS4_Na5AlS4",
    "JVASP-12818_Na5InS4",
]

NEB_INCAR = """\
SYSTEM = NEB_refined

! Electronic structure
ENCUT  = 400
EDIFF  = 1E-05
PREC   = Normal
ALGO   = Fast
LREAL  = Auto
ISYM   = 0

! NEB settings
IBRION = 3
POTIM  = 0
NSW    = 500
EDIFFG = -0.05
IMAGES = {n_images}
SPRING = -5
LCLIMB = .TRUE.

! Read wavefunctions from endpoint relaxations for faster convergence
ISTART = 0
LWAVE  = .FALSE.
LCHARG = .FALSE.

! Parallelisation
NPAR   = {n_images}
"""

KPOINTS = """\
Gamma-point only
 0
Gamma
 1 1 1
 0 0 0
"""


def submit_sh(job_name, n_images):
    cores_per_node = 32
    cores_per_image = cores_per_node // n_images
    ntasks = cores_per_image * n_images
    return f"""\
#!/bin/bash -l
#SBATCH --job-name={job_name}_neb2
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


def interpolate_images(p0, p1, n_images):
    """Linear interpolation with PBC-aware shortest path, n_images+2 total."""
    images = []
    for k in range(n_images + 2):
        t = k / (n_images + 1)
        pos = p0["positions"].copy()
        diff = p1["positions"] - p0["positions"]
        diff -= np.round(diff)   # shortest PBC path per atom
        pos = (p0["positions"] + t * diff) % 1.0
        images.append(dict(
            lattice=p0["lattice"],
            elements=p0["elements"],
            counts=p0["counts"],
            positions=pos,
        ))
    return images


for phase in PHASES:
    neb_hop    = CALC_DIR / phase / "neb" / "hop1"
    ep0_contcar = neb_hop / "00_relax" / "CONTCAR"
    ep6_contcar = neb_hop / "06_relax" / "CONTCAR"
    potcar      = neb_hop / "POTCAR"

    missing = [p for p in [ep0_contcar, ep6_contcar] if not p.exists()]
    if missing:
        print(f"WARN  {phase}: missing {[str(m.name) for m in missing]} — skip")
        print(f"       (run endpoint relaxations first)")
        continue

    p0 = read_poscar(ep0_contcar)
    p1 = read_poscar(ep6_contcar)

    # Verify atom counts match
    if p0["counts"] != p1["counts"]:
        print(f"ERR  {phase}: endpoint atom counts differ — skip")
        continue

    images = interpolate_images(p0, p1, N_IMAGES)

    out_dir = neb_hop.parent / "hop1_v2"
    out_dir.mkdir(exist_ok=True)

    for k, img in enumerate(images):
        img_dir = out_dir / f"{k:02d}"
        img_dir.mkdir(exist_ok=True)
        write_poscar(img, img_dir / "POSCAR",
                     comment=f"{phase} NEB_v2 image {k:02d}")

    (out_dir / "INCAR").write_text(NEB_INCAR.format(n_images=N_IMAGES))
    (out_dir / "KPOINTS").write_text(KPOINTS)
    (out_dir / "submit.sh").write_text(
        submit_sh(phase.split("_")[0], N_IMAGES)
    )
    if potcar.exists():
        shutil.copy(potcar, out_dir / "POTCAR")

    print(f"{phase}")
    print(f"  → {out_dir.relative_to(CALC_DIR)}  ({N_IMAGES+2} images: 00…{N_IMAGES+1:02d})")

print("\nDone. Submit with: python batch_submit_neb_v2.py")
