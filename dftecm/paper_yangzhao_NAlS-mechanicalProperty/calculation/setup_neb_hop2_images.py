"""
Step 2: after hop2 endpoints converge, interpolate images and write NEB inputs.
Run from calculation/:  python setup_neb_hop2_images.py
"""

import shutil
import numpy as np
from pathlib import Path

CALC_DIR  = Path(__file__).parent
PHASE_DIR = CALC_DIR / "mp-560538_Na3AlS3"
N_IMAGES  = 9

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


def submit_neb_sh(n_images):
    cores_per_node = 32
    cores_per_image = cores_per_node // n_images
    ntasks = cores_per_image * n_images
    return f"""\
#!/bin/bash -l
#SBATCH --job-name=Na3AlS3_neb_hop2
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


hop_dir = PHASE_DIR / "neb" / "hop2"
ep0_contcar  = hop_dir / "00_relax" / "CONTCAR"
ep10_contcar = hop_dir / "10_relax" / "CONTCAR"
potcar       = PHASE_DIR / "relax" / "POTCAR"

missing = [p for p in [ep0_contcar, ep10_contcar] if not p.exists()]
if missing:
    print(f"ERROR: missing {[str(m) for m in missing]}")
    print("Run endpoint relaxations first (submit.sh in 00_relax and 10_relax).")
    raise SystemExit(1)

p0  = read_poscar(ep0_contcar)
p10 = read_poscar(ep10_contcar)

neb_dir = hop_dir / "hop1_v2"
neb_dir.mkdir(exist_ok=True)

# Write interpolated images
for k in range(N_IMAGES + 2):
    t = k / (N_IMAGES + 1)
    diff = p10["positions"] - p0["positions"]
    diff -= np.round(diff)
    pos = (p0["positions"] + t * diff) % 1.0
    img = dict(lattice=p0["lattice"], elements=p0["elements"],
               counts=p0["counts"], positions=pos)
    img_dir = neb_dir / f"{k:02d}"
    img_dir.mkdir(exist_ok=True)
    write_poscar(img, img_dir / "POSCAR",
                 comment=f"Na3AlS3 hop2 image {k:02d}")

(neb_dir / "INCAR").write_text(NEB_INCAR.format(n_images=N_IMAGES))
(neb_dir / "KPOINTS").write_text(KPOINTS)
(neb_dir / "submit.sh").write_text(submit_neb_sh(N_IMAGES))
if potcar.exists():
    shutil.copy(potcar, neb_dir / "POTCAR")

print(f"Written {N_IMAGES + 2} images to {neb_dir.relative_to(CALC_DIR)}")
print(f"Submit: cd {neb_dir.relative_to(CALC_DIR)} && sbatch submit.sh")
