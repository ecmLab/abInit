"""
Step 2: after hop_bottleneck endpoint relaxations finish, generate NEB images
using IDPP interpolation (avoids atomic clashes from linear interpolation).

IDPP = Image Dependent Pair Potential — optimises the initial path so that
interatomic distances change smoothly, preventing atoms from passing through
each other.  This is critical for longer hops (≥ 3.5 Å).

Run from calculation/:  python setup_neb_bottleneck_images.py
Then:  python batch_submit_neb_bottleneck.py
"""

import shutil
import numpy as np
from pathlib import Path

import ase.io
from ase.mep import NEB

CALC_DIR = Path(__file__).parent
N_IMAGES = 9

PHASES = [
    {"dir": "mp-560538_Na3AlS3",    "name": "Na3AlS3",  "idx_hop": 14, "idx_vac": 18},
    {"dir": "user_Na5AlS4_Na5AlS4", "name": "Na5AlS4",  "idx_hop": 10, "idx_vac": 12},
    {"dir": "JVASP-12818_Na5InS4",  "name": "Na5InS4",  "idx_hop": 1,  "idx_vac": 8},
]

NEB_INCAR = """\
SYSTEM = NEB_bottleneck

ENCUT  = 500
EDIFF  = 1E-06
PREC   = Accurate
ALGO   = Fast
LREAL  = Auto
ISYM   = 0

IBRION = 3
POTIM  = 0
NSW    = 1000
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


def submit_neb_sh(job_name, n_images):
    cores_per_node = 32
    cores_per_image = cores_per_node // n_images
    ntasks = cores_per_image * n_images
    return f"""\
#!/bin/bash -l
#SBATCH --job-name={job_name}_bn_neb
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


for phase in PHASES:
    phase_dir = CALC_DIR / phase["dir"]
    hop_dir   = phase_dir / "neb" / "hop_bottleneck"
    ep0_con   = hop_dir / "00_relax" / "CONTCAR"
    ep10_con  = hop_dir / "10_relax" / "CONTCAR"
    potcar    = phase_dir / "relax" / "POTCAR"

    missing = [p for p in [ep0_con, ep10_con] if not p.exists()]
    if missing:
        print(f"SKIP {phase['name']}: missing {[m.name for m in missing]}")
        continue

    # Read endpoints with ASE
    start = ase.io.read(str(ep0_con),  format="vasp")
    end   = ase.io.read(str(ep10_con), format="vasp")

    # Fix atom ordering in end so the hopping atom is at the same index as in start.
    # In start (idx_vac removed): hopper sits at local index idx_hop.
    # In end   (idx_hop removed): hopper (formerly idx_vac) sits at local index idx_vac-1.
    # IDPP requires atom i in start ↔ atom i in end; without this fix it interpolates
    # the wrong atoms and the hop never occurs.
    assert phase["idx_hop"] < phase["idx_vac"], "assume idx_hop < idx_vac"
    hop_in_end   = phase["idx_vac"] - 1   # current local index of hopper in end
    hop_in_start = phase["idx_hop"]        # target local index (must match start)
    n = len(end)
    src = list(range(n))
    src.pop(hop_in_end)
    src.insert(hop_in_start, hop_in_end)
    end = end[src]

    # Build NEB image list: endpoint copies + N_IMAGES intermediates
    images = [start.copy() for _ in range(N_IMAGES + 2)]
    images[-1] = end.copy()

    # IDPP interpolation — avoids atomic clashes
    neb = NEB(images, method="improvedtangent")
    neb.interpolate(method="idpp")
    print(f"{phase['name']}: IDPP interpolation done")

    # Write image POSCARs
    neb_dir = hop_dir / "neb_images"
    neb_dir.mkdir(exist_ok=True)

    for k, img in enumerate(images):
        img_dir = neb_dir / f"{k:02d}"
        img_dir.mkdir(exist_ok=True)
        ase.io.write(str(img_dir / "POSCAR"), img, format="vasp", direct=True)

    # Write VASP inputs
    (neb_dir / "INCAR").write_text(NEB_INCAR.format(n_images=N_IMAGES))
    (neb_dir / "KPOINTS").write_text(KPOINTS)
    (neb_dir / "submit.sh").write_text(
        submit_neb_sh(phase["name"][:6], N_IMAGES))
    if potcar.exists():
        shutil.copy(potcar, neb_dir / "POTCAR")

    # Check min interatomic distance in saddle image to verify no clashes
    saddle = images[N_IMAGES // 2 + 1]
    pos    = saddle.get_positions()
    cell   = saddle.get_cell()
    min_d  = np.inf
    for i in range(len(pos)):
        for j in range(i+1, len(pos)):
            d = pos[j] - pos[i]
            # Apply MIC
            d -= cell.T @ np.round(np.linalg.solve(cell.T, d))
            dist = np.linalg.norm(d)
            if dist < min_d:
                min_d = dist
    print(f"  Min interatomic distance at saddle image: {min_d:.3f} Å "
          f"({'OK' if min_d > 1.5 else 'WARNING: possible clash'})")
    print(f"  → {neb_dir.relative_to(CALC_DIR)}  ({N_IMAGES+2} images)")

print("\nSubmit with: python batch_submit_neb_bottleneck.py")
