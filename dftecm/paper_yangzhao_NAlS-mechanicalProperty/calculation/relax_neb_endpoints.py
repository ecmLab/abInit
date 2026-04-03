"""
Step 1 of NEB refinement: relax the vacancy-containing endpoint structures
(images 00 and 06 from the existing NEB setup) with ionic positions free
but cell shape/volume fixed (ISIF=2).

Creates:
  <phase>/neb/hop1/00_relax/  ← relax starting endpoint
  <phase>/neb/hop1/06_relax/  ← relax ending endpoint

Run from calculation/:  python relax_neb_endpoints.py
After jobs finish, run:  python restart_neb_from_endpoints.py
"""

import shutil
from pathlib import Path

CALC_DIR = Path(__file__).parent

PHASES = [
    "mp-560538_Na3AlS3",
    "user_Na5AlS4_Na5AlS4",
    "JVASP-12818_Na5InS4",
]

INCAR_ENDPOINT = """\
SYSTEM = endpoint_relax

! Electronic structure
ENCUT  = 500
EDIFF  = 1E-05
PREC   = Normal
ALGO   = Fast
LREAL  = Auto
ISYM   = 0

! Ionic relaxation (cell fixed, ions free)
IBRION = 2
ISIF   = 2
NSW    = 300
EDIFFG = -0.02

! Output
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


def submit_sh(job_name):
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


for phase in PHASES:
    neb_hop = CALC_DIR / phase / "neb" / "hop1"
    potcar = neb_hop / "POTCAR"

    for img_id in ["00", "06"]:
        img_poscar = neb_hop / img_id / "POSCAR"
        if not img_poscar.exists():
            print(f"WARN  {phase}  img{img_id}/POSCAR not found — skip")
            continue

        out_dir = neb_hop / f"{img_id}_relax"
        out_dir.mkdir(exist_ok=True)

        shutil.copy(img_poscar, out_dir / "POSCAR")
        (out_dir / "INCAR").write_text(INCAR_ENDPOINT)
        (out_dir / "KPOINTS").write_text(KPOINTS)
        (out_dir / "submit.sh").write_text(
            submit_sh(f"{phase.split('_')[0]}_ep{img_id}")
        )
        if potcar.exists():
            shutil.copy(potcar, out_dir / "POTCAR")

        print(f"  {phase}  img{img_id} → {out_dir.relative_to(CALC_DIR)}")

print("\nDone. Submit with: python batch_submit_endpoint_relax.py")
