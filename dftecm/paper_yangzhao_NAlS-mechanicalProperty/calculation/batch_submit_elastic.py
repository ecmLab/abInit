"""
Submit all elastic jobs.
Run from the calculation/ directory on the cluster AFTER all relax jobs are done.

Step 1 — copy relaxed structures:
  python batch_submit_elastic.py --prepare

Step 2 — submit:
  python batch_submit_elastic.py
"""

import shutil
import subprocess
import sys
from pathlib import Path

CALC_DIR = Path(__file__).parent
PREPARE  = "--prepare" in sys.argv

submitted, failed, skipped = [], [], []

for job_dir in sorted(CALC_DIR.iterdir()):
    if not job_dir.is_dir() or job_dir.name.startswith("."):
        continue
    relax_dir   = job_dir / "relax"
    elastic_dir = job_dir / "elastic"

    if PREPARE:
        contcar = relax_dir / "CONTCAR"
        if contcar.exists() and contcar.stat().st_size > 0:
            shutil.copy(contcar, elastic_dir / "POSCAR")
            print(f"Copied CONTCAR → elastic/POSCAR  ({job_dir.name})")
        else:
            print(f"WARN  no CONTCAR yet: {job_dir.name}")
            skipped.append(job_dir.name)
            continue

    if not (elastic_dir / "POSCAR").exists():
        skipped.append(job_dir.name)
        continue

    result = subprocess.run(
        ["sbatch", "submit.sh"], cwd=elastic_dir, capture_output=True, text=True
    )
    if result.returncode == 0:
        msg = f"{job_dir.name:45s}  {result.stdout.strip()}"
        submitted.append(msg)
        print("OK  " + msg)
    else:
        msg = f"{job_dir.name}: {result.stderr.strip()}"
        failed.append(msg)
        print("ERR " + msg)

print(f"\nSubmitted: {len(submitted)}  |  Failed: {len(failed)}  |  Skipped: {len(skipped)}")
