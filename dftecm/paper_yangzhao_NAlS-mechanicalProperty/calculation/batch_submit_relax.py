"""
Submit all relax jobs.
Run from the calculation/ directory on the cluster:
  python batch_submit_relax.py
"""

import subprocess
from pathlib import Path

CALC_DIR = Path(__file__).parent

submitted, failed, skipped = [], [], []

for job_dir in sorted(CALC_DIR.iterdir()):
    if not job_dir.is_dir() or job_dir.name.startswith("."):
        continue
    relax_dir = job_dir / "relax"
    if not (relax_dir / "POSCAR").exists():
        skipped.append(job_dir.name)
        continue

    result = subprocess.run(
        ["sbatch", "submit.sh"], cwd=relax_dir, capture_output=True, text=True
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
