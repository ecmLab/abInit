"""
Submit all AIMD jobs for the three diffusivity phases.
Run from calculation/:  python batch_submit_aimd.py
"""

import subprocess
from pathlib import Path

CALC_DIR = Path(__file__).parent
PHASES   = [
    "mp-560538_Na3AlS3",
    "user_Na5AlS4_Na5AlS4",
    "JVASP-12818_Na5InS4",
]
TEMPS = [600, 800, 1000]

submitted, failed, skipped = [], [], []

for phase in PHASES:
    for T in TEMPS:
        job_dir = CALC_DIR / phase / "aimd" / f"{T}K"
        if not (job_dir / "POSCAR").exists():
            skipped.append(str(job_dir))
            continue

        result = subprocess.run(
            ["sbatch", "submit.sh"], cwd=job_dir,
            capture_output=True, text=True
        )
        label = f"{phase}/{T}K"
        if result.returncode == 0:
            msg = f"{label:55s}  {result.stdout.strip()}"
            submitted.append(msg)
            print("OK  " + msg)
        else:
            msg = f"{label}: {result.stderr.strip()}"
            failed.append(msg)
            print("ERR " + msg)

print(f"\nSubmitted: {len(submitted)}  |  Failed: {len(failed)}  |  Skipped: {len(skipped)}")
