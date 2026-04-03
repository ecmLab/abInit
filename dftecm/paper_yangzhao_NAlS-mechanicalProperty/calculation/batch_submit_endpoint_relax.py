"""
Submit endpoint relaxation jobs for NEB refinement.
Run from calculation/:  python batch_submit_endpoint_relax.py
"""

import subprocess
from pathlib import Path

CALC_DIR = Path(__file__).parent
PHASES   = [
    "mp-560538_Na3AlS3",
    "user_Na5AlS4_Na5AlS4",
    "JVASP-12818_Na5InS4",
]

submitted, failed, skipped = [], [], []

for phase in PHASES:
    for img_id in ["00", "06"]:
        job_dir = CALC_DIR / phase / "neb" / "hop1" / f"{img_id}_relax"
        if not (job_dir / "INCAR").exists():
            skipped.append(str(job_dir))
            print(f"SKIP  {phase} ep{img_id}  (run relax_neb_endpoints.py first)")
            continue

        result = subprocess.run(
            ["sbatch", "submit.sh"], cwd=job_dir,
            capture_output=True, text=True
        )
        label = f"{phase}/neb/hop1/{img_id}_relax"
        if result.returncode == 0:
            msg = f"{label:60s}  {result.stdout.strip()}"
            submitted.append(msg)
            print("OK  " + msg)
        else:
            msg = f"{label}: {result.stderr.strip()}"
            failed.append(msg)
            print("ERR " + msg)

print(f"\nSubmitted: {len(submitted)}  |  Failed: {len(failed)}  |  Skipped: {len(skipped)}")
