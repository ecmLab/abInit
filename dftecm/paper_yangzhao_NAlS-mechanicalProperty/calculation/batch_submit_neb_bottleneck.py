"""
Submit bottleneck NEB jobs for all three phases.
Run from calculation/:  python batch_submit_neb_bottleneck.py
"""

import subprocess
from pathlib import Path

CALC_DIR = Path(__file__).parent
PHASES = [
    "mp-560538_Na3AlS3",
    "user_Na5AlS4_Na5AlS4",
    "JVASP-12818_Na5InS4",
]

submitted, failed, skipped = [], [], []

for phase in PHASES:
    job_dir = CALC_DIR / phase / "neb" / "hop_bottleneck" / "neb_images"
    if not (job_dir / "INCAR").exists():
        skipped.append(phase)
        print(f"SKIP {phase}  (run setup_neb_bottleneck_images.py first)")
        continue

    result = subprocess.run(
        ["sbatch", "submit.sh"], cwd=job_dir,
        capture_output=True, text=True
    )
    label = f"{phase}/neb/hop_bottleneck/neb_images"
    if result.returncode == 0:
        msg = f"{label}  {result.stdout.strip()}"
        submitted.append(msg); print("OK  " + msg)
    else:
        msg = f"{label}: {result.stderr.strip()}"
        failed.append(msg); print("ERR " + msg)

print(f"\nSubmitted: {len(submitted)}  Failed: {len(failed)}  Skipped: {len(skipped)}")
