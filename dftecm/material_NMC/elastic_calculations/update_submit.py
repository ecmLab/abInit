#!/usr/bin/env python3
"""
Update all submit.sh files based on sid.sh template.
"""

import os
from pathlib import Path

CALC_DIR = Path(__file__).parent

def generate_submit_sh(job_name):
    """Generate submit.sh content based on sid.sh template."""
    return f"""#!/bin/bash -l
# The -l above is required to get the full environment with modules

#SBATCH --job-name={job_name}
##SBATCH --account=purewater
#SBATCH --partition=tier3
#SBATCH --partition=debug

# displays outputs/err
##SBATCH --output=%x_%j.out
##SBATCH --error=%x_%j.err
#SBATCH --mail-user=slack:@qhteme

#  wall-clock time for tier3
#SBATCH -t 0-48:00:00

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem-per-cpu=500M

#spack load vasp@6.3.2 /vpzm2zw
spack load vasp@6.3.2/iriumba
## Commands to run:

# srun -n vasp_std
srun vasp_std
"""


def main():
    # Find all submit.sh files
    submit_files = list(CALC_DIR.glob("*/Li_*/*/submit.sh"))

    print(f"Found {len(submit_files)} submit.sh files to update")

    updated = 0
    for submit_path in sorted(submit_files):
        # Extract job name from path: NMC-811/Li_100/relax -> NMC811_Li100_relax
        parts = submit_path.relative_to(CALC_DIR).parts
        nmc_type = parts[0].replace("-", "")  # NMC-811 -> NMC811
        li_content = parts[1]  # Li_100
        calc_type = parts[2]  # relax or elastic

        job_name = f"{nmc_type}_{li_content}_{calc_type}"

        # Write updated submit.sh
        content = generate_submit_sh(job_name)
        with open(submit_path, 'w') as f:
            f.write(content)

        print(f"  Updated: {submit_path.relative_to(CALC_DIR)} -> {job_name}")
        updated += 1

    print(f"\nUpdated {updated}/{len(submit_files)} submit.sh files")


if __name__ == "__main__":
    main()
