#!/bin/bash
#SBATCH --job-name=NMC-411-mono_Li100_relax
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=24:00:00
#SBATCH --partition=normal

module load vasp/6.3.0

cd $SLURM_SUBMIT_DIR

mpirun -np $SLURM_NTASKS vasp_std > vasp.out 2>&1

echo "Job completed at $(date)"
