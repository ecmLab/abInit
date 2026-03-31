#!/bin/bash -l
#SBATCH --job-name=JVASP-12818_neb
#SBATCH --account=membrane
#SBATCH --partition=tier3

#SBATCH --mail-user=slack:@qhteme

#SBATCH -t 10:00:00

#SBATCH --ntasks=30
#SBATCH --ntasks-per-node=30
#SBATCH --mem-per-cpu=500M

spack load vasp@6.3.2
srun vasp_std
