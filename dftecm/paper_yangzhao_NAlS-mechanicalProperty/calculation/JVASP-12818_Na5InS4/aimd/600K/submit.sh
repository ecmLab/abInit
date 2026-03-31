#!/bin/bash -l
#SBATCH --job-name=JVASP-12818_600K
#SBATCH --account=membrane
#SBATCH --partition=tier3

#SBATCH --mail-user=slack:@qhteme

#SBATCH -t 3-00:00:00

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=500M

spack load vasp@6.3.2
srun vasp_std
