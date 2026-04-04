#!/bin/bash -l
#SBATCH --job-name=Na5AlS_bn_neb
#SBATCH --account=membrane
#SBATCH --partition=tier3
#SBATCH --mail-user=slack:@qhteme
#SBATCH -t 10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=27
#SBATCH --ntasks-per-node=27
#SBATCH --mem-per-cpu=2000M

spack load vasp@6.3.2
srun vasp_std
