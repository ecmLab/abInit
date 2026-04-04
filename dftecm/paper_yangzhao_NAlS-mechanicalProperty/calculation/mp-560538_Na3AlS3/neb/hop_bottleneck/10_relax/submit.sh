#!/bin/bash -l
#SBATCH --job-name=Na3AlS3_bn_10
#SBATCH --account=membrane
#SBATCH --partition=tier3
#SBATCH --mail-user=slack:@qhteme
#SBATCH -t 5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=2000M

spack load vasp@6.3.2
srun vasp_std
