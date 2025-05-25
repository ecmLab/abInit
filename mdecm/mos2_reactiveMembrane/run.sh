#!/bin/bash -l
# The -l above is required to get the full environment with modules

#SBATCH --job-name=MoS2
#SBATCH --account=purewater
##SBATCH --partition=tier3
#SBATCH --partition=debug

# displays outputs/err
#SBATCH --output=%x_%j.out
#SBATCH --mail-user=slack:@qhteme
#SBATCH --error=%x_%j.err

#  wall-clock time for tier3
##SBATCH -t 0-60:00:00
#SBATCH -t 0-02:00:00

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem-per-cpu=20g

#load the lammps module
spack load lammps@20230802.2

srun lmp -in lmp.in
