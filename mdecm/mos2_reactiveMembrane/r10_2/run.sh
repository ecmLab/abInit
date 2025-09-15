#!/bin/bash -l
# The -l above is required to get the full environment with modules

#SBATCH --job-name=MoS2
#SBATCH --account=purewater
#SBATCH --partition=tier3
##SBATCH --partition=debug

# displays outputs/err
##SBATCH --output=%x_%j.out
##SBATCH --error=%x_%j.err
#SBATCH --mail-user=slack:@qhteme

#  wall-clock time for tier3
#SBATCH -t 0-12:00:00
##SBATCH -t 0-02:00:00

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem-per-cpu=500M

#load the lammps module
spack load lammps@20230802.2

srun lmp -in lmp.in
