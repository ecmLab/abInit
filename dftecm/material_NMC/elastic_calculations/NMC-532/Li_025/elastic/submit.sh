#!/bin/bash -l
# The -l above is required to get the full environment with modules

#SBATCH --job-name=NMC532_Li_025_elastic
#SBATCH --account=purewater
#SBATCH --partition=tier3
##SBATCH --partition=debug

# displays outputs/err
##SBATCH --output=%x_%j.out
##SBATCH --error=%x_%j.err
#SBATCH --mail-user=slack:@qhteme

#  wall-clock time for tier3
#SBATCH -t 0-48:00:00

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem-per-cpu=500M

spack load vasp@6.3.2
## Commands to run:

# srun -n vasp_std
srun vasp_std
