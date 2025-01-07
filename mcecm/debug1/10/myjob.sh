#!/bin/bash
# Job name:
#SBATCH --job-name=MC
#
# Partition:
#SBATCH --partition=tier3
#
# Processors:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#
# Wall clock limit:
#SBATCH --time=20:00:00
#SBATCH --mem=50g


## Commands to run:
#spack load lammps@20230208 /cuxhkce
#mpirun --mca btl ^sm -n 4 python mc_mpi.py
mpirun -n 32 python sporc_mpi6.py
