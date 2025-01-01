#!/bin/bash
# Job name:
#SBATCH --job-name=MC
#
# Partition:
#SBATCH --partition=tier3
#
# Processors:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#
# Wall clock limit:
#SBATCH --time=01:00:00
#SBATCH --mem=10g


## Commands to run:
#spack load lammps@20230208 /cuxhkce
#srun lmp -log none -in mc_run.in
#srun python mc_mpi.py
#mpirun --mca btl ^sm -n 4 python mc_mpi.py
mpirun -n 16 python mc_mpi.py
