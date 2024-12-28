#!/bin/bash
# Job name:
#SBATCH --job-name=NMC_LPS
#
# Partition:
#SBATCH --partition=tier3
#
# Processors:
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=8
#
# Wall clock limit:
#SBATCH --time=01:00:00
#SBATCH --mem=10g


## Commands to run:
spack load lammps@20230208 /cuxhkce
srun lmp -log none -in mc_run.in
