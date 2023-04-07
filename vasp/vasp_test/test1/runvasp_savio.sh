#!/bin/bash
# Job name:
#SBATCH --job-name=test-vasp
#SBATCH --account=fc_ceder
# Partition:
#SBATCH --partition=savio
#
# Processors:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#
# Wall clock limit:
#SBATCH --time=00:30:00

ulimit -s unlimited
module load intel/2016.4.072
module load mkl/2016.4.072
module load openmpi/2.0.2-intel
module load fftw/3.3.7
mpirun /path/to/your/vasp/binaries