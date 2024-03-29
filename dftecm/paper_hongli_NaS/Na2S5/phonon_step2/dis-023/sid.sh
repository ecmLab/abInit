#!/bin/bash
# Job name:
#SBATCH --job-name=test-vasp
# Partition:
#SBATCH --partition=tier3 #debug

#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=slack:@maseme
#SBATCH --mail-type=ALL

# Processors:
##SBATCH --nodes=5
##SBATCH --ntasks-per-node=5
##SBATCH --mem-per-cpu=4g

# Sid's changes single node
##SBATCH --ntasks=5
##SBATCH --cpus-per-task=16
##SBATCH --mem-per-cpu=8g

# Sid's changes for MPI
# --nodes should be equal to the total number of cores you want
# divided by --ntasks-per-node if you want to specify cores.
# But I would recommend you let the scheduler do the work for you

#SBATCH --ntasks=32
#SBATCH --mem-per-cpu=8g

# Wall clock limit:
#SBATCH --time=0-12:00:00

#spack load vasp@6.3.2 /vpzm2zw
#spack load vasp /4qy7axx
spack load vasp@6.3.2/iriumba
## Commands to run:
##ulimit -s unlimited
##export I_MPI_ADJUST_REDUCE=3
##export MPIR_CVAR_COLL_ALIAS_CHECK=0

# srun -n vasp_std
srun vasp_std
