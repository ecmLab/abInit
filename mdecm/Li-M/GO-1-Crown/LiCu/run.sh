#!/bin/bash -l
# The -l above is required to get the full environment with modules



#SBATCH --job-name=DesalinationModeling
#SBATCH --account=waterdesal
#SBATCH --partition=tier3

# The name of the script is run.sh
#SBATCH -J run.sh

# displays outputs/err
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

#SBATCH --mail-user=m.shahbabaei@gmail.com
#SBATCH --mail-user=slack:@maseme
#SBATCH --mail-type=ALL

#  wall-clock time for tier3
#SBATCH -t 0-60:00:00

#SBATCH --nodes=10
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=20g

#load the lammps module

spack load lammps@20220623.2 /vyejl4g
#module load intelmpi python/3.8.6 cuda/11.1.1
#module load intelmpi/2021.3.0-intel2021.3.0
#module load cuda/11.1.1
#module load python/3.8.6
#module load LAMMPS/3Mar20 
#module load gcc/10.2.0 openmpi/4.0.5-gcc10.2.0 hdf5/1.10.7-gcc10.2.0 LAMMPS/3Mar20

#echo "SLURM_NTASKS: " $SLURM_NTASKS

#ulimit -n 2048
#export OMP_NUM_THREADS=1


srun lmp -in in.meam

