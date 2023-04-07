#!/bin/bash -l
#

#SBATCH --mail-user tek5268@rit.edu
#SBATCH --mail-type=ALL

#SBATCH -t 0-12:0:0

#SBATCH -A membrane -p tier3

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem-per-cpu=3g

export ESPRESSO_PSEUDO=/home/tek5268/qe_potentials
printenv ESPRESSO_PSEUDO

srun -n 2 pw.x < co2.relax.in
