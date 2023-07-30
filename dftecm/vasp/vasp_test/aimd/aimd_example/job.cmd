#!/bin/bash
#SBATCH -J test_tuto 
#SBATCH -N 1
#SBATCH -p debug
#SBATCH -t 00:30:00
#SBATCH --mail-user=vlacivita@lbl.gov
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

echo "Ncores: 96"
echo "Runtime: 30 min"

#run the application:
srun -n 24 -c 2 --cpu_bind=cores /global/homes/l/lacivita/bin/edison/5.4.4/intel/17.0.2.174/mzaaw5t/bin/vasp_gam
