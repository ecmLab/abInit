#!/bin/bash
#SBATCH --job-name=charge_defect
#SBATCH --qos=regular
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=32
#SBATCH --constraint=knl
#SBATCH --mail-user=howardtu@lbl.gov
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes

#export PATH=$PATH:/global/homes/h/howietu/VASP_bin/Cori/KNL
vasp_std=/global/homes/h/howietu/VASP_bin/Cori/KNL/vasp_std
vasp_gam=/global/homes/h/howietu/VASP_bin/Cori/KNL/vasp_gam

#OpenMP settings:
export OMP_NUM_THREADS=4
export OMP_PLACES=threads
export OMP_PROC_BIND=true

## Commands to run:
srun -n 32 -c 2 --cpu_bind=cores $vasp_gam   # For cori: 32 cores per node
