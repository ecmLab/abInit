#!/bin/bash -l
#

prefix="pfoa"
jobfile="slurm_rlx.sh"

jobname=$prefix-rlx
outfile=$prefix-rlx.o
errfile=$prefix-rlx.e

spack load quantum-espresso/7hk6bwf

sbatch -J $jobname -o $outfile -e $errfile $jobfile
