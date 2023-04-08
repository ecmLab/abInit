#!/bin/bash -l
#$ -N opt
#$ -o defect.out
#$ -pe impi 16
#$ -V
#$ -cwd

printf "STARTED\n" > STATUS

mpiexec.hydra pvasp.5.4.4.intel.gamma
grep "reached required accuracy" OUTCAR
if [ $? -ne 0 ] ; then printf "FAILED TO gamma_band\n" >> ../STATUS ; exit ; fi
cd ..
printf "gamma_band\n" >> STATUS

printf "DONE\n" >> STATUS
