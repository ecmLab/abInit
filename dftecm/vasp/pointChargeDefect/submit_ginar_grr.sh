#!/bin/bash -l
#$ -N opt
#$ -o vasp.out
#$ -pe impi 16
#$ -V
#$ -cwd

printf "STARTED\n" > STATUS

mkdir 1.gamma_relax
cp INCAR 1.gamma_relax/INCAR
cp POSCAR 1.gamma_relax/POSCAR
cp POTCAR 1.gamma_relax/POTCAR
cd 1.gamma_relax
printf "Automatic\n0\nGamma\n1 1 1" > KPOINTS
sed -i 's/LCHARG = True/LCHARG = False/' INCAR
sed -i 's/EDIFFG = -0.02//' INCAR
mpiexec.hydra pvasp.5.4.4.intel.gamma
grep "reached required accuracy" OUTCAR
if [ $? -ne 0 ] ; then printf "FAILED TO GAMMA_RELAX\n" >> ../STATUS ; exit ; fi
cd ..
printf "GAMMA_RELAXED\n" >> STATUS

mkdir 2.full_relax
cp INCAR 2.full_relax/INCAR
cp 1.gamma_relax/CONTCAR 2.full_relax/POSCAR
cp KPOINTS 2.full_relax/KPOINTS
cp POTCAR 2.full_relax/POTCAR
cd 2.full_relax
sed -i 's/EDIFFG = -0.02//' INCAR
sed -i 's/LCHARG = True/LCHARG = False/' INCAR
mpiexec.hydra pvasp.5.4.4.intel
grep "reached required accuracy" OUTCAR
if [ $? -ne 0 ] ; then printf "FAILED TO FULL_RELAX\n" >> ../STATUS ; exit ; fi
cd ..
printf "FULL_RELAXED\n" >> STATUS

mkdir 3.double_relax
cp INCAR 3.double_relax/INCAR
cp 2.full_relax/CONTCAR 3.double_relax/POSCAR
cp KPOINTS 3.double_relax/KPOINTS
cp POTCAR 3.double_relax/POTCAR
cd 3.double_relax
sed -i 's/IBRION = 2/IBRION = 1/' INCAR
mpiexec.hydra pvasp.5.4.4.intel
grep "reached required accuracy" OUTCAR
if [ $? -ne 0 ] ; then printf "FAILED TO DOUBLE_RELAX\n" >> ../STATUS ; exit ; fi
cd ..
printf "DOUBLE_RELAXED\n" >> STATUS

printf "DONE\n" >> STATUS
