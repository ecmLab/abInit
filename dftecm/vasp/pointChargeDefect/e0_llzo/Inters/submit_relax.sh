
for pthi in */; do
  cd $pthi/
  for pthj in */; do
    cd $pthj/

     mkdir -p 1.gamma_relax
     cp INCAR 1.gamma_relax/INCAR
     cp POSCAR 1.gamma_relax/POSCAR
     cp KPOINTS 1.gamma_relax/KPOINTS
     cp POTCAR 1.gamma_relax/POTCAR
     cp ../../qsub_gamma.sh 1.gamma_relax/
     cd 1.gamma_relax
      printf "Automatic\n0\nGamma\n1 1 1" > KPOINTS
      sbatch qsub_gamma.sh
     cd ..

    cd ..
   done
   cd ..
done
   
