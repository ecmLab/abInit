
for pthi in */; do
  cd $pthi/
  for pthj in */; do
    cd $pthj/

     cd 1.gamma_relax
     mv CONTCAR POSCAR
     cp ../../../qsub_gamma.sh ./
     sbatch qsub_gamma.sh
     cd ..

    cd ..
   done
   cd ..
done
   
