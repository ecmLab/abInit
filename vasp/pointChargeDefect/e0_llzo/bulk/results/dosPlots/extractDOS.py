import os
import shutil

# U values for La and Zr orbitals
LaU = [0.0,4.0,4.5,5.0,5.5,6.5,7.0,7.5,8.0,8.5,9.5,10.0];
ZrU = [0.0,0.5,1.5,2.0,2.5,3.5,4.0];


Prefix=os.getcwd();

for i1 in LaU:
  for i2 in ZrU:
     pth_atom = '../'+'LaU_'+str(i1)+'_ZrU_'+str(i2)+'/3.gamma_band/PDOS_atom.pdf';
     out_atom = 'La'+str(i1)+'Zr'+str(i2)+'_atom.pdf'
     try: shutil.copy(pth_atom,out_atom);
     except: print('Not done'); continue;
    
     pth_obt = '../'+'LaU_'+str(i1)+'_ZrU_'+str(i2)+'/3.gamma_band/PDOS_orbital.pdf';
     out_obt = 'La'+str(i1)+'Zr'+str(i2)+'_orbital.pdf'
     shutil.copy(pth_obt,out_obt);     

