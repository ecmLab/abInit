# Input bulk structure
from pymatgen.io.vasp.inputs import Poscar
_poscar = Poscar.from_file("POSCAR")
_poscar.structure


# Generate VASP input files
import os

Prefix=os.getcwd();

FullDir=os.path.join(Prefix,'3.gamma_band');
try: os.chdir(FullDir);
except: print('{} is not done'.format(FullDir)); 
###Extract atomic projection of DOS
try:
   os.system('DOSProcess.py --aqn 4 --atom --elim -1.8 7.8 --figName PDOS_atom.pdf');
except: print('{} is not done'.format(FullDir)); 
###Extract orbital projection of DOS
try:
   os.system('DOSProcess.py --aqn 4 --orb --elim -1.8 7.8 --figName PDOS_orbital.pdf');
except: print('{} is not done'.format(FullDir));
