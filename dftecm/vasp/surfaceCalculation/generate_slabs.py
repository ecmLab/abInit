# Input bulk structure
from pymatgen.io.vasp.inputs import Poscar
_poscar = Poscar.from_file("llzo/POSCAR_hseRelax")
_poscar.structure

# Generation Slab Structures
from pymatgen.core.surface import SlabGenerator
slab_structs = SlabGenerator(_poscar.structure, [0,1,0],20,20).get_slab()

# Generate VASP input files
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.io.vasp.sets import MPHSERelaxSet
#user_settings = {'potcar_functional': 'PBE_52', 'user_incar_settings': {'EDIFFG': 1e-5, 'LHFCALC':True, 'HFSCREEN':0.2}}
#vasp_input = MPRelaxSet(_poscar.structure, **user_settings)
user_settings = {'potcar_functional': 'PBE_52', 'user_incar_settings': {'EDIFF': 1e-6,'LCHARG':False, 'AEXX':0.3}}
vasp_input = MPHSERelaxSet(slab_structs, **user_settings)
vasp_input.write_input('.')

