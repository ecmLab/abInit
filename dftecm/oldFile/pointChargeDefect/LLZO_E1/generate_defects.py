# Input bulk structure
from pymatgen.io.vasp.inputs import Poscar
_poscar = Poscar.from_file("bulk/Li7La3Zr2O12_0_POSCAR")
_poscar.structure

# Generation Defect Structures
from pycdt.core.defectsmaker import ChargedDefectsStructures
charged_defects_structs = ChargedDefectsStructures(_poscar.structure)
from pycdt.utils.vasp import DefectRelaxSet
defect_struct = charged_defects_structs.get_ith_supercell_of_defect_type(0, 'vacancies')

# Generate VASP input files
import os
for key in charged_defects_structs.defects.keys():
    if key == "bulk":
        continue
    os.makedirs(f"{key}", exist_ok=True)
    for vac_dict in charged_defects_structs.defects[key]:
#         print(vac_dict, type(vac_dict))
        os.makedirs(f'{key}/{vac_dict["name"]}', exist_ok=True)
        user_settings = {'potcar_functional': 'PBE_52', 'user_incar_settings': {'EDIFFG': -1e-3, 'EDIFF': 1e-8, 'LASPH': True}}
        vasp_input = DefectRelaxSet(vac_dict['supercell']['structure'], **user_settings) # defect_struct is defect supercell
        vasp_input.write_input(f'{key}/{vac_dict["name"]}')
