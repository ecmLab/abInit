#from pymatgen.ext.matproj import MPRester
from mp_api.client import MPRester
from pymatgen.io.vasp.sets import MPScanRelaxSet

# Initialize the REST API interface. You may need to put your own API key in as an arg.
with MPRester("H9GD9SySZmkHxPyxI3T5ylATO9Z8XFAd") as mpr:

  # Chemical formula
   a = mpr.get_structure_by_material_id('mp-77')
#   b = MPScanRelaxSet(a,user_potcar_functional='PBE_64')
   b = MPScanRelaxSet(a)
   b.write_input('.')
