#!/usr/bin/env python
""" 
Generate vasp input for geometry relaxation using standard settings as defined in pymatgen.
Default setting is MPRelaxSet (for Materials Project compatibility).
Note that these settings can be changed by users.

Author: Valentina Lacivita 

Date: August 2018
"""

import os
import argparse
import numpy as np

from pymatgen.core import Structure
from pymatgen.io.vasp.sets import MPRelaxSet, MITRelaxSet, MVLScanRelaxSet, MPHSERelaxSet

parser = argparse.ArgumentParser()
parser.add_argument("input", type=str, help="path to input structure, e.g. cif file")
parser.add_argument("--output", type=str, help="output path", default="./")
parser.add_argument("--primitive", help="use primitive cell", action="store_true")
parser.add_argument("--mit", help="use MITRelaxSet", action="store_true")
parser.add_argument("--scan", help="use MVLScanRelaxSet", action="store_true")
parser.add_argument("--hse", help="use MPHSERelaxSet", action="store_true")

args = parser.parse_args()
print(args)
####

structure = Structure.from_file(args.input, primitive=args.primitive)

if args.mit:
	relaxset = MITRelaxSet(structure)
elif args.scan:
	relaxset = MVLScanRelaxSet(structure)
elif args.hse:
	relaxset = MPHSERelaxSet(structure)
else:
	relaxset = MPRelaxSet(structure)

relaxset.write_input(args.output)
