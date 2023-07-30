#!/usr/bin/env python
""" 
Generate vasp input for MD using standard settings as defined in pymatgen.
Default setting is MITMDSet (NVT ensemble). 
Note that these settings can be changed by users.

Author: Valentina Lacivita 

Date: August 2018
"""

import os
import argparse
import numpy as np

from pymatgen.core import Structure
from pymatgen.io.vasp.sets import MITMDSet, MVLNPTMDSet

parser = argparse.ArgumentParser()
parser.add_argument("input", type=str, help="path to input structure, e.g. cif or POSCAR file")
parser.add_argument("--output", type=str, help="output path", default="./")
parser.add_argument("start_temp", help="initial temperature of the simulation")
parser.add_argument("end_temp", help="final temperature of the simulation")
parser.add_argument("nsteps", help="total number of steps of the simulation")
parser.add_argument("--npt", help="use MVLNPTMDSet for MD run in NPT ensemble", action="store_true")

args = parser.parse_args()
print(args)
####

structure = Structure.from_file(args.input)

if args.npt:
	mdset = MVLNPTMDSet(structure, args.start_temp, args.end_temp, args.nsteps)
else:
	mdset = MITMDSet(structure, args.start_temp, args.end_temp, args.nsteps)

mdset.write_input(args.output)
