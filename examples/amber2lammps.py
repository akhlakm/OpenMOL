#!/usr/bin/env python3

# oleylamine example

import sys
sys.path.append("..")

import openmol
import amber_parm7 as parm 
import lammps_full as lmp 


# read amber parm files
p = parm.read('oleylamine.prmtop', 'oleylamine.rst7')

# calculate necessary lammps data items
p = lmp.build(p)

# write lammps data
lmp.write(p, 'data.oleylamine')

# save as openmol json file
openmol.write_json(p, 'oleylamine.json')

