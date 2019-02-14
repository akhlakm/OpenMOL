#!/usr/bin/env python3

# oleylamine example

import sys
sys.path.append("..")

import openmol
import amber_parm7 as parm 
import lammps_qmag as lammps
import tripos_mol2 as mol2

# read amber parm files
p = parm.read('oleylamine.prmtop', 'oleylamine.rst7')

# calculate necessary lammps items and write
p = lammps.build(p)
lmp = lammps.Writer(p, 'data.oleylamine')
lmp.write()

# calculate necessary mol2 items and write
p = mol2.build(p)
mol = mol2.Writer(p, 'oleylamine.mol2')
mol.write()

# save everything as openmol json file
openmol.write_json(p, 'oleylamine.json')

