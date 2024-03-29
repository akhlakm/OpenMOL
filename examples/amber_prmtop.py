#!/usr/bin/env python3

""" Simple example using OpenMOL to convert from
	AMBER PARM files into LAMMPS dat file.

	This file is a part of OpenMOL python module.
	License GPLv3.0 Copyright (c) 2019 Akhlak Mahmood """


# oleylamine example

import sys
sys.path.append("..")

import openmol
from openmol import amber_parm7 as parm 
from openmol import lammps_qmag as lammps
from openmol import tripos_mol2 as mol2

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

