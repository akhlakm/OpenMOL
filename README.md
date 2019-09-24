# OpenMOL

OpenMOL is a python3 package that attempts to convert between popular molecular dynamics data file formats.

## Features

    • Convert AMBER PARM7 and Restart files into LAMMPS data files.
    • Read, manipulate and write MOL2 files.
    • Save data files into portable openmol json format without losing any properties.
    • Build complex mol2 systems using Discovery Studio and import into tleap.
    • Extensible.

## Current Implementations

Parser | Reader | Writer
-------|--------|-------
Amber  |  Yes   |  x
LAMMPS |   x    | Yes
MOL2   |  Yes   | Yes

## Usage Example: Amber to LAMMPS

```python
import sys
# add openmol to PYTHONPATH
sys.path.append("/path/to/openmol")

import openmol
import amber_parm7 as parm 
import lammps_full as lammps

# read amber parm files
p = parm.read('oleylamine.prmtop', 'oleylamine.rst7')

# calculate necessary lammps items and write
p = lammps.build(p)
lmp = lammps.Writer(p, 'data.oleylamine')
lmp.write()

# save everything as openmol json file
openmol.write_json(p, 'oleylamine.json')
```

## Usage Example: Fix VMD mol2

```python
import tripos_mol2 as mol2

# read vmd output mol2 file
p = mol2.read('vmd_saved.mol2')

# calculate necessary mol2 residue info
p = mol2.build(p)

# write the fixed mol2 file
mol2.Writer(p, 'fixed_vmd_saved.mol2').write()
```

## Usage Example: Fix DSV mol2 residue names

```python
import tripos_mol2 as mol2

p = mol2.read('dsv_output.mol2')

for i, atom in enumerate(p.atom_name):
	p.atom_resname[i] = 'FeO'

p = mol2.build(p)

# write the fixed mol2 file
mol2.Writer(p, 'fixed_dsv_output.mol2').write()
```

## License

GNU General Public License v3.0 (GPLv3.0)

Copyright (c) 2019 Akhlak Mahmood
