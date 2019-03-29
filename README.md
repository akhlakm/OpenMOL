# OpenMOL

OpenMOL is a python3 package that attempts to convert between popular molecular dynamics data file formats.

## Features

    • Convert AMBER PARM7 and Restart files into LAMMPS data file.
    • Read and write MOL2 file.
    • Save data files into portable openmol json format without losing any properties.
    • Extensible.

## Current Implementations

Parser | Reader | Writer
-------|--------|-------
Amber  |  Yes   |  x
LAMMPS |   x    | Yes
MOL2   |  Yes   | Yes

## Usage: Amber to LAMMPS

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

Please refer to the User Guide for more examples and details.

## License

GNU General Public License v3.0 (GPLv3.0)
Copyright (c) 2019 Akhlak Mahmood
