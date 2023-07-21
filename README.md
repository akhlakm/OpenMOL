# OpenMOL
OpenMOL is a Python package to convert between popular molecular dynamics data file formats. It exposes a pythonic way to easily manipulate MOL2 and LAMMPS data files.

## Features
- Convert AMBER PARM7 and restart files into LAMMPS data files.
- Read, manipulate and write MOL2 files.
- Save data files into portable openmol json format without losing any properties.
- Build complex mol2 systems using Discovery Studio and import into *tleap*.
- Extensible.

## Installation
You can install `openmol` directly from [PyPI](https://pypi.org/project/openmol) using `pip`.

```sh
pip install openmol
```

## Usage Examples

### AMBER to LAMMPS
Convert AMBER **prmtop** and **restart** files to a LAMMPS **data** file.
```python
from openmol import amber_parm7 as amber
from openmol import lammps_full as lammps

# read amber parm files
p = amber.read('system.prmtop', 'system.rst7')

# calculate necessary lammps items and write
p = lammps.build(p)
lmp = lammps.Writer(p, 'system.data')
lmp.write()
```

### Fix VMD MOL2
The *SUBSTRUCTURE* section of the MOL2 files created by VMD are not always properly formated. Use the following to fix the format.

```python
from openmol import tripos_mol2 as mol2

# read vmd output mol2 file
p = mol2.read('vmd_saved.mol2')

# calculate necessary mol2 residue info
p = mol2.build(p)

# write the fixed mol2 file
mol2.Writer(p, 'fixed_vmd_saved.mol2').write()
```

### Fix DSV mol2 residue names
Discovery Studio Visualizer generated MOL2 files can contain digits as suffixes in the residue names. Update the residue names for all atoms as below.

```python
from openmol import tripos_mol2 as mol2

p = mol2.read('dsv_output.mol2')

# set all atom's resname to FeO
for i, atom in enumerate(p.atom_name):
	p.atom_resname[i] = 'FeO'

p = mol2.build(p)

# write the fixed mol2 file
mol2.Writer(p, 'fixed_dsv_output.mol2').write()
```

## Current Implementations

Parser | Reader | Writer
-------|--------|-------
Amber  |  Yes   |  x
LAMMPS |   x    | Yes
MOL2   |  Yes   | Yes

*Contributions welcome!*

For the latest development version you can clone the [Git repository](https://github.com/akhlakm/OpenMOL).

```sh
git clone https://github.com/akhlakm/OpenMOL
cd OpenMOL
pip install -e .
```

## License
GNU General Public License v3.0 (GPLv3.0)
