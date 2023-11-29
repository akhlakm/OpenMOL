#!/usr/bin/env python3

""" LAMMPS dat file writer for atom_style 'qmag'.

	This file is a part of OpenMOL python module.
	License GPLv3.0 Copyright (c) 2023 Akhlak Mahmood """

from openmol import core
from . import lammps_full as lmp 


def initialize():
	""" Initialize an empty openmol object with LAMMPS qmag
		specific items. """

	MOL = core.initialize()
	MOL = dict(lmp.initialize(), **MOL)
	MOL['source_format'] = "LAMMPS QMAG"
	MOL['atom_qm'] = []
	return core.AttrDict(MOL)


def build(MOL):
	MOL = lmp.build(MOL)
	MOL = dict(initialize(), **MOL)

	if len(MOL['atom_qm']) != MOL['no_atoms']:
		MOL['atom_qm'] = [0.0 for i in range(MOL['no_atoms'])]

	MOL['_lammps_qmag_built'] = True
	return core.AttrDict(MOL)


def qm_for_index(MOL, ix, qm):
	if not MOL.get('_lammps_qmag_built', False):
		MOL = build(MOL)

	MOL['atom_qm'][ix] = qm
	return MOL


def print_qm(MOL):
	if not MOL.get('atom_qm', False):
		print("-- No qm data set for MOL.")
		return

	if not MOL.get('_lammps_qmag_built', False):
		print('-- Warning: MOL not build() for QMAG likely to fail while writing.')

	for i, qm in enumerate(MOL['atom_qm']):
		if qm != 0.0:
			atom = {
				'id' : i + 1,
				'name': MOL['atom_name'][i],
				'type': MOL['atom_type_index'][i] + 1,
				'resid': MOL['atom_resid'][i] + 1,
				'resname': MOL['atom_resname'][i],
				'charge': MOL['atom_q'][i],
				'qm': MOL['atom_qm'][i]
			}

			atomstr = "{id:>7d} {name:>3} {type:>3} {resid:>4d} {resname:>4} " \
					  "{charge:>11.6f} {qm:>7.4f}\n"

			print(atomstr.format(**atom))


class Writer(lmp.Writer):

	def __init__(self, MOL, data_file):
		super(Writer, self).__init__(MOL, data_file)

	def title(self):
		self.fp.write("%s \n\n" %self.MOL['title'])

	def atoms(self):
		self.fp.write("\nAtoms # atom_style_qmag\n\n")

		for i in range(self.MOL['no_atoms']):
			resid = self.MOL['atom_resid'][i]
			atom = {
				'id' : i + 1,
				'name': self.MOL['atom_name'][i],
				'x': self.MOL['atom_x'][i],
				'y': self.MOL['atom_y'][i],
				'z': self.MOL['atom_z'][i],
				'typeid': self.MOL['atom_type_index'][i] + 1,
				'resid': self.MOL['atom_resid'][i] + 1,
				'charge': self.MOL['atom_q'][i],
				'type': self.MOL['atom_type'][i],
				'qm': self.MOL['atom_qm'][i]
			}

			atomstr =	"{id:>7d} {resid:>4d} {typeid:>3} {charge:>10.6f}  " \
						"{x:>8.4f}  {y:>8.4f}  {z:>8.4f}   {qm:>7.4f} # {type}\n"

			self.fp.write(atomstr.format(**atom))

	def write(self):
		if not self.MOL.get('_lammps_qmag_built', False):
			print('-- Warning: MOL not build() for QMAG likely to fail while writing.')

		if not self.MOL.get('atom_qm', False):
			print('-- Error: No qm data found for the atoms.')
			return False

		if len(self.MOL['atom_qm']) != self.MOL['no_atoms']:
			print('-- Error: Not enough qm data for all the atoms.')
			return False

		super(Writer, self).write()
