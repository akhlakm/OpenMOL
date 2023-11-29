#!/usr/bin/env python3

""" LAMMPS dat file reader and writer for atom_style 'full'.

	This file is a part of OpenMOL python module.
	License GPLv3.0 Copyright (c) 2023 Akhlak Mahmood """

import re
import math
from openmol import core

# If box info is not set in openmol object
# we may need to estimate form the max and min value
# of the coordinates. This is the buffer for that.
BOX_BUFFER = 3.0	# A

def initialize(new_items : dict = {}):
	""" Initialize an empty openmol object with LAMMPS
		specific items. """

	MOL = core.initialize()
	MOL['source_format'] = "LAMMPS FULL"
	MOL['no_bond_types'] = 0
	MOL['no_angle_types'] = 0
	MOL['no_dihed_types'] = 0
	MOL['unique_atom_mass'] = []
	MOL['pair_ff_index'] = []

	MOL.update(new_items)

	return MOL


def build(MOL):
	""" Go through the openmol object and see if LAMMPS 
		specific items are properly calculated. If not,
		attemt to calculate them. """

	# add/update with default lammps items
	MOL = initialize(MOL)

	# build residue id and name list
	if len(MOL['atom_resname']) == 0:
		for r, st in enumerate(MOL['residue_start']):
			res = MOL['residue_name'][r]

			# assume final atom is the last atom of residue
			end = MOL['no_atoms']

			# if there are more residues defined, update last atom
			if len(MOL['residue_start']) > r + 1:
				end = MOL['residue_start'][r+1]

			for i in range(st, end):
				MOL['atom_resid'].append(r)
				MOL['atom_resname'].append(res)

	# if we have individual atom masses list, build type's masses
	if len(MOL['unique_atom_mass']) == 0 and len(MOL['atom_mass']):
		for unique_atom in MOL['unique_atom_types']:
			i = MOL['atom_type'].index(unique_atom)
			MOL['unique_atom_mass'].append(MOL['atom_mass'][i])

	if len(MOL['unique_atom_mass']) != MOL['no_atom_types']:
		print('-- LAMMPS Build Error: fail to build mass list, length mismatch.')

	# build indices of types
	if len(MOL['atom_type_index']) != MOL['no_atoms']:
		MOL['atom_type_index'] = []
		for i in range(MOL['no_atoms']):
			atom_type = MOL['atom_type'][i]
			MOL['atom_type_index'].append(MOL['unique_atom_types'].index(atom_type))

	if len(MOL['atom_type_index']) != MOL['no_atoms']:
		print('-- LAMMPS Build Error: fail to build atom type indices, length mismatch.')

	# assuming we have parm7 epsilon sigma lists built, build the lj params of each type
	if len(MOL['FF_lj_epsilon']) != MOL['no_atom_types'] or \
			len(MOL['FF_lj_sigma']) != MOL['no_atom_types']:

		MOL['FF_lj_epsilon'] = []
		MOL['FF_lj_sigma'] = []
		for i in range(MOL['no_atom_types']):
			# get the first atom of this type
			aix = MOL['atom_type_index'].index(i)

			if len(MOL.pair_ff_index) > i:
				# get parm7 pair ff index of that atom
				pfx = MOL['pair_ff_index'][aix]
				MOL['FF_lj_epsilon'].append(MOL['parm7_lj_epsilon'][pfx])
				MOL['FF_lj_sigma'].append(MOL['parm7_lj_sigma'][pfx])

	if len(MOL['FF_lj_epsilon']) != MOL['no_atom_types'] or len(MOL['FF_lj_sigma']) != MOL['no_atom_types']:
		print('-- LAMMPS Build Error: fail to build pair coeffs, length mismatch.')

	MOL['_lammps_built'] = True
	print("Build done")
	return MOL


class Reader(core.Reader):
	def __init__(self):
		super(Reader, self).__init__()
		self.Mol['source_format'] = "LAMMPS FULL"
		self.Mol['unique_atom_mass'] = []
		self.Mol['_lammps_built'] = False

	def read(self, lammps_data_file : str):
		super(Reader, self).read_file(lammps_data_file)
		self._process_lines()
		print("-- WARN: FF params reading is not currently implemented")
		core.check(self.Mol)

	def build(self):
		self.Mol = build(self.Mol)

	def _parse_str_as_type(self, string : str, dtype : callable, line, i):
		errstr  = f"-- Read Error: failed to parse {string} as {dtype}, "
		errstr += f"line {i}: {line}"
		try:
			value = dtype(string)
		except TypeError:
			raise TypeError(errstr)
		return value

	def _section_starts(self, line, i, startswith) -> int:
		""" Return count if a line denotes a section start.
		Ex. atom count, or atom type count.
		"""
		items = re.findall(rf'^(\d+)\s+{startswith}$', line)
		if items:
			return self._parse_str_as_type(items[0], int, line, i+1)
		else:
			return -1

	def _process_lines(self):
		line_no = 0
		section_line_no = 0
		section = None

		for i, line in enumerate(self.lines):
			line_no += 1
			section_line_no += 1
			line = line.strip()

			if len(line) == 0:
				continue

			# comments
			if line.startswith("#"):
				continue

			first_word = line.split()[0]

			# title, can be optional
			if i == 0:
				self.Mol.title = line
				continue

			# count section start
			if section is None:
				at = self._section_starts(line, i, "atoms")
				if at >= 0:
					self.Mol.no_atoms = at
					section = 'counts'
					continue

			if first_word == "Masses":
				section = 'mass'
				print('Reading mass list ...')
				continue
			elif first_word == "Atoms":
				section = 'atom'
				print('Reading atom list ...')
				continue
			elif first_word == "Bonds":
				section = 'bond'
				print('Reading bond list ...')
				continue
			elif first_word == "Angles":
				section = 'angle'
				print('Reading angle list ...')
				continue
			elif first_word == "Dihedrals":
				section = 'dihed'
				print('Reading dihedral list ...')
				continue
			elif first_word == "Impropers":
				section = 'improper'
				print('Reading improper torsion list ...')
				continue

			elif line.endswith('Coeffs'):
				print('-- Warning: Coeffs reading not implemented yet:', line)
				section = None
				continue

			if section is None:
				# ignore the coeff sections for now.
				continue

			elif section == 'counts':
				counts = re.findall(r'^(\d+)\s+([a-z]+)s$', line)
				if counts:
					counts = counts[0] # get the tuple/first match
					number = self._parse_str_as_type(counts[0], int, line, i)
					item = counts[1]
					if item == 'bond':
						self.Mol.no_bonds = number
					elif item == 'angle':
						self.Mol.no_angles = number
					elif item == 'dihedral':
						self.Mol.no_diheds = number
					elif item == 'improper':
						self.Mol.no_improper = number
					else:
						print('-- Read Error: unknown count item %s (line %d)' %(item, i+1))
						return
				else:
					at = self._section_starts(line, i, "atom types")
					if at >= 0:
						section = 'types'
						self.Mol.no_atom_types = at

			elif section == 'types':
				counts = re.findall(r'^(\d+)\s+([a-z]+)\s+types$', line)
				if counts:
					counts = counts[0] # get the tuple/first match
					number = self._parse_str_as_type(counts[0], int, line, i)
					item = counts[1]
					if item == 'bond':
						self.Mol.no_bond_types = number
					elif item == 'angle':
						self.Mol.no_angle_types = number
					elif item == 'dihedral':
						self.Mol.no_dihed_types = number
					elif item == 'improper':
						self.Mol.no_improper_types = number
					else:
						print('-- Read Error: unknown count item %s (line %d)' %(item, i+1))
						return

				else:
					boxsize = re.findall(r'([+-]?[0-9]*[.]?[0-9]+)\s+([+-]?[0-9]*[.]?[0-9]+)\s+xlo\s+xhi$', line)
					if boxsize:
						section = 'boxsize'
						boxsize = boxsize[0] # get the tuple/first match
						low = self._parse_str_as_type(boxsize[0], float, line, i)
						high = self._parse_str_as_type(boxsize[1], float, line, i)
						self.Mol.box_x_low = low
						self.Mol.box_x_high = high
						self.Mol.box_x = high - low

			elif section == 'boxsize':
				parts = line.split()
				if parts[0] == "Masses":
					section = 'mass'
					continue

				assert len(parts) >= 4, \
					"Invalid box info, line %d: %s" %(i+1, line)

				low = self._parse_str_as_type(parts[0], float, line, i)
				high = self._parse_str_as_type(parts[1], float, line, i)
				if parts[2] == 'ylo':
					self.Mol.box_y = high - low
					self.Mol.box_y_high = high
					self.Mol.box_y_low = low
				elif parts[2] == 'zlo':
					self.Mol.box_z = high - low
					self.Mol.box_z_high = high
					self.Mol.box_z_low = low
				else:
					errstr  = f"-- Read Error: failed to parse boxsize. "
					errstr += f"line {i}: {line}"
					ValueError(errstr)

			elif section == 'mass':
				parts = line.split("#")
				info = parts[0].split()
				comment = parts[1].strip() if len(parts) > 1 else None

				if parts[0] == "Atoms":
					section = 'atom'
					continue

				assert len(info) >= 2, \
					"Invalid mass info, line %d: %s" %(i+1, line)

				type_id = self._parse_str_as_type(info[0], int, line, i)
				mass = self._parse_str_as_type(info[1], float, line, i)

				self.Mol.unique_atom_mass.append(mass)

				if comment:
					self.Mol.unique_atom_types.append(comment)
				else:
					self.Mol.unique_atom_types.append(type_id - 1)

			elif section == 'atom':
				parts = line.split("#")
				info = parts[0].split()
				comment = parts[1].strip() if len(parts) > 1 else None
				
				if parts[0] == "Bonds":
					section = 'bond'
					continue

				assert len(info) >= 7, \
					"Invalid atom info, line %d: %s" %(i+1, line)

				res_id = self._parse_str_as_type(info[1], int, line, i)
				at_type = self._parse_str_as_type(info[2], int, line, i)
				
				at_q = self._parse_str_as_type(info[3], float, line, i)
				at_x = self._parse_str_as_type(info[4], float, line, i)
				at_y = self._parse_str_as_type(info[5], float, line, i)
				at_z = self._parse_str_as_type(info[6], float, line, i)

				self.Mol.atom_x.append(at_x)
				self.Mol.atom_y.append(at_y)
				self.Mol.atom_z.append(at_z)
				self.Mol.atom_q.append(at_q)
				self.Mol.atom_type_index.append(at_type - 1)
				self.Mol.atom_resid.append(res_id - 1)
				self.Mol.atom_resname.append(res_id)

				if comment:
					comment = comment.strip()
					self.Mol.atom_name.append(comment)
					self.Mol.atom_type.append(comment)

			elif section == 'bond':
				parts = line.split("#")
				info = parts[0].split()
				comment = parts[1].strip() if len(parts) > 1 else None
				
				if parts[0] == "Angles":
					section = 'angle'
					continue

				assert len(info) >= 4, \
					"Invalid bond info, line %d: %s" %(i+1, line)

				bond_type = self._parse_str_as_type(info[1], int, line, i)
				bond_from_atom = self._parse_str_as_type(info[2], int, line, i)
				bond_to_atom = self._parse_str_as_type(info[3], int, line, i)

				self.Mol.bond_ff_index.append(bond_type - 1)
				self.Mol.bond_from.append(bond_from_atom - 1)
				self.Mol.bond_to.append(bond_to_atom - 1)

			elif section == 'angle':
				parts = line.split("#")
				info = parts[0].split()
				comment = parts[1].strip() if len(parts) > 1 else None
				
				if parts[0] == "Dihedrals":
					section = 'dihed'
					continue

				assert len(info) >= 4, \
					"Invalid angle info, line %d: %s" %(i+1, line)

				angle_type = self._parse_str_as_type(info[0], int, line, i)
				angle_a = self._parse_str_as_type(info[1], int, line, i)
				angle_b = self._parse_str_as_type(info[2], int, line, i)
				angle_c = self._parse_str_as_type(info[3], int, line, i)

				self.Mol.angle_a.append(angle_a - 1)
				self.Mol.angle_b.append(angle_b - 1)
				self.Mol.angle_c.append(angle_c - 1)
				self.Mol.angle_ff_index.append(angle_type - 1)

			elif section == 'dihed':
				parts = line.split("#")
				info = parts[0].split()
				comment = parts[1].strip() if len(parts) > 1 else None
				
				if parts[0] == "Impropers":
					section = 'improper'
					continue

				assert len(info) >= 5, \
					"Invalid angle info, line %d: %s" %(i+1, line)

				dihed_type = self._parse_str_as_type(info[0], int, line, i)
				dihed_a = self._parse_str_as_type(info[1], int, line, i)
				dihed_b = self._parse_str_as_type(info[2], int, line, i)
				dihed_c = self._parse_str_as_type(info[3], int, line, i)
				dihed_d = self._parse_str_as_type(info[4], int, line, i)

				self.Mol.dihed_a.append(dihed_a - 1)
				self.Mol.dihed_b.append(dihed_b - 1)
				self.Mol.dihed_c.append(dihed_c - 1)
				self.Mol.dihed_d.append(dihed_d - 1)
				self.Mol.dihed_ff_index.append(dihed_type - 1)

			elif section == 'improper':
				parts = line.split("#")
				info = parts[0].split()
				comment = parts[1].strip() if len(parts) > 1 else None
				
				if len(parts) == 1:
					section = None
					continue

				assert len(info) >= 5, \
					"Invalid angle info, line %d: %s" %(i+1, line)

				improper_type = self._parse_str_as_type(info[0], int, line, i)
				improper_a = self._parse_str_as_type(info[1], int, line, i)
				improper_b = self._parse_str_as_type(info[2], int, line, i)
				improper_c = self._parse_str_as_type(info[3], int, line, i)
				improper_d = self._parse_str_as_type(info[4], int, line, i)

				self.Mol.improper_a.append(improper_a - 1)
				self.Mol.improper_b.append(improper_b - 1)
				self.Mol.improper_c.append(improper_c - 1)
				self.Mol.improper_d.append(improper_d - 1)
				self.Mol.improper_ff_index.append(improper_type - 1)
			else:
				print("-- WARN: Unknown section, line %d: %s" %(i+1, line))

		print("Read OK")


class Writer(core.Writer):
	def __init__(self, MOL, data_file):
		# open the file for writing
		super(Writer, self).__init__(MOL, data_file)

	def title(self):
		self.MOL['title'] = self.MOL['title'].replace("\n", " ")
		self.fp.write("%s (by OpenMOL)\n\n" %self.MOL['title'])

	def counts(self):
		self.fp.write("%d atoms\n" %self.MOL['no_atoms'])
		self.fp.write("%d bonds\n" %self.MOL['no_bonds'])
		self.fp.write("%d angles\n" %self.MOL['no_angles'])
		self.fp.write("%d dihedrals\n" %self.MOL['no_diheds'])
		self.fp.write("0 impropers\n\n")	# todo: fix it

	def types(self):
		self.fp.write("%d atom types\n" %self.MOL['no_atom_types'])
		self.fp.write("%d bond types\n" %self.MOL['no_bond_types'])
		self.fp.write("%d angle types\n" %self.MOL['no_angle_types'])
		self.fp.write("%d dihedral types\n\n" %self.MOL['no_dihed_types'])

	def box_info(self):
		# @todo: handle situations where coords have -ve values
		if self.MOL['box_x'] == 0.0:
			xlo = min(self.MOL['atom_x']) - BOX_BUFFER
			xhi = max(self.MOL['atom_x']) + BOX_BUFFER
		else:
			xlo = 0.0
			xhi = self.MOL['box_x']

		if self.MOL['box_y'] == 0.0:
			ylo = min(self.MOL['atom_y']) - BOX_BUFFER
			yhi = max(self.MOL['atom_y']) + BOX_BUFFER
		else:
			ylo = 0.0
			yhi = self.MOL['box_y']

		if self.MOL['box_z'] == 0.0:
			zlo = min(self.MOL['atom_z']) - BOX_BUFFER
			zhi = max(self.MOL['atom_z']) + BOX_BUFFER
		else:
			zlo = 0.0
			zhi = self.MOL['box_z']

		self.fp.write("%8.4f %8.4f xlo xhi\n" %(xlo, xhi))
		self.fp.write("%8.4f %8.4f ylo yhi\n" %(ylo, yhi))
		self.fp.write("%8.4f %8.4f zlo zhi\n" %(zlo, zhi))

	def masses(self):
		self.fp.write("\nMasses\n\n")
		for i in range(self.MOL['no_atom_types']):
			self.fp.write('%3d  %6.3f   # %s\n'
				%(i+1, self.MOL['unique_atom_mass'][i], self.MOL['unique_atom_types'][i]))

	def pair_coeffs(self):
		self.fp.write("\nPair Coeffs\n\n")
		for i in range(self.MOL['no_atom_types']):
			self.fp.write('%3d  %10.4f   %10.4f   # %s\n'
				%(i+1, self.MOL['FF_lj_epsilon'][i], self.MOL['FF_lj_sigma'][i], self.MOL['unique_atom_types'][i]))

	def bond_coeffs(self):
		self.fp.write("\nBond Coeffs\n\n")
		for i in range(self.MOL['no_bond_types']):
			self.fp.write('%3d  %6.3f   %6.3f\n'
				%(i+1, self.MOL['FF_bond_k'][i], self.MOL['FF_bond_eq'][i]))

	def angle_coeffs(self):
		self.fp.write("\nAngle Coeffs\n\n")
		for i in range(self.MOL['no_angle_types']):
			self.fp.write('%3d  %6.3f  %6.3f\n'
				%(i+1, self.MOL['FF_angle_k'][i], math.degrees(self.MOL['FF_angle_eq'][i])))

	def dihed_coeffs(self):
		self.fp.write("\nDihedral Coeffs\n\n")
		for i in range(self.MOL['no_dihed_types']):
			phase = int(math.cos(self.MOL['FF_dihed_phase'][i]))

			if phase == 0:
				phase = 1
			else: phase = -1

			period = int(self.MOL['FF_dihed_periodicity'][i])
			self.fp.write('%3d  %6.3f  %2d  %d\n'
				%(i+1, self.MOL['FF_dihed_k'][i], phase, period))

	def atoms(self):
		self.fp.write("\nAtoms # atom_style_full\n\n")

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
				'type': self.MOL['atom_type'][i]
			}

			atomstr =	"{id:>7d} {resid:>4d} {typeid:>3} {charge:>10.6f}  " \
						"{x:>8.4f}  {y:>8.4f}  {z:>8.4f} # {type}\n"

			self.fp.write(atomstr.format(**atom))

	def bonds(self):
		self.fp.write("\nBonds\n\n")

		for i in range(self.MOL['no_bonds']):
			bond = {
				'id' : i + 1,
				'from': self.MOL['bond_from'][i] + 1,
				'to': self.MOL['bond_to'][i] + 1,
				'type': self.MOL['bond_ff_index'][i] + 1,
			}
			bondstr =	"{id:>7d}  {type:>5d}  {from:>7d}  {to:>7d} \n"
			self.fp.write(bondstr.format(**bond))

	def angles(self):
		self.fp.write("\nAngles\n\n")

		for i in range(self.MOL['no_angles']):
			angle = {
				'id' : i + 1,
				'a': self.MOL['angle_a'][i] + 1,
				'b': self.MOL['angle_b'][i] + 1,
				'c': self.MOL['angle_c'][i] + 1,
				'type': self.MOL['angle_ff_index'][i] + 1,
			}
			anglestr =	"{id:>7d}  {type:>3d}  {a:>5d}  {b:>5d}  {c:>5d} \n"
			self.fp.write(anglestr.format(**angle))


	def diheds(self):
		self.fp.write("\nDihedrals\n\n")

		for i in range(self.MOL['no_diheds']):
			dihed = {
				'id' : i + 1,
				'a': self.MOL['dihed_a'][i] + 1,
				'b': self.MOL['dihed_b'][i] + 1,
				'c': self.MOL['dihed_c'][i] + 1,
				'd': self.MOL['dihed_d'][i] + 1,
				'type': self.MOL['dihed_ff_index'][i] + 1,
			}
			dihedstr = "{id:>7d}  {type:>3d}  {a:>5d}  {b:>5d}  {c:>5d}  {d:>5d} \n"
			self.fp.write(dihedstr.format(**dihed))


	def write(self):
		if not self.MOL.get('_lammps_built', False):
			print('-- Warning: MOL not getting build() for LAMMPS likely to fail while writing.')

		self.title()
		self.counts()
		self.types()
		self.box_info()
		self.masses()
		self.pair_coeffs()
		self.bond_coeffs()
		self.angle_coeffs()
		self.dihed_coeffs()
		self.atoms()
		self.bonds()
		self.angles()
		self.diheds()

		# close the file
		super(Writer, self).write()
