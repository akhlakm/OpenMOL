#!/usr/bin/env python3

""" LAMMPS dat file writer for atom_style 'full'.

	This file is a part of OpenMOL python module.
	License GPLv3.0 Copyright (c) 2023 Akhlak Mahmood """

import re
import math
import openmol

# If box info is not set in openmol object
# we may need to estimate form the max and min value
# of the coordinates. This is the buffer for that.
BOX_BUFFER = 3.0	# A

def initialize():
	""" Initialize an empty openmol object with LAMMPS
		specific items. """

	MOL = openmol.initialize()
	MOL['source_format'] = "LAMMPS FULL"

	MOL['no_bond_types'] = 0
	MOL['no_angle_types'] = 0
	MOL['no_dihed_types'] = 0
	MOL['unique_atom_mass'] = []

	MOL['parm7_lj_epsilon'] = []
	MOL['parm7_lj_sigma'] = []

	return MOL


def build(MOL):
	""" Go through the openmol object and see if LAMMPS 
		specific items are properly calculated. If not,
		attemt to calculate them. """

	# add/update with default lammps items
	MOL = dict(initialize(), **MOL)

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

	# using parm7 data
	MOL['no_bond_types'] = len(MOL['FF_bond_k'])
	MOL['no_angle_types'] = len(MOL['FF_angle_k'])
	MOL['no_dihed_types'] = len(MOL['FF_dihed_k'])

	# if we have individual atom masses list, build type's masses
	if len(MOL['unique_atom_mass']) == 0 and len(MOL['atom_mass']):
		for unique_atom in MOL['unique_atom_types']:
			i = MOL['atom_type'].index(unique_atom)
			MOL['unique_atom_mass'].append(MOL['atom_mass'][i])

	if len(MOL['unique_atom_mass']) != MOL['no_atom_types']:
		print('-- LAMMPS Build Error: fail to build mass list, length mismatch.')

	# use PARM7: if we have A, B coeffs, build epsilon, sigma of parm7
	# @todo: move this to amber_parm7.py
	if len(MOL['parm7_lj_acoeff']) and len(MOL['parm7_lj_sigma']) != len(MOL['unique_atom_types']):
		if not MOL.get('parm7_lj_index', False):
			print('-- LAMMPS Build Error: non bonded parm7 indices not found.')

		else:
			# Amber specific way of finding out these values
			# See http://ambermd.org/formats.html
			for i in range(MOL['PARM_NTYPES']):
				j = MOL['parm7_lj_index'][i * (MOL['PARM_NTYPES'] + 1)] - 1

				A = MOL['parm7_lj_acoeff'][j]
				B = MOL['parm7_lj_bcoeff'][j]

				if A == 0.0:
					eps = 0.0
				else:
					eps = 0.25 * B**2 / A

				if B == 0.0:
					sigma = 0.0
				else:
					sigma = (A / B)**(1.0/6.0)

				MOL['parm7_lj_epsilon'].append(eps)
				MOL['parm7_lj_sigma'].append(sigma)

	if len(MOL['parm7_lj_epsilon']) != MOL['PARM_NTYPES'] or len(MOL['parm7_lj_sigma']) != MOL['PARM_NTYPES']:
		print('-- LAMMPS Build Error: fail to build parm7 lj coeffs, length mismatch.')

	# build indices of types
	if len(MOL['atom_type_index']) != MOL['no_atoms']:
		MOL['atom_type_index'] = []
		for i in range(MOL['no_atoms']):
			atom_type = MOL['atom_type'][i]
			MOL['atom_type_index'].append(MOL['unique_atom_types'].index(atom_type))

	if len(MOL['atom_type_index']) != MOL['no_atoms']:
		print('-- LAMMPS Build Error: fail to build atom type indices, length mismatch.')

	# assuming we have parm7 epsilon sigma lists built, build the lj params of each type
	if len(MOL['FF_lj_epsilon']) != MOL['no_atom_types'] or len(MOL['FF_lj_sigma']) != MOL['no_atom_types']:
		MOL['FF_lj_epsilon'] = []
		MOL['FF_lj_sigma'] = []
		for i in range(MOL['no_atom_types']):
			# get the first atom of this type
			aix = MOL['atom_type_index'].index(i)
			# get parm7 pair ff index of that atom
			pfx = MOL['pair_ff_index'][aix]
			MOL['FF_lj_epsilon'].append(MOL['parm7_lj_epsilon'][pfx])
			MOL['FF_lj_sigma'].append(MOL['parm7_lj_sigma'][pfx])

	if len(MOL['FF_lj_epsilon']) != MOL['no_atom_types'] or len(MOL['FF_lj_sigma']) != MOL['no_atom_types']:
		print('-- LAMMPS Build Error: fail to build pair coeffs, length mismatch.')

	MOL['_lammps_built'] = True
	return openmol.AttrDict(MOL)


class Reader(openmol.Reader):
	def __init__(self):
		super(Reader, self).__init__()
		self.Mol['source_format'] = "LAMMPS FULL"
		self.Mol['unique_atom_mass'] = []
		self.Mol['_lammps_built'] = False

	def read(self, lammps_data_file : str):
		super(Reader, self).read(lammps_data_file)
		self._process_lines()

	def _parse_str_as_type(self, string : str, dtype : callable, line, i):
		errstr  = f"-- Read Error: failed to parse {string} as {dtype}. "
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
			return self._parse_str_as_type(items[0], int, line, i)
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
				continue
			elif first_word == "Atoms":
				section = 'atom'
				continue
			elif first_word == "Bonds":
				section = 'bonds'
				continue
			elif first_word == "Angles":
				section = 'angles'
				continue
			elif first_word == "Dihedrals":
				section = 'diheds'
				continue
			elif first_word == "Impropers":
				section = 'impropers'
				continue

			if section == 'counts':
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
				masses = parts[0].split()
				comment = parts[1].strip() if len(parts) > 1 else None

				if parts[0] == "Atoms":
					section = 'atom'
					continue

				type_id = self._parse_str_as_type(masses[0], int, line, i)
				mass = self._parse_str_as_type(masses[1], float, line, i)

				self.Mol.unique_atom_mass.append(mass)

				if comment:
					self.Mol.unique_atom_types.append(comment)
				else:
					self.Mol.unique_atom_types.append(type_id - 1)

			elif section == 'atom':
				parts = line.split("#")
				info = parts[0].split()
				comment = parts[1].strip() if len(parts) > 1 else None
				
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



class Writer(openmol.Writer):
	def __init__(self, MOL, data_file):
		# open the file for writing
		super(Writer, self).__init__(MOL, data_file)

	def title(self):
		self.fp.write("%s (generated by OpenMOL lammps_full)\n\n" %self.MOL['title'])

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
