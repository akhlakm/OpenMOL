#!/usr/bin/env python3

""" TRIPOS MOL2 file reader and writer.

	This file is a part of OpenMOL python module.
	License GPLv3.0 Copyright (c) 2023 Akhlak Mahmood """

from openmol import core


def initialize(new_items : dict = {}):
	""" Initialize an openmol object with TRIPOS MOL2
		specific items. """

	MOL = core.initialize()

	MOL['source_format'] = 'TRIPOS MOL2'

	MOL['atom_status_bit'] = []
	MOL['residue_dict_type'] = []
	MOL['residue_chain'] = []
	MOL['residue_sub_type'] = []
	MOL['residue_comment'] = []
	MOL['residue_inter_bonds'] = []
	MOL['residue_status_bits'] = []

	MOL.update(new_items)

	return MOL


def build(MOL):
	""" Go through the openmol object and see if MOL2
		specific items are properly calculated.
		If not, attemt to determine/guess them. """

	# add/update with default mol2 items
	MOL = initialize(MOL)

	if not MOL['type']:
		MOL['type'] = 'SMALL'

	if not MOL['charge_type']:
		MOL['charge_type'] = 'USER_CHARGES'

	# Fix residue list
	# case when resnames are defined in substruct section
	# but not in atoms
	if len(MOL['atom_resname']) == 0:
		# print("-- Warning: atom residue info missing. Building from residue list.\n")
		for r, st in enumerate(MOL['residue_start']):
			res = MOL['residue_name'][r]

			end = MOL['no_atoms']
			if len(MOL['residue_start']) > r + 1:
				end = MOL['residue_start'][r+1]

			for i in range(st, end):
				MOL['atom_resid'].append(r)
				MOL['atom_resname'].append(res)

	unique_resids = list(set(MOL['atom_resid']))

	# Fix atom resids order
	if len(unique_resids) > 1:
		old_id = 0
		id_num = 0
		for i, resid in enumerate(MOL['atom_resid']):
			if old_id > resid:
				# resid should increase all the times
				print('-- Error: resid order invalid at atom %d' %(i+1))
				print('Previous atom', old_id+1, 'Current atom', resid+1)
				return None

			# new residue begins
			if resid != old_id:
				old_id = resid
				id_num += 1
			MOL['atom_resid'][i] = id_num

	# if no type set, use single bond
	if len(MOL['bond_type']) == 0:
		print("-- Warning: bond type info missing. Assuming single bonds.\n")
		for i in range(MOL['no_bonds']):
			MOL['bond_type'].append("1")

	# Fix residue list
	# case when resnames/ids are defined in atoms section
	if len(MOL['residue_start']) != len(unique_resids):
		# print("-- Warning: residue list trancated. Building from atoms residue info.\n")
		current_id = 0
		MOL['residue_start'] = []
		MOL['residue_type'] = []
		MOL['residue_name'] = []

		# add the first atom, 0 based indexing
		MOL['residue_start'].append(0)
		MOL['residue_name'].append(MOL['atom_resname'][0])

		for i, resid in enumerate(MOL['atom_resid']):
			# new residue begins
			if resid != current_id:
				current_id = resid
				MOL['residue_start'].append(i)
				MOL['residue_name'].append(MOL['atom_resname'][i])

	MOL['no_residues'] = len(MOL['residue_name'])
	MOL['no_atoms'] = len(MOL['atom_name'])
	MOL['no_bonds'] = len(MOL['bond_from'])

	if len(unique_resids) != MOL['no_residues']:
		print("-- Error: unique_resids", len(unique_resids), 'no_residues', MOL['no_residues'])
		return None

	if len(MOL['residue_type']) == 0:
		for i in range(MOL['no_residues']):
			MOL['residue_type'].append("RESIDUE")

	MOL['_mol2_built'] = True
	print("Build done")
	return MOL


def check_last_section(section, MOL):
	""" Parse and process the last read section """

	if section == 'ATOM':
		if not core.check_atoms_ok(MOL):
			# @todo: do this check here
			return False
		else:
			unique_atom_types = set(MOL['atom_type'])
			MOL['unique_atom_types'] = list(unique_atom_types)

	elif section == 'BOND':
		if not core.check_bonds_ok(MOL):
			# @todo: do this check here
			return False

	elif section == 'SUBSTRUCTURE':
		if not core.check_residues_ok(MOL):
			# @todo: do this check here
			return False

	# print('OK')
	return True


def read(mol2_file):
	""" Read a TRIPOS MOL2 file and store as OpenMOL object.
		All indices are decremented to use 0 base indexing. """

	MOL = initialize()

	line_no = 0
	section_line_no = 0
	section = None

	# variables for sanity checking
	type_ok = False
	charge_ok = False
	name_ok = False
	summary_ok = False

	print("\nReading:", mol2_file, end=' ... ')

	for line in open(mol2_file, 'r'):

		line_no += 1
		section_line_no += 1

		line = line.strip()

		if len(line) == 0:
			continue

		# comments
		if line.startswith('#'):
			MOL['description'] += "%s\n" %line
			continue

		# new section definition
		if line.startswith('@<'):
			parts = line.split('>')

			if len(parts) < 2:
				print('-- Error: Invalid MOL2 [line %d]:\n%s' %(line_no, line))
				return None
			else:
				if not check_last_section(section, MOL):
					return False

				section = parts[1]
				# print('Reading %s ...' %section, end=' ')
				section_line_no = 0
				continue

		# handle sections
		if section == 'MOLECULE':
			if section_line_no == 1:
				MOL['title'] = line 
				name_ok = True

			elif section_line_no == 2:
				parts = line.split()
				if len(parts) > 0:
					MOL['no_atoms'] = int(parts[0])

				if len(parts) > 1:
					MOL['no_bonds'] = int(parts[1])

				if len(parts) > 2:
					MOL['no_residues'] = int(parts[2])

				if len(parts) > 3:
					MOL['no_features'] = int(parts[3])

				if len(parts) > 4:
					MOL['no_sets'] = int(parts[4])

				summary_ok = True

			elif section_line_no == 3:
				MOL['type'] = line 
				type_ok = True

			elif section_line_no == 4:
				MOL['charge_type'] = line
				charge_ok = True

			elif section_line_no == 5:
				MOL['mol2_status_bits'] = line

			elif section_line_no == 6:
				MOL['mol2_comment'] = line

		elif section == 'ATOM':
			parts = line.split()

			if len(parts) < 6:
				print('-- Error: Invalid MOL2 [line %d]:\n%s' %(line_no, line))
				return None

			# mandatory items
			MOL['atom_name'].append(parts[1])
			MOL['atom_x'].append(float(parts[2]))
			MOL['atom_y'].append(float(parts[3]))
			MOL['atom_z'].append(float(parts[4]))
			MOL['atom_type'].append(parts[5])

			# optional items
			if len(parts) > 6:
				MOL['atom_resid'].append(int(parts[6]) - 1)

			if len(parts) > 7:
				MOL['atom_resname'].append(parts[7])

			if len(parts) > 8:
				charge = round(float(parts[8]), 4)
				MOL['atom_q'].append(charge)

			if len(parts) > 9:
				MOL['atom_status_bit'].append(parts[9])


		elif section == 'BOND':
			parts = line.split()

			if len(parts) < 4:
				print('-- Error: Invalid MOL2 [line %d]:\n%s' %(line_no, line))
				return None

			MOL['bond_from'].append(int(parts[1]) - 1)
			MOL['bond_to'].append(int(parts[2]) - 1)
			MOL['bond_type'].append(parts[3])

			if len(parts) > 4:
				MOL['bond_status_bit'].append(parts[4])

		elif section == 'SUBSTRUCTURE':
			parts = line.split()

			if len(parts) < 3:
				print('-- Error: Invalid MOL2 [line %d]:\n%s' %(line_no, line))
				return None

			MOL['residue_name'].append(parts[1])
			MOL['residue_start'].append(int(parts[2]) - 1)

			if len(parts) > 3:
				MOL['residue_type'].append(parts[3])
			if len(parts) > 4:
				MOL['residue_dict_type'].append(parts[3])
			if len(parts) > 5:
				MOL['residue_chain'].append(parts[3])
			if len(parts) > 6:
				MOL['residue_sub_type'].append(parts[3])
			if len(parts) > 7:
				MOL['residue_inter_bonds'].append(parts[3])
			if len(parts) > 8:
				MOL['residue_status_bits'].append(parts[3])
			if len(parts) > 9:
				MOL['residue_comment'].append(parts[3])

		else:
			# unknown section
			# @todo: extend here if needed
			print('-- Ignored unknown MOL2 section: %s ' %section)

	if not check_last_section(section, MOL):
		return False

	print('Done.')
	return MOL


class Writer(core.Writer):
	""" 	Write a TRIPOS MOL2 file using the OpenMOL object data.
			All indices are incremented by 1 to follow 1 based
			indexing in MOL2 format. """

	def __init__(self, MOL, data_file):
		# open the file
		super(Writer, self).__init__(MOL, data_file)

	def molecule(self):
		self.fp.write('@<TRIPOS>MOLECULE\n')
		molecstr = 	"{title}\n" \
					"{no_atoms:>5d} {no_bonds:>5d} {no_residues:>5d} 0 0\n" \
					"{type}\n" \
					"{charge_type}\n\n"

		self.fp.write(molecstr.format(**self.MOL))

	def atoms(self):
		self.fp.write('@<TRIPOS>ATOM\n')

		for i in range(self.MOL['no_atoms']):
			resid = self.MOL['atom_resid'][i]
			atom = {
				'id' : i + 1,
				'name': self.MOL['atom_name'][i],
				'x': self.MOL['atom_x'][i],
				'y': self.MOL['atom_y'][i],
				'z': self.MOL['atom_z'][i],
				'type': self.MOL['atom_type'][i],
				'resid': self.MOL['atom_resid'][i] + 1,
				'resname': self.MOL['residue_name'][resid],
				'charge': self.MOL['atom_q'][i]
			}
			atomstr =	"{id:>7d} {name:<5}  " \
						"{x:>8.4f}  {y:>8.4f}  {z:>8.4f}   {type:>3} " \
						"{resid:>3} {resname:<5}   {charge:>11.4f}\n"
			self.fp.write(atomstr.format(**atom))

	def bonds(self):
		self.fp.write('@<TRIPOS>BOND\n')
		for i in range(self.MOL['no_bonds']):
			bond = {
				'id' : i + 1,
				'from': self.MOL['bond_from'][i] + 1,
				'to': self.MOL['bond_to'][i] + 1,
				'type': self.MOL['bond_type'][i],
			}

			# Aromatic bonds
			if bond['type'] in [1.5, 'ar']:
				bond['type'] = 'ar'

			bondstr =	"{id:>7d}  {from:>7d}  {to:>7d}   {type:>3} \n"
			self.fp.write(bondstr.format(**bond))

	def substructures(self):
		self.fp.write('@<TRIPOS>SUBSTRUCTURE\n')
		for i in range(self.MOL['no_residues']):
			resid = self.MOL['atom_resid'][i]
			subst = {
				'id' : i + 1,
				'name': self.MOL['residue_name'][i],
				'root': self.MOL['residue_start'][i] + 1,
				'type': self.MOL['residue_type'][i],
			}
			resstr =	"{id:>7d}  {name:>7}  {root:>7d}   {type:>7} \n"
			self.fp.write(resstr.format(**subst))

	# @extend: add additional sections if needed

	def write(self):
		# if not self.MOL.get('_mol2_built', False):
		# 	print('-- Warning: call tripos_mol2.build(MOL) before writing. Continuing anyway ...')

		self.molecule()
		self.atoms()
		self.bonds()
		self.substructures()

		# close the file
		super(Writer, self).write()
