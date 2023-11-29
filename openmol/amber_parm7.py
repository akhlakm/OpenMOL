#!/usr/bin/env python3

""" AMBER PARM7 and RESTART file reader.

	This file is a part of OpenMOL python module.
	License GPLv3.0 Copyright (c) 2023 Akhlak Mahmood """

from openmol import core

# AMBER PARM7 pointers list
# See http://ambermd.org/formats.html

pointers = [
	'NATOM', 	'NTYPES', 'NBONH',  'MBONA',  'NTHETH', 'MTHETA',
	'NPHIH',    'MPHIA',  'NHPARM', 'NPARM',  'NNB',    'NRES',
	'NBONA',    'NTHETA', 'NPHIA',  'NUMBND', 'NUMANG', 'NPTRA',
	'NATYP',    'NPHB',   'IFPERT', 'NBPER',  'NGPER',  'NDPER',
	'MBPER',    'MGPER',  'MDPER',  'IFBOX',  'NMXRS',  'IFCAP',
	'NUMEXTRA', 'NCOPY'
]

def initialize():
	""" Initialize an openmol object with Amber
		specific items. """

	MOL = core.initialize()
	MOL['source_format'] = "AMBER PARM7"

	MOL['parm_version_string'] = None
	MOL['atom_no_excluded'] = []

	MOL['parm7_lj_acoeff'] = []
	MOL['parm7_lj_bcoeff'] = []
	MOL['parm7_lj_index'] = []

	MOL['parm7_lj_epsilon'] = []
	MOL['parm7_lj_sigma'] = []

	for i in pointers:
		MOL['PARM_%s' %i] = 0

	return MOL


def build(MOL):
	""" Go through the openmol object and see if PARM7 
		specific items are properly calculated. If not,
		attemt to calculate them. """

	# add/update with default parm7 items
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

	MOL['no_bond_types'] = len(MOL['FF_bond_k'])
	MOL['no_angle_types'] = len(MOL['FF_angle_k'])
	MOL['no_dihed_types'] = len(MOL['FF_dihed_k'])

	# Build the unique atom types
	MOL['unique_atom_types'] = list(set(MOL['atom_type']))

	# If we have A, B coeffs, build epsilon, sigma of parm7
	if len(MOL['parm7_lj_acoeff']) and \
			len(MOL['parm7_lj_sigma']) != len(MOL['unique_atom_types']):

		if not MOL.get('parm7_lj_index', False):
			print('-- PARM7 Build Error: non bonded parm7 indices not found.')

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

	if len(MOL['parm7_lj_epsilon']) != MOL['PARM_NTYPES'] or \
			len(MOL['parm7_lj_sigma']) != MOL['PARM_NTYPES']:
		print('-- PARM7 Build Error: fail to build parm7 lj coeffs, length mismatch.')

	# build indices of types
	if len(MOL['atom_type_index']) != MOL['no_atoms']:
		MOL['atom_type_index'] = []
		for i in range(MOL['no_atoms']):
			atom_type = MOL['atom_type'][i]
			MOL['atom_type_index'].append(MOL['unique_atom_types'].index(atom_type))

	if len(MOL['atom_type_index']) != MOL['no_atoms']:
		print('-- PARM7 Build Error: fail to build atom type indices, length mismatch.')

	# assuming we have parm7 epsilon sigma lists built, build the lj params of each type
	if len(MOL['FF_lj_epsilon']) != MOL['no_atom_types'] or \
			len(MOL['FF_lj_sigma']) != MOL['no_atom_types']:
		
		MOL['FF_lj_epsilon'] = []
		MOL['FF_lj_sigma'] = []
		for i in range(MOL['no_atom_types']):
			# get the first atom of this type
			aix = MOL['atom_type_index'].index(i)
			# get parm7 pair ff index of that atom
			pfx = MOL['pair_ff_index'][aix]
			MOL['FF_lj_epsilon'].append(MOL['parm7_lj_epsilon'][pfx])
			MOL['FF_lj_sigma'].append(MOL['parm7_lj_sigma'][pfx])

	if len(MOL['FF_lj_epsilon']) != MOL['no_atom_types'] or \
			len(MOL['FF_lj_sigma']) != MOL['no_atom_types']:
		print('-- PARM7 Build Error: fail to build pair coeffs, length mismatch.')

	MOL['_parm7_built'] = True
	return core.AttrDict(MOL)


def process_last_section(MOL, section, lines, sformat):
	""" Parse and process the last read section of PARM7 file """

	items = []
	for line in lines:
		if sformat == "10I8":
			n = 8
			for i in range(0, len(line), n):
				item = line[i:i+n].strip()
				if len(item) > 0:
					items.append(item)
		else:
			items += line.strip().split()

	if not section:
		# no previous section
		return True

	if section == 'TITLE':
		MOL['title'] = lines[0]

	elif section == 'POINTERS':
		if len(items) > len(pointers):
			print('\n-- Warning: unknown PRMTOP pointer found. Ignoring ...', end=' ')

		for i, v in enumerate(items):
			if i < len(pointers):
				MOL['PARM_%s' %pointers[i]] = int(v)

		MOL['no_atoms'] = MOL['PARM_NATOM']
		MOL['no_bonds'] = MOL['PARM_NBONA'] + MOL['PARM_NBONH']
		MOL['no_angles'] = MOL['PARM_NTHETH'] + MOL['PARM_MTHETA']
		MOL['no_diheds'] = MOL['PARM_NPHIH'] + MOL['PARM_MPHIA']
		MOL['no_residues'] = MOL['PARM_NRES']
		MOL['no_atom_types'] = MOL['PARM_NATYP']

	elif section == 'ATOM_NAME':
		if len(items) != MOL['no_atoms']:
			print('\n-- Error: no_atoms and ATOM_NAME section mismatch')
			print('-- no of atoms: %d, record found: %d' %(MOL['no_atoms'], len(items)))
			print('-- this may happen if the system is huge.')
			print('-- reparsing section lines for atom names with 20A4 formatting ...')
			items = []
			for l in lines:
				n = 4
				items += [l[i:i+n].strip() for i in range(0, len(l), n)]

			if len(items) != MOL['no_atoms']:
				print('\n-- Error: no_atoms and ATOM_NAME section mismatch after reparsing')
				print('-- no of atoms: %d, record found: %d' %(MOL['no_atoms'], len(items)))
				print('-- nothing else to do. :( ')
				return False

		for name in items:
			MOL['atom_name'].append(name)

	elif section == 'CHARGE':
		if len(items) != MOL['no_atoms']:
			print('-- Error: no_atoms and CHARGE section mismatch')
			return False

		for q in items:
			# in electionic units, see http://ambermd.org/formats.html
			MOL['atom_q'].append(float(q)/18.2223)

	elif section == 'ATOMIC_NUMBER':
		if len(items) != MOL['no_atoms']:
			print('-- Error: no_atoms and ATOMIC_NUMBER section mismatch')
			return False

		for A in items:
			MOL['atom_atomic_no'].append(int(A))

	elif section == 'MASS':
		if len(items) != MOL['no_atoms']:
			print('-- Error: no_atoms and MASS section mismatch')
			return False

		for m in items:
			MOL['atom_mass'].append(float(m))

	# lj_ff_index for each atom
	elif section == 'ATOM_TYPE_INDEX':
		if len(items) != MOL['no_atoms']:
			print('-- Error: no_atoms and ATOM_TYPE_INDEX section mismatch')
			return False

		for t in items:
			# 1 based indexing, subtract 1
			MOL['pair_ff_index'].append(int(t) - 1)

	elif section == 'NUMBER_EXCLUDED_ATOMS':
		if len(items) != MOL['no_atoms']:
			print('-- Error: no_atoms and NUMBER_EXCLUDED_ATOMS section mismatch')
			return False

		for ex in items:
			MOL['atom_no_excluded'].append(int(ex))

	# lj parm index, needed to find epsilon, sigma
	elif section == 'NONBONDED_PARM_INDEX':
		if len(items) != MOL['PARM_NTYPES']**2:
			print('-- Error: PARM_NTYPES and NONBONDED_PARM_INDEX section mismatch')
			return False

		for ix in items:
			MOL['parm7_lj_index'].append(int(ix))

	elif section == 'RESIDUE_LABEL':
		if len(items) != MOL['no_residues']:
			print('-- Error: no_residues and RESIDUE_LABEL section mismatch')
			return False

		for i in items:
			MOL['residue_name'].append(i)

	elif section == 'RESIDUE_POINTER':
		if len(items) != MOL['no_residues']:
			print('-- Error: no_residues and RESIDUE_POINTER section mismatch')
			return False

		for i in items:
			MOL['residue_start'].append(int(i) - 1)

	elif section in ['BONDS_INC_HYDROGEN', 'BONDS_WITHOUT_HYDROGEN']:
		for i in range(0, len(items), 3):
			index0 = int(abs(int(items[i]))/3)
			index1 = int(abs(int(items[i+1]))/3)
			index2 = int(items[i+2]) - 1
			# we do not distinguish between H or other atoms for now
			# store as regular bond info
			MOL['bond_from'].append(index0)
			MOL['bond_to'].append(index1)
			MOL['bond_ff_index'].append(index2)

	elif section in ['ANGLES_INC_HYDROGEN', 'ANGLES_WITHOUT_HYDROGEN']:
		for i in range(0, len(items), 4):
			index0 = int(abs(int(items[i]))/3)
			index1 = int(abs(int(items[i+1]))/3)
			index2 = int(abs(int(items[i+2]))/3)
			index3 = int(items[i+3]) - 1
			# we do not distinguish between H or other atoms for now
			# store as regular angle info
			MOL['angle_a'].append(index0)
			MOL['angle_b'].append(index1)
			MOL['angle_c'].append(index2)
			MOL['angle_ff_index'].append(index3)

	elif section in ['DIHEDRALS_INC_HYDROGEN', 'DIHEDRALS_WITHOUT_HYDROGEN']:
		for i in range(0, len(items), 5):
			index0 = int(abs(int(items[i]))/3)
			index1 = int(abs(int(items[i+1]))/3)
			index2 = int(abs(int(items[i+2]))/3)
			index3 = int(abs(int(items[i+3]))/3)
			index4 = int(items[i+4]) - 1
			# we do not distinguish between H or other atoms for now
			# store as regular dihedral info
			MOL['dihed_a'].append(index0)
			MOL['dihed_b'].append(index1)
			MOL['dihed_c'].append(index2)
			MOL['dihed_d'].append(index3)
			MOL['dihed_ff_index'].append(index4)

	elif section == 'AMBER_ATOM_TYPE':
		if len(items) != MOL['no_atoms']:
			print('-- Error: no_atoms and AMBER_ATOM_TYPE section mismatch')
			return False

		for t in items:
			MOL['atom_type'].append(t)
			if t not in MOL['unique_atom_types']:
				MOL['unique_atom_types'].append(t)

	elif section == 'BOND_FORCE_CONSTANT':
		for i in items:
			MOL['FF_bond_k'].append(float(i))

	elif section == 'BOND_EQUIL_VALUE':
		for i in items:
			MOL['FF_bond_eq'].append(float(i))

	elif section == 'ANGLE_FORCE_CONSTANT':
		for i in items:
			MOL['FF_angle_k'].append(float(i))

	elif section == 'ANGLE_EQUIL_VALUE':
		for i in items:
			MOL['FF_angle_eq'].append(float(i))

	elif section == 'DIHEDRAL_FORCE_CONSTANT':
		for i in items:
			MOL['FF_dihed_k'].append(float(i))

	elif section == 'DIHEDRAL_PERIODICITY':
		for i in items:
			MOL['FF_dihed_periodicity'].append(float(i))

	elif section == 'DIHEDRAL_PHASE':
		for i in items:
			MOL['FF_dihed_phase'].append(float(i))

	elif section == 'LENNARD_JONES_ACOEF':
		for i in items:
			MOL['parm7_lj_acoeff'].append(float(i))

	elif section == 'LENNARD_JONES_BCOEF':
		for i in items:
			MOL['parm7_lj_bcoeff'].append(float(i))

	else:
		print('IGNORED')
		return True

	print('OK')
	return True

def read_prmtop(prmtop):
	global section_lines

	MOL = initialize()

	line_no = 0					# global line number
	section_line_no = 0			# line number within current section
	section = None
	section_format = None
	section_lines = []			# lines of current section

	for line in open(prmtop, 'r'):
		line_no += 1
		section_line_no += 1

		if len(line) == 0:
			continue

		if line.startswith('%VERSION'):
			parts = line.split()
			if len(parts) < 2:
				print('-- Error: Invalid PRMTOP [line %d]:\n%s' %(line_no, line))
				return None
			else:
				MOL['parm_version_string'] = ' '.join(parts[1:])

		# new section
		elif line.startswith('%FLAG'):
			parts = line.split()
			if len(parts) < 2:
				print('-- Error: Invalid PRMTOP [line %d]:\n%s' %(line_no, line))
				return None
			else:
				# new section found, first process the previous section if any
				if not process_last_section(MOL, section, section_lines, section_format):
					return False

				section = parts[1]
				print('Reading %s ...' %section, end=' ')
				section_line_no = 0
				section_lines = []

		elif line.startswith('%FORMAT'):
			parts = line.strip().split('(')
			if len(parts) < 2:
				print('-- Error: Invalid PRMTOP [line %d]:\n%s' %(line_no, line))
				return None
			else:
				section_format = parts[1][:-1]

		else:
			section_lines.append(line)

	# process the final section
	if not process_last_section(MOL, section, section_lines, section_format):
		return False

	print('Reading Done')
	return MOL

def read_rst7(MOL, rst_file):
	line_no = 0
	no_atoms = 0
	items = []
	last_line = None

	for line in open(rst_file, 'r'):
		line_no += 1
		line = line.strip()

		if len(line) == 0:
			continue

		last_line = line

		if line_no == 1:
			rst_title = line
			continue
		elif line_no == 2:
			parts = line.split()
			no_atoms = int(parts[0])
			if len(parts) > 1:
				# optional time string
				MOL['time'] = float(parts[1])
			if len(parts) > 2:
				# optional temperature string
				MOL['temp'] = float(parts[2])

			if no_atoms != MOL['no_atoms']:
				print('-- Error: RST7 no_atoms mismatch with the PARM7 file.')
				return False
		else:
			# rest is all the atomic coordinates in 3D
			items += line.split()

	if len(items) < no_atoms * 3:
		print('-- Error: RST7 no_atoms, coordinate items mismatch.')
		print('coordinates: ', len(items), 'no of atoms: ', no_atoms)
		return False

	print('Reading coordinates ...', end=' ')

	for i in range(0, no_atoms*3, 3):
		MOL['atom_x'].append(float(items[i]))
		MOL['atom_y'].append(float(items[i+1]))
		MOL['atom_z'].append(float(items[i+2]))

	print('OK')

	# read optional velocities
	if len(items) >= no_atoms*6:
		print('Reading velocities ...', end=' ')

		for i in range(no_atoms*3, no_atoms*6, 3):
			MOL['atom_vx'].append(float(items[i]))
			MOL['atom_vy'].append(float(items[i+1]))
			MOL['atom_vz'].append(float(items[i+2]))

		print('OK')

	box_size_conditions = [
		len(items) == no_atoms*3 + 3,	# x,y,z coordinates + a,b,c
		len(items) == no_atoms*6 + 3,	# x,y,z coordinates and velocities + a,b,c
		len(items) == no_atoms*3 + 6,	# x,y,z coordinates + a,b,c + angles
		len(items) == no_atoms*6 + 6,	# x,y,z coordinates and velocities + a,b,c + angles
	]

	if any(box_size_conditions):
		print('Reading box dimension ...', end=' ')
		parts = last_line.split()
		MOL['box_x'] = float(parts[0])
		MOL['box_y'] = float(parts[1])
		MOL['box_z'] = float(parts[2])

		print('OK')
	else:
		print('-- No PBC box information found.')

	box_angle_conditions = [
		len(items) == no_atoms*3 + 6,
		len(items) == no_atoms*6 + 6,
	]

	if any(box_angle_conditions):
		print('Reading box angles ...', end=' ')
		parts = last_line.split()
		MOL['box_alpha'] = float(parts[3])
		MOL['box_beta'] = float(parts[4])
		MOL['box_gamma'] = float(parts[5])

		print('OK')

	print('Reading Done')
	return MOL

def read(prmtop, rst7):
	print(f"\nReading {prmtop} and {rst7}")
	MOL = read_prmtop(prmtop)
	if not MOL:
		return False
	
	MOL = read_rst7(MOL, rst7)
	if not core.check(MOL):
		return False

	return MOL
