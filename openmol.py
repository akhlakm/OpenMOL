import json

def initialize():
	MOL = {}

	MOL['title'] = None
	MOL['description'] = ""
	MOL['source_format'] = None
	MOL['type'] = None
	MOL['charge_type'] = None

	MOL['no_atoms'] = 0
	MOL['no_bonds'] = 0
	MOL['no_atom_types'] = 0
	MOL['no_residues'] = 0
	MOL['no_angles'] = 0
	MOL['no_dihedrals'] = 0
	MOL['unique_atom_types'] = []

	MOL['atom_name'] = []
	MOL['atom_x'] = []
	MOL['atom_y'] = []
	MOL['atom_z'] = []
	MOL['atom_vx'] = []
	MOL['atom_vy'] = []
	MOL['atom_vz'] = []
	MOL['atom_q'] = []
	MOL['atom_type'] = []
	MOL['atom_resname'] = []
	MOL['atom_resid'] = []
	MOL['atom_mass'] = []
	MOL['atom_atomic_no'] = []

	MOL['bond_from'] = []
	MOL['bond_to'] = []
	MOL['bond_type'] = []

	MOL['angle_a'] = []
	MOL['angle_b'] = []
	MOL['angle_c'] = []

	MOL['dihed_a'] = []
	MOL['dihed_b'] = []
	MOL['dihed_c'] = []
	MOL['dihed_d'] = []

	MOL['residue_name'] = []
	MOL['residue_start'] = []
	MOL['residue_end'] = []
	MOL['residue_type'] = []

	MOL['box_x'] = 0.0
	MOL['box_y'] = 0.0
	MOL['box_z'] = 0.0
	MOL['box_alpha'] = 0.0
	MOL['box_beta'] = 0.0
	MOL['box_gamma'] = 0.0

	return MOL


def write_json(MOL, json_file, compress=False):
	with open(json_file, 'w+') as fp:
		if compress:
			json.dump(MOL,fp, indent=None, separators=(',', ':'))
		else:
			json.dump(MOL, fp, indent=4)

	print('%s written' %json_file)


def load_json(json_file):
	with open(json_file, 'r') as fp:
		MOL = json.load(fp)

	MOL['json_file'] = json_file
	print('%s loaded' %json_file)

	return MOL


def check_atoms_ok(MOL):
	conditions_fail = [
		not MOL['no_atoms'],
		len(MOL['atom_q']) and MOL['no_atoms'] != len(MOL['atom_y']),
	]

	conditions_ok = [
		MOL['no_atoms'] == len(MOL['atom_name']),
		MOL['no_atoms'] == len(MOL['atom_x']),
		MOL['no_atoms'] == len(MOL['atom_y']),
		MOL['no_atoms'] == len(MOL['atom_z']),
		MOL['no_atoms'] == len(MOL['atom_type']),
	]

	if any(conditions_fail) or not all(conditions_ok):
		print('Error: Summary and atom attribute list mismatch.')
		return False

	return True

def check_bonds_ok(MOL):
	conditions_fail = [
		MOL['no_bonds'] == None and len(MOL['bond_from']) > 0,
	]

	conditions_ok = [
		MOL['no_bonds'] == len(MOL['bond_from']),
		MOL['no_bonds'] == len(MOL['bond_to']),
	]

	if any(conditions_fail) or not all(conditions_ok):
		print('Error: Summary and bond attribute list mismatch.')
		return False

	return True

def check_residues_ok(MOL):
	conditions_fail = [
		MOL['no_residues'] == None and len(MOL['residue_name']) > 0,
	]

	conditions_ok = [
		MOL['no_residues'] == len(MOL['residue_name']),
		MOL['no_residues'] == len(MOL['residue_start']),
	]

	if any(conditions_fail) or not all(conditions_ok):
		print('Error: Summary and residue attribute list mismatch.')
		return False

	return True


def check(MOL):
	return check_atoms_ok(MOL) and check_bonds_ok(MOL) and check_residues_ok(MOL)

def update_summary(MOL, overwrite=False):
	if overwrite or MOL['no_atoms'] is None:
		MOL['no_atoms'] = len(MOL['atom_x'])

	if overwrite or MOL['no_bonds'] is None:
		MOL['no_bonds'] = len(MOL['bond_from'])

	if overwrite or MOL['no_residues'] is None:
		MOL['no_residues'] = len(MOL['residue_name'])

	if overwrite or MOL['no_atom_types'] is None:
		MOL['unique_atom_types'] = list(set(MOL['atom_type']))
		MOL['no_atom_types'] = len(MOL['unique_atom_types'])
