#!/usr/bin/env python3

""" Define main OpenMOL object definition and helper functions.

    This file is a part of OpenMOL python module.
    License GPLv3.0 Copyright (c) 2023 Akhlak Mahmood """

import json
import pandas as pd

from .utils import AttrDict

def initialize():
    """ Generate empty OpenMOL dictionary object with
        the default properties.
        Note: openmol uses 0 based indexing. """

    MOL = {}

    MOL['title'] = None
    MOL['description'] = ""
    MOL['source_format'] = None
    MOL['type'] = None
    MOL['charge_type'] = None

    MOL['no_atoms'] = 0
    MOL['no_bonds'] = 0
    MOL['no_angles'] = 0
    MOL['no_diheds'] = 0
    MOL['no_improper'] = 0
    MOL['no_residues'] = 0

    MOL['no_atom_types'] = 0
    MOL['no_bond_types'] = 0
    MOL['no_angle_types'] = 0
    MOL['no_dihed_types'] = 0
    MOL['no_improper_types'] = 0

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
    MOL['atom_type_index'] = []
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

    MOL['improper_a'] = []
    MOL['improper_b'] = []
    MOL['improper_c'] = []
    MOL['improper_d'] = []

    MOL['residue_name'] = []
    MOL['residue_start'] = []
    MOL['residue_end'] = []
    MOL['residue_type'] = []

    MOL['box_x'] = 0.0
    MOL['box_y'] = 0.0
    MOL['box_z'] = 0.0
    MOL['box_alpha'] = 90.0
    MOL['box_beta'] = 90.0
    MOL['box_gamma'] = 90.0

    # FF parameters
    MOL['FF_lj_epsilon'] = []			# epsilon of each atom type
    MOL['FF_lj_sigma'] = []				# sigma of each atom type
    MOL['pair_ff_index'] = []			# index connecting atom_type_index

    MOL['FF_bond_k'] = []
    MOL['FF_bond_eq'] = []
    MOL['bond_ff_index'] = []

    MOL['FF_angle_k'] = []
    MOL['FF_angle_eq'] = []
    MOL['angle_ff_index'] = []

    MOL['FF_dihed_k'] = []
    MOL['FF_dihed_phase'] = []			# in radians
    MOL['FF_dihed_periodicity'] = []
    MOL['dihed_ff_index'] = []

    MOL['FF_improper_k'] = []
    MOL['FF_improper_phase'] = []		# in radians
    MOL['FF_improper_periodicity'] = []
    MOL['improper_ff_index'] = []

    return AttrDict(MOL)


def write_json(MOL, json_file, compress=False):
    """ Write the openmol object as openmol JSON file
        Optional compress argument can be used to save
        without any indentation. """

    with open(json_file, 'w+') as fp:
        if compress:
            json.dump(MOL,fp, indent=None, separators=(',', ':'))
        else:
            json.dump(MOL, fp, indent=4)

    print('Write OK: %s' %json_file)


def load_json(json_file):
    """ Load a openmol type JSON file and return the openmol object """
    with open(json_file, 'r') as fp:
        MOL = json.load(fp)

    MOL['source_json'] = json_file
    print('Load OK: %s' %json_file)

    return AttrDict(MOL)


def check_atoms_ok(MOL):
    conditions_fail = [
        not MOL['no_atoms'],
        len(MOL['atom_q']) and MOL['no_atoms'] != len(MOL['atom_q']),
    ]

    conditions_ok = [
        MOL['no_atoms'] == len(MOL['atom_name']),
        MOL['no_atoms'] == len(MOL['atom_x']),
        MOL['no_atoms'] == len(MOL['atom_y']),
        MOL['no_atoms'] == len(MOL['atom_z']),
        MOL['no_atoms'] == len(MOL['atom_type']),
    ]

    if any(conditions_fail) or not all(conditions_ok):
        print('-- Error: Summary and atom attribute list mismatch.')
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
        print('-- Error: Summary and bond attribute list mismatch.')
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
        print('-- Error: Summary and residue attribute list mismatch.')
        return False

    return True


def check(MOL):
    return check_atoms_ok(MOL) and check_bonds_ok(MOL) and check_residues_ok(MOL)

def update_summary(MOL, overwrite=False):
    """ Attempts to update the atom, type, bond and residue counts
        from the length of their list. Should be called if
        they have been changed manually. """

    if overwrite or MOL['no_atoms'] is None:
        MOL['no_atoms'] = len(MOL['atom_x'])

    if overwrite or MOL['no_bonds'] is None:
        MOL['no_bonds'] = len(MOL['bond_from'])

    if overwrite or MOL['no_residues'] is None:
        MOL['no_residues'] = len(MOL['residue_name'])

    if overwrite or MOL['no_atom_types'] is None:
        MOL['unique_atom_types'] = list(set(MOL['atom_type']))
        MOL['no_atom_types'] = len(MOL['unique_atom_types'])

    return MOL


class Reader(object):
    """ Base file reader interface to implement in different
        Reader classes. """

    def __init__(self):
        self.Mol = initialize()
        self.in_file = None
        self.lines = []
        self.section = None
        self.section_start = 0
        self.section_lines = []
        self.section_format = None


    def read_file(self, input_file : str = None):
        if input_file is not None:
            self.in_file = input_file
        self.lines = open(self.in_file).readlines()
        print('Reading:', input_file)
        self._process_lines()
        print('Read OK:', input_file)
        return self.Mol


    def as_df(self, prefix = 'atom_', fields = []):
        """ Return the Mol info section as a Pandas dataframe. """
        if fields:
            return pd.DataFrame({
                k : v for k, v in self.Mol.items() if k in fields
            })
        else:
            return pd.DataFrame({
                k : v for k, v in self.Mol.items() if k.startswith(prefix)
            })


    def _process_lines(self):
        for line_no, line in enumerate(self.lines):
            line = line.strip()
            try:
                next_line = self.lines[line_no+1]
                next_line = next_line.strip()
            except:
                next_line = None
            if len(line) == 0:
                continue
            self._identify_section(line_no+1, line, next_line)
            self.section_lines.append(line)

        self._process_last_section(self.section, self.section_lines,
                                   self.section_format)


    def _identify_section(self, line_no : int, line : str, next_line : str) -> None:
        raise NotImplementedError()


    def _new_section(self, section : str, line_no : int, section_format = None):
        """ Declare the start of a new section. """
        if section == self.section:
            return
        if self.section is not None:
            self._process_last_section(self.section, self.section_lines,
                                    self.section_format)
        self.section = section
        self.section_start = line_no
        self.section_format = section_format
        self.section_lines = []
        print('- Section: %s ...' %section, end=' ')


    def _process_last_section(self, section : str, lines : list, sformat : str):
        """ Parse and process the last read section of file """
        print('IGNORED')


    def _str_to_type(self, string : str, dtype : callable, line, line_no):
        errstr  = f"failed to parse {string} as {dtype}, "
        errstr += f"line {line_no}: {line}"
        try:
            value = dtype(string)
        except TypeError:
            raise TypeError(errstr)
        return value


class Writer(object):
    """ Base file writer interface to implement in different
        Writer classes. """

    def __init__(self, MOL, out_file):
        self.out_file = out_file
        self.fp = open(out_file, 'w+')
        self.MOL = MOL

    def title(self):
        self.fp.write("%s (by OpenMOL)\n\n" %self.MOL['title'])

    def close(self):
        self.fp.close()
        print('Write OK: %s' %self.out_file)


    def write(self):
        self.close()

    def save(self):
        self.write()

