#!/usr/bin/env python3

""" LAMMPS dat file reader and writer for atom_style 'full'.

    This file is a part of OpenMOL python module.
    License GPLv3.0 Copyright (c) 2023 Akhlak Mahmood """

from openmol import core

class Reader(core.Reader):
    def __init__(self):
        super(Reader, self).__init__()
        self.Mol = core.AttrDict()
        self.Mol['source_format'] = "PSF"
        self.Mol['_psf_built'] = False
        self.Mol['atom_q'] = []
        self.Mol['atom_name'] = []
        self.Mol['atom_type'] = []
        self.Mol['atom_mass'] = []
        self.Mol['atom_resid'] = []
        self.Mol['atom_resname'] = []
        self.Mol['atom_molecule'] = []

    def read(self, psf_file : str):
        super(Reader, self).read_file(psf_file)
        print("-- Warning: only atoms section is implemented")
        self._process_lines()

    def _parse_str_as_type(self, string : str, dtype : callable, line, i):
        errstr  = f"-- Read Error: failed to parse {string} as {dtype}, "
        errstr += f"line {i}: {line}"
        try:
            value = dtype(string)
        except TypeError:
            raise TypeError(errstr)
        return value

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
            if line.startswith("REMARK "):
                continue
            
            words = line.split()

            if len(words) < 2:
                # possibly comment or title
                continue

            if words[1] == "!NATOM" and section is None:
                # atom list
                section = 'atom'
                continue

            if words[1] == "!BOND" and section is None:
                # bond list
                section = 'bond'
                continue

            elif section == 'atom':
                if words[1][0] == "!":
                    print("-- WARN: Ignore section at line %d: %s" %(i+1, line))
                    section = None
                    continue

                assert len(words) >= 8, \
                    "Invalid atom info, line %d: %s" %(i+1, line)

                mol_name = words[1]
                res_id = self._parse_str_as_type(words[2], int, line, i)
                res_name = words[3]
                at_name = words[4] # atom name
                at_type = words[5] # type name
                at_q = self._parse_str_as_type(words[6], float, line, i)
                at_mass = self._parse_str_as_type(words[7], float, line, i)

                self.Mol.atom_q.append(at_q)
                self.Mol.atom_name.append(at_name)
                self.Mol.atom_type.append(at_type)
                self.Mol.atom_mass.append(at_mass)
                self.Mol.atom_resid.append(res_id - 1)
                self.Mol.atom_resname.append(res_name)
                self.Mol.atom_molecule.append(mol_name)

            else:
                pass

        print("Read OK")

