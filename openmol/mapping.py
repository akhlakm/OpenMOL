import os
import json

import openmol
import openmol.tripos_mol2 as mol2


class FFMap:
    """ Create a map between two force fields using atom names as
    the common identifier. """
    def __init__(self, ff1_name : str, ff2_name : str) -> None:
        self.ff1_name = ff1_name
        self.ff2_name = ff2_name

        self.ff1_type2name = {}
        self.ff1_name2type = {}
        self.ff1_name2charge = {}

        self.ff2_type2name = {}
        self.ff2_name2type = {}
        self.ff2_name2charge = {}


    def _read(self, filepath):
        if filepath.endswith(".json"):
            return openmol.load_json(filepath)
        elif filepath.endswith(".mol2"):
            molecule = mol2.read(filepath)
            molecule = mol2.build(molecule)
            return molecule
        else:
            raise ValueError("Unsupported file format", filepath)

    def save_mapping(self, filename="ff.map.json"):
        with open(filename, "w+") as fp:
            json.dump(self.__dict__, fp)

    def load_mapping(self, filename="ff.map.json"):
        with open(filename, "r") as fp:
            d = json.load(fp)
            self.__dict__.update(d)


    def load_ff1_file(self, filepath):
        """ Load the type one system file (json/mol2). """
        mol = self._read(filepath)

        # Build name list
        for i, atom_name in enumerate(mol.atom_name):
            atom_type = mol.atom_type[i]
            atom_charge = mol.atom_q[i]

            if atom_type not in self.ff1_type2name:
                self.ff1_type2name[atom_type] = atom_name

            if atom_name not in self.ff1_name2type:
                self.ff1_name2type[atom_name] = atom_type

            if atom_name not in self.ff1_name2charge:
                self.ff1_name2charge[atom_name] = atom_charge


    def load_ff2_file(self, filepath):
        """ Load the type two system file (json/mol2). """
        mol = self._read(filepath)

        # Build name list
        for i, atom_name in enumerate(mol.atom_name):
            atom_type = mol.atom_type[i]
            atom_charge = mol.atom_q[i]

            if atom_type not in self.ff2_type2name:
                self.ff2_type2name[atom_type] = atom_name

            if atom_name not in self.ff2_name2type:
                self.ff2_name2type[atom_name] = atom_type

            if atom_name not in self.ff2_name2charge:
                self.ff2_name2charge[atom_name] = atom_charge
