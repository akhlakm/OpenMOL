#!/usr/bin/env python
"""
Create a library of atomic charges and 
use the library to estimate charges of new molecules.
This allows extremely fast estimation of the charges of systems
with repeating units.

Max standard error of the estimated atomic charges = 0.05 e.

"""

__author__  = "Akhlak Mahmood"
__version__ = 23.1128

import os
import json
import numpy as np
from openmol import tripos_mol2 as mol2

def find_connections(mol, atom_index) -> list[tuple]:
    """ Find the indices of all atoms bonded to the specified atom
        and their types.
    """
    bonded_atoms = []

    for jx, bond_type in enumerate(mol.bond_type):
        if mol.bond_from[jx] == atom_index:
            ot = mol.bond_to[jx]
            bonded_atoms.append((ot, bond_type))
        elif mol.bond_to[jx] == atom_index:
            ot = mol.bond_from[jx]
            bonded_atoms.append((ot, bond_type))

    return bonded_atoms


def find_neighbors(mol, atom_index) -> list:
    """ Find the indices and bond types of the atoms
        connected to all the bonded atoms.
    """
    neighbor_list = []
    bonded_atoms = find_connections(mol, atom_index)
    for ix, it in bonded_atoms:
        neighbors = find_connections(mol, ix)
        for ni, nt in neighbors:
            neighbor_list.append((it, ix, nt, ni))
    return neighbor_list


def unique_identifier(mol, atom_index) -> str:
    """
    Calculate a unique identifier of an atom based on it's connected atoms
    and atoms connected to those atoms (neighbors).
    Use atom types and a comma separated string for the neighbors.
    Output format is ATOM(BONDTYPE-BONDEDATOM,BONDTYPE-NEIGHBOR1,),
    """
    nlist = find_neighbors(mol, atom_index)
    bonded = {}
    bondtypes = {}

    # dict keys must be index, if type is used, two different atoms of
    # the same types will be merged together.
    for bt, bi, nt, ni in nlist:
        if bi not in bonded:
            bonded[bi] = []
        if ni != atom_index:
            bonded[bi].append((nt, ni))
        bondtypes[bi] = bt

    for bi, v in bonded.items():
        assert len(v) <= 4, "Too many boned atoms"
        bat = mol.atom_type[bi]
        bty = bondtypes[bi]
        bond_str = f"{bty}-{bat}"
        neighbor_str = ",".join([
            nt +"-"+ mol.atom_type[ni] for nt, ni in sorted(v)
        ])
        bonded[bi] = f"{bond_str}:{neighbor_str}" if neighbor_str else bond_str

    unique = mol.atom_type[atom_index]
    for v in sorted(bonded.values()):
        unique += f"({v})"

    return unique


def unique_charge_dict(mol, charge_dict) -> dict:
    """
    Calculate a unique identifier of an atom based on it's connected atoms
    and atoms connected to those atoms (neighbors).
    Then make a dict of charges for those atoms.

    We could consider neighbors of the neighbors to uniquely identify the atoms,
    but it would make the generated library less general.
    """

    cc = { k : v for k, v in charge_dict.items() }

    for i in range(len(mol.atom_name)):
        ch = mol.atom_q[i]
        ids = unique_identifier(mol, i)
        if ids not in cc:
            cc[ids] = []
            print("New signature:", ids)
        cc[ids].append(ch)

    return cc
    



def create_lib_files(lib_prefix, charge_dict = None):
    """ Create the output lib and err files.
        If a charge_dict is given, use that, otherwise load the json file.
    """
    if charge_dict is None:
        # Load the charges from json file.
        outputjson = lib_prefix + ".json"
        with open(outputjson) as fp:
            charge_dict = json.load(fp)
    
    assert type(charge_dict) == dict 
    print("Total atom signatures:", len(charge_dict))

    outputlib = lib_prefix + ".lib"
    outputerr = lib_prefix + ".err"

    # Use the mean values for prediction.
    charge_lib = {}
    for k, v in charge_dict.items():
        charge_lib[k] = np.mean(v)

    # If std err is too high, the charges will not be reliable.
    charge_err = {}
    for k, v in charge_dict.items():
        charge_err[k] = np.std(v)

    # Save as a separate lib and err files.
    with open(outputlib, "w") as fp:
        json.dump({k : v for k, v in sorted(charge_lib.items())}, fp, indent=4)
    print("Saved:", outputlib)

    with open(outputerr, "w") as fp:
        json.dump({k : v for k, v in sorted(charge_err.items())}, fp, indent=4)
    print("Saved:", outputerr)


def create_charge_dir(libdir, lib_prefix):
    """ Read all the mol2 files in a directory.
        Create unique identifier of each atom.
        Create a map of atomic charges using the identifiers.
        Save to an output json file to load in the next run.
        Also save the mean values in a lib file for current use.

        There are normally slight variations (~ 0.05 e max) in the
        charges for each atoms. So we average them out.

    """
    outputjson = lib_prefix + ".json"
    
    if os.path.isfile(outputjson):
        with open(outputjson) as fp:
            charge_dict = json.load(fp)
    else:
        charge_dict = {}

    for f in os.listdir(libdir):
        path = os.path.join(libdir, f)
        if f.endswith(".mol2"):
            print("Reading", path)
            mol = mol2.read(path)
            unique_charge_dict(mol, charge_dict)
            print("Total:", len(charge_dict))

    with open(outputjson, "w") as fp:
        json.dump({k: v for k, v in charge_dict.items()}, fp)

    print("Saved:", outputjson)
    create_lib_files(lib_prefix, charge_dict)



def create_charge_file(mol2file, lib_prefix):
    """ Read the mol2 file.
        Create unique identifier of each atom.
        Create a map of atomic charges using the identifiers.
        Save to an output json file to load in the next run.
        Also save the mean values in a lib file for current use.

        There are normally slight variations (~ 0.05 e max) in the
        charges for each atoms. So we average them out.

    """
    outputjson = lib_prefix + ".json"

    if os.path.isfile(outputjson):
        with open(outputjson) as fp:
            charge_dict = json.load(fp)
    else:
        charge_dict = {}
    print("Total signatures loaded:", len(charge_dict))

    mol = mol2.read(mol2file)

    charge_dict = unique_charge_dict(mol, charge_dict)
    print("New total signatures:", len(charge_dict))

    with open(outputjson, "w") as fp:
        json.dump({k: v for k, v in charge_dict.items()}, fp)

    print("Saved:", outputjson)
    create_lib_files(lib_prefix, charge_dict)


def calculate_charges(mol2file, outfile = None, lib_prefix="charges"):
    """ Calculate charges based on atomic connectivity using the
        charges.lib file. Output as a mol2 file.
        The lib file must exist in the working directory.

        Raises ValueError if the atom signature is not in the library.
    """
    # Overwrite the mol2file
    outfile = outfile if outfile else mol2file+".lib.mol2"
    libfile = lib_prefix + ".lib"

    # Load the charge library
    with open(libfile) as fp:
        charge_dict = json.load(fp)

    print("Reading", mol2file)
    mol = mol2.read(mol2file)

    for i in range(len(mol.atom_name)):
        atom_id = unique_identifier(mol, i)

        # check if we have it in the library
        if atom_id in charge_dict:
            mol.atom_q[i] = charge_dict[atom_id]
        else:
            raise ValueError(
                f"Error: Atom {i+1} with signature {atom_id} "
                f"is not in charge lib.")

    # Save the output file.
    wr = mol2.Writer(mol2.build(mol), outfile)
    wr.write()

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-d', '--dir', default=None,
        help="Directory containing mol2 files to build/update charge library.")

    parser.add_argument(
        '-a', '--add', default=None,
        help="MOL2 file to add to the charge library.")

    parser.add_argument(
        '-f', '--mol2', default=None,
        help="MOL2 file to calculate charges using the charge library.")

    parser.add_argument(
        '-o', '--outmol2', default=None,
        help="Output MOL2 file to store the charges.")

    parser.add_argument(
        '-b', '--build', default=False, action='store_true',
        help="Specify to build the charge.lib and charge.err files.")

    parser.add_argument(
        '-l', '--lib_prefix', default="charges",
        help="Path prefix to charges.json and charges.lib files.")

    args = parser.parse_args()

    if args.build:
        create_lib_files(args.lib_prefix)

    elif args.dir:
        # Update the charges.lib by parsing exisiting mol2 files.
        create_charge_dir(args.dir, args.lib_prefix)

    elif args.add:
        # Calculate charges for a mol2 file using the charges.lib file.
        create_charge_file(args.add, args.lib_prefix)

    elif args.mol2:
        # Calculate charges for a mol2 file using the charges.lib file.
        calculate_charges(args.mol2, args.outmol2, args.lib_prefix)
