import numpy as np
from openmol.array import Molecule, MolReader, Info, Crds, Strc

class Mol2Reader(MolReader):
    def __init__(self):
        super().__init__()

        self.meta.comment = ""
        self.meta.atom_status_bit = []
        self.meta.bond_from = []
        self.meta.bond_to = []
        self.meta.bond_type = []
        self.meta.bond_status_bit = []

    def _identify_section(self, line_no: int, line: str, next_line : str):
        words = line.split()

        if line.startswith('#'):
            self._new_section('COMMENT', line_no)

        elif line.startswith('@<'):
            parts = line.split('>')
            section = parts[1]
            self._new_section(section, line_no)

    def _process_last_section(self, section: str, lines: list, sformat: str):
        if section == 'COMMENT':
            self.meta.comment += " ".join(line[7:].strip() for line in lines)
            print("OK")

        elif section == 'MOLECULE':
            for i, line in enumerate(lines):
                if i == 0: continue
                elif i == 1:
                    self.meta.title = line
                elif i == 2:
                    parts = line.split()
                    if len(parts) > 0:
                        self.meta.no_atoms = int(parts[0])

                    if len(parts) > 1:
                        self.meta.no_bonds = int(parts[1])

                    if len(parts) > 2:
                        self.meta.no_residues = int(parts[2])

                    if len(parts) > 3:
                        self.meta.no_features = int(parts[3])

                    if len(parts) > 4:
                        self.meta.no_sets = int(parts[4])
                elif i == 3:
                    self.meta.mol2_type = line
                elif i == 4:
                    self.meta.charge_type = line
                elif i == 5:
                    self.meta.mol2_status_bits = line
                elif i == 6:
                    self.meta.comment += line
            print('OK')
                    
        elif section == 'ATOM':
            for i, line in enumerate(lines):
                if i == 0: continue
                parts = line.split()

                if len(parts) < 6:
                    raise ValueError('Invalid MOL2 [line %d]:\n%s' %(i, line))

                info, crds, strc = self.new_atom()

                info[Info.Name]  = parts[1]
                crds[Crds.X]     = float(parts[2])
                crds[Crds.Y]     = float(parts[3])
                crds[Crds.Z]     = float(parts[4])
                info[Info.Type]  = parts[5]
                strc[Strc.SegId] = 0

                # optional items
                if len(parts) > 6:
                    strc[Strc.ResId] = int(parts[6]) - 1

                if len(parts) > 7:
                    info[Info.ResName] = parts[7]

                if len(parts) > 8:
                    info[Info.Q] = round(float(parts[8]), 4)

                if len(parts) > 9:
                    self.meta.atom_status_bit.append(parts[9])

                self.add_atom(info, crds, strc)
            print('OK')

        elif section == 'BOND':
            for i, line in enumerate(lines):
                if i == 0: continue
                parts = line.split()

                if len(parts) < 4:
                    raise ValueError('Invalid MOL2 [line %d]:\n%s' %(i, line))

                at1 = int(parts[1]) - 1
                at2 = int(parts[2]) - 1
                bond_type = int(parts[3])

                if at1 >= len(self.strc):
                    raise ValueError("Atom index too high in bond")

                if at2 >= len(self.strc):
                    raise ValueError("Atom index too high in bond")

                for i in range(bond_type):
                    for bn in Strc.Bonds:
                        if self.strc[at1][bn] is np.nan:
                            self.strc[at1][bn] = at2
                            break

                for i in range(bond_type):
                    for bn in Strc.Bonds:
                        if self.strc[at2][bn] is np.nan:
                            self.strc[at2][bn] = at1
                            break

                self.meta.bond_from.append(at1)
                self.meta.bond_to.append(at2)
                self.meta.bond_type.append(parts[3])

                if len(parts) > 4:
                    self.meta.bond_status_bit.append(parts[4])
            print('OK')
        else:
            print('SKIP')


class Mol2Writer(Molecule):
    def __init__(self):
        super().__init__()

    def write_file(self, mol2file):
        with open(mol2file, 'w') as fp:
            self.molecule(fp)
            self.atoms(fp)
            self.bonds(fp)
            self.substructures(fp)

        print('Write OK:', mol2file)

    def total_bonds(self) -> int:
        return int(np.sum([
            len(np.unique(strc[Strc.Bonds][~np.isnan(strc[Strc.Bonds])]))
            for strc in self.strc
        ]) / 2)

    def molecule(self, fp):
        fp.write('@<TRIPOS>MOLECULE\n')
        molecstr = 	"{title}\n" \
                    "{no_atoms:>5d} {no_bonds:>5d} {no_residues:>5d} 0 0\n" \
                    "{type}\n" \
                    "{charge_type}\n\n"

        fp.write(molecstr.format(**dict(
            title = self.meta.get('title', 'Molecule'),
            no_atoms = self.info.shape[0],
            no_bonds = self.total_bonds(),
            no_residues = len(np.unique(self.strc[:,Strc.ResId])),
            type = self.meta.get('mol2_type', 'SMALL'),
            charge_type = self.meta.get('charge_type', 'USER_CHARGES'),
        )))

    def atoms(self, fp):
        fp.write('@<TRIPOS>ATOM\n')

        for i in range(self.info.shape[0]):
            info = self.info[i]
            crds = self.crds[i]
            strc = self.strc[i]

            atom = {
                'id' : i + 1,
                'name': self.info[i, Info.Name],
                'x':    self.crds[i, Crds.X],
                'y':    self.crds[i, Crds.Y],
                'z':    self.crds[i, Crds.Z],
                'type': self.info[i, Info.Type],
                'resid': int(self.strc[i, Strc.ResId]) + 1,
                'resname': self.info[i, Info.ResName],
                'charge': float(self.info[i, Info.Q])
            }
            atomstr =	"{id:>7d} {name:<5}  " \
                        "{x:>8.4f}  {y:>8.4f}  {z:>8.4f}   {type:>3} " \
                        "{resid:>3} {resname:<5}   {charge:>11.4f}\n"
            fp.write(atomstr.format(**atom))

    def bonds(self, fp):
        fp.write('@<TRIPOS>BOND\n')
        n = 0
        added = []
        for i in range(self.info.shape[0]):
            strc = self.strc[i, Strc.Bonds]
            for j, c in zip(*np.unique(strc[~np.isnan(strc)], return_counts=True)):
                if [i, j] in added or [j, i] in added:
                    continue
                n += 1
                bond = {
                    'id' : n,
                    'from': i + 1,
                    'to': int(j + 1),
                    'type': c,
                }

                added.append([i, j])

                bondstr =	"{id:>7d}  {from:>7d}  {to:>7d}   {type:>3} \n"
                fp.write(bondstr.format(**bond))

    def substructures(self, fp):
        fp.write('@<TRIPOS>SUBSTRUCTURE\n')

        start = np.nan
        n = 0
        for i in range(self.info.shape[0]):
            if self.strc[i, Strc.ResId] != start:
                n += 1
                start = self.strc[i, Strc.ResId]
                subst = {
                    'id' : n,
                    'name': self.info[i, Info.ResName],
                    'root': i + 1,
                    'type': 'RESIDUE',
                }
                resstr =	"{id:>7d}  {name:>7}  {root:>7d}   {type:>7} \n"
                fp.write(resstr.format(**subst))

    # @extend: add additional sections if needed

def read(mol2file : str) -> Molecule:
    mol = Mol2Reader()
    mol.read_file(mol2file)
    return Molecule.copy(mol)

def write(mol : Molecule, mol2file : str):
    mol = Mol2Writer.copy(mol)
    mol.write_file(mol2file)

if __name__ == '__main__':
    mol1 = read('gma.mol2')
    mol2 = Molecule.copy(mol1)
    copy = Molecule.copy(mol2)
    copy.translate(0, 0, 5)
    mol2.insert(copy)
    write(mol2, "test.mol2")
    breakpoint()
