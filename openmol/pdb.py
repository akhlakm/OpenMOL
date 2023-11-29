from openmol import core

class PDBReader(core.Reader):
    def __init__(self):
        super().__init__()
        self.Mol['source_format'] = "PDB"
        self.Mol['_pdb_built'] = False
        self.Mol['comment'] = ""
        self.Mol['atom_segment'] = []
        self.Mol['atom_chain'] = []
        self.Mol['atom_altloc'] = []
        self.Mol['atom_icode'] = []
        self.Mol['atom_occupancy'] = []
        self.Mol['atom_temp_factor'] = []
        self.Mol['bond_tacticity'] = []

    def read_file(self, pdb_file : str):
        super().read_file(pdb_file)

    def _identify_section(self, line_no: int, line: str, next_line : str):
        words = line.split()
        nextwords = next_line.split() if next_line else []
        if words[0] == 'REMARK' and self.section is None:
            # header/title
            self._new_section('REMARK', line_no)

        elif words[0] == 'CRYST1':
            self._new_section('CRYST1', line_no)

        elif words[0] == 'ATOM':
            self._new_section('ATOM', line_no)

        elif words[0] == 'REMARK':
            # look ahead
            if len(nextwords) > 1 and nextwords[1] == "CONECT":
                self._new_section('CONECT', line_no)
            elif len(nextwords) > 1 and nextwords[1] == "RESCON":
                self._new_section('RESCON', line_no)
            elif len(nextwords) > 1 and nextwords[1] == "SEGRNG":
                self._new_section('SEGRNG', line_no)
        
        elif words[0] == 'END':
            self._new_section('END', line_no)


    def _process_last_section(self, section: str, lines: list, sformat: str):
        if section == 'REMARK':
            self.Mol.comment += " ".join(line[7:].strip() for line in lines)
            print("OK")

        elif section == 'CRYST1':
            words = lines[0].split()
            self.Mol['box_x'] = self._str_to_type(
                                words[1], float, lines[0], self.section_start)
            self.Mol['box_y'] = self._str_to_type(
                                words[2], float, lines[0], self.section_start)
            self.Mol['box_z'] = self._str_to_type(
                                words[3], float, lines[0], self.section_start)
            self.Mol['box_alpha'] = self._str_to_type(
                                words[4], float, lines[0], self.section_start)
            self.Mol['box_beta'] = self._str_to_type(
                                words[5], float, lines[0], self.section_start)
            self.Mol['box_gamma'] = self._str_to_type(
                                words[6], float, lines[0], self.section_start)
            print("OK")

        elif section == 'ATOM':
            for i, line in enumerate(lines):
                ln = self.section_start + i
                # Try to follow the wwpdb standard.
                # https://files.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_Letter.pdf
                try:
                    x = self._str_to_type(line[30:38].strip(), float, line, ln)
                    y = self._str_to_type(line[38:46].strip(), float, line, ln)
                    z = self._str_to_type(line[46:54].strip(), float, line, ln)

                    name = line[12:16].strip()
                    resname = line[17:20].strip()
                    resid = self._str_to_type(line[22:26].strip(),
                                              int, line, ln) - 1
                    
                    at_type = line[76:78].strip()
                    if " " in at_type:
                        raise ValueError('invalid element format')
                    
                    try:
                        q = self._str_to_type(
                            line[78:80].strip(), float, line, ln)
                    except:
                        q = 0.0

                    self.Mol.atom_name.append(name)
                    self.Mol.atom_type.append(at_type)
                    self.Mol.atom_resname.append(resname)
                    self.Mol.atom_resid.append(resid)
                    self.Mol.atom_x.append(x)
                    self.Mol.atom_y.append(y)
                    self.Mol.atom_z.append(z)
                    self.Mol.atom_q.append(q)

                    try:
                        self.Mol.atom_occupancy.append(self._str_to_type(
                            line[54:60].strip(), float, line, ln))
                    except:
                        print("failed to parse occupancy")

                    try:
                        self.Mol.atom_temp_factor.append(self._str_to_type(
                            line[60:66].strip(), float, line, ln))
                    except:
                        print("failed to parse tempFactor")

                    self.Mol.atom_chain.append(line[21].strip())
                    self.Mol.atom_altloc.append(line[16].strip())
                    self.Mol.atom_icode.append(line[26].strip())
                    self.Mol.atom_segment.append(line[66:77].strip())

                except:
                    raise ValueError("Non-standard PDB format")
                    # words = line.split()
                    # self.Mol['atom_name'] = words[2]
                    # self.Mol['atom_resname'] = words[3]
                    # self.Mol['atom_chainid'] = words[4]

            print("OK")

        elif section == 'CONECT':
            unique = []
            for i, line in enumerate(lines):
                ln = self.section_start + i
                words = line.strip().split()
                if words[1] != 'CONECT': continue
                bn_from = self._str_to_type(words[2], int, line, ln) - 1
                for other in words[3:]:
                    bn_type = 1
                    tacticity = 0
                    if ":" in other and other.count(":") == 1:
                        bn_to, bn_type = other.split(":")
                    elif ":" in other and other.count(":") == 2:
                        bn_to, bn_type, tacticity = other.split(":")
                    else:
                        bn_to = other

                    bn_to = self._str_to_type(bn_to, int, line, ln) - 1
                    if bn_type == "1.5":
                        bn_type = 'ar'
                    else:
                        bn_type = self._str_to_type(bn_type, int, line, ln) # type: ignore
                    forward = (bn_from, bn_to, bn_type, tacticity)
                    reverse = (bn_to, bn_from, bn_type, tacticity)
                    if reverse not in unique:
                        unique.append(forward)

            for item in unique:
                self.Mol.bond_from.append(item[0])
                self.Mol.bond_to.append(item[1])
                self.Mol.bond_type.append(item[2])
                self.Mol.bond_tacticity.append(item[3])

            print("OK")

        elif section == 'END':
            print("OK")

        else:
            print("IGNORED")


class PDB:
    def __init__(self) -> None:
        self.Mol = core.AttrDict()

    def build(self):
        unique_resids = list(set(self.Mol.atom_resid))

        # Fix atom resids order.
        # resid should increase all the times.
        if len(unique_resids) > 1:
            offset = 0
            old_id = 0
            for i, resid in enumerate(self.Mol.atom_resid):
                if old_id > resid and resid == 0:
                    print('- ResID restarted at atom %d' %(i+1))
                    offset += old_id - resid + 1
                elif old_id > resid and resid != 0:
                    ValueError("Invalid ResID order at atom %d" %(i+1))

                old_id = resid
                self.Mol.atom_resid[i] = resid + offset

        self.Mol['_pdb_built'] = True


    def read(self, file_path : str):
        reader = PDBReader()
        reader.read_file(file_path)
        self.Mol = reader.Mol
        self.build()


if __name__ == '__main__':
    r = PDB()
    r.read("test.pdb")
    # for k, v in r.Mol.items():
    #     if v: print(k, v)

    from openmol.tripos_mol2 import Writer, build
    wr = Writer(build(r.Mol), "test.mol2")
    wr.write()
