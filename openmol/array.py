import numpy as np

from .utils import AttrDict

## Everything is per atom basis.
class Info:
    _n_cols =  4

    Name    =  0
    Type    =  1
    Q       =  2
    ResName =  3

class Crds:
    _n_cols =  3

    X       =  0
    Y       =  1
    Z       =  2

class Strc:
    _n_cols =  10

    ResId   =  0
    Bonds   = [1, 2, 3, 4, 5, 6, 7, 8]  # Max 8 bonds allowed.
    SegId   =  9


class MolReader:
    """ Base file reader interface to implement in different
        Reader classes. """

    def __init__(self):
        self.info = []
        self.crds = []
        self.strc = []
        self.meta = AttrDict()

        self.in_file = None
        self.lines = []
        self.section = None
        self.section_start = 0
        self.section_lines = []
        self.section_format = None

    def new_atom(self) -> tuple:
        return (
            [np.nan] * Info._n_cols,
            [np.nan] * Crds._n_cols,
            [np.nan] * Strc._n_cols
        )

    def add_atom(self, info, crds, strc):
        assert len(info) == Info._n_cols
        assert len(crds) == Crds._n_cols
        assert len(strc) == Strc._n_cols

        self.info.append(info)
        self.crds.append(crds)
        self.strc.append(strc)

    def finalize(self):
        self.info = np.array(self.info)
        self.crds = np.array(self.crds)
        self.strc = np.array(self.strc)

    def read_file(self, input_file : str = None):
        if input_file is not None:
            self.in_file = input_file
        self.lines = open(self.in_file).readlines()
        print('Reading:', input_file)
        self._process_lines()
        print('Read OK:', input_file)
        self.finalize()

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


class Molecule:
    def __init__(self):
        self.info : np.array = None
        self.crds : np.array = None
        self.strc : np.array = None
        self.meta : AttrDict = None

    @classmethod
    def copy(cls, reader : MolReader):
        self = cls()
        self.info = np.copy(reader.info)
        self.crds = np.copy(reader.crds)
        self.strc = np.copy(reader.strc)
        self.meta = AttrDict(dict(reader.meta))

        assert len(self.info) == len(self.crds)
        assert len(self.info) == len(self.strc)
        assert len(self.crds) == len(self.strc)
        return self

    def write_file(self, input_file : str):
        raise NotImplementedError()

    def insert(self, molecule : 'Molecule'):
        molecule = Molecule.copy(molecule)
        
        incr = np.max(self.strc[:, Strc.SegId])
        molecule.strc[:, Strc.SegId] += incr + 1

        incr = np.max(self.strc[:, Strc.ResId])
        molecule.strc[:, Strc.ResId] += incr + 1

        incr = len(self.info)
        molecule.strc[:, Strc.Bonds] += incr

        # Concat
        self.info = np.concatenate([self.info, molecule.info], axis=0)
        self.crds = np.concatenate([self.crds, molecule.crds], axis=0)
        self.strc = np.concatenate([self.strc, molecule.strc], axis=0)

    def translate(self, x, y, z):
        self.crds[:, Crds.X] += x
        self.crds[:, Crds.Y] += y
        self.crds[:, Crds.Z] += z

    def get_atoms(self, name : str = None, atype : str = None, resname : str = None,
                    resid : int = None, segid : int = None) -> np.array:
        mask = [True] * len(self.info)
        if name:
            mask = mask & (self.info[:, Info.Name] == name)
        if atype:
            mask = mask & (self.info[:, Info.Type] == atype)
        if resname:
            mask = mask & ((self.info[:, Info.ResName] == resname))
        if resid:
            mask = mask & ((self.strc[:, Strc.ResId] == resid))
        if segid:
            mask = mask & ((self.strc[:, Strc.SegId] == segid))

        return np.where(mask)
