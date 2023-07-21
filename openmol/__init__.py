__version__ = "1.2.0"
__author__ 	= "Akhlak Mahmood, Yingling Group, MSE, NCSU"

from .utils import AttrDict
from .openmol import initialize, Writer
from .openmol import check, update_summary
from .openmol import check_atoms_ok, check_bonds_ok, check_residues_ok

del openmol, utils
