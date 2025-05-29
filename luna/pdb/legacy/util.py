import tempfile
import logging
from io import StringIO
from os.path import exists
from shutil import move as rename_pdb_file


from openbabel.pybel import readfile
from openbabel.pybel import Molecule as PybelWrapper

from rdkit.Chem import (MolFromMolBlock, SanitizeFlags,
                        SanitizeMol, MolToMolFile)

from luna.util.file import (is_directory_valid, new_unique_filename,
                            remove_files)
from luna.util.default_values import ENTRY_SEPARATOR, OPENBABEL
from luna.pdb.selector import AtomSelector
from luna.pdb.PDB.PDBList import PDBList
from luna.pdb.PDB.PDBIO import Select
from luna.pdb.PDB import Selection
from luna.pdb.PDB.PDBIO import PDBIO
from luna.pdb.PDB.PDBParser import PDBParser
from luna.mol.validator import MolValidator
from luna.mol.standardiser import Standardizer
from luna.wrappers.base import MolWrapper
from luna.wrappers.obabel import convert_molecule
from luna.wrappers.rdkit import read_mol_from_file

from luna.util.exceptions import (IllegalArgumentError, MoleculeNotFoundError,
                                  ChainNotFoundError, FileNotCreated,
                                  PDBNotReadError, MoleculeSizeError,
                                  MoleculeObjectError)


logger = logging.getLogger()


ENTITY_LEVEL_NAME = {"A": "Atom",
                     "R": "Residue",
                     "C": "Chain",
                     "M": "Model",
                     "S": "Structure"}







