from input.entry_validator import *
from input.util import entry_format
from input.existance_validator import *

from os.path import (basename, splitext)
from MyBio.PDB.PDBParser import PDBParser

parser = PDBParser(PERMISSIVE=1)

pdbFile = "../tmp/icode/1IL3.pdb"
pdbId = splitext(basename(pdbFile))[0]

structure = parser.get_structure(pdbId, pdbFile)
model = structure[0]

pdbId = pdbId
chain = 'A'
resname = "7DG"
resnum = 301
icode = ' '

# Entry format
# print(entry_format([pdbId, chain]))

# PLI validator
# entry = ":".join([str(x) for x in (pdbId, chain, resname, resnum)])
# entryValidator = PliEntryValidator()

# PPI validator
# entry = ":".join([str(x) for x in (pdbId, chain)])
# entryValidator = PpiEntryValidator()

# PLI validator - Upload files
# pdbId = "filename-bem-grande_e_com_underscore"
# entry = ":".join([str(x) for x in (pdbId, chain, resname, resnum)])
# entryValidator = PliEntryValidator(True)

# PPI validator - Upload files
# pdbId = "filename-bem-grande_e_com_underscore"
# entry = ":".join([str(x) for x in (pdbId, chain)])
# entryValidator = PpiEntryValidator(True)

# print(entryValidator.is_valid(entry))



# Check molecule existance
# Not exists
# mol = ("W", 298, " ")
# mol = ("H_7DG", 298, " ")
# mol = ("H_7DG", 298, "A")

# Exists
# mol = ("H_7DG", 301, " ")
# mol = ("W", 302, " ")

# ICODE for 3V0W
# PDB with ICODE
# pdbFile = "../tmp/icode/3V0W.pdb"
# pdbId = splitext(basename(pdbFile))[0]
# structure = parser.get_structure(pdbId, pdbFile)
# model = structure[0]
# chain = "H"
# mol = (" ", 52, "D")

# exists, message = check_molecule_existance(model, chain, mol)
# print(exists)
# print(message)
