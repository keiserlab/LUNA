from util import logging_ini

from collections import namedtuple
from os.path import (basename, splitext)
from interaction.contact import *

from bio.pdb_parser import PDBParser

parser = PDBParser(PERMISSIVE=1)

pdbFile = "../tmp/icode/1IL3.pdb"
pdbId = splitext(basename(pdbFile))[0]

structure = parser.get_structure(pdbId, pdbFile)
model = structure[0]

chain = 'A'
resname = "7DG"
resnum = 301
icode = ' '


# CHAIN
# entities = get_contacts_for_entity(model, model[chain], level='R')

# Ligand = model[chain][("H_%s" % resname, resnum, icode)]
# entities = get_contacts_for_entity(model, Ligand, level='A')

# Atom = model[chain][("H_%s" % resname, resnum, icode)]['N9']
# entities = get_contacts_for_entity(model, Atom, level='C')

# print(all_contacts_nh_search(model, level='C'))
