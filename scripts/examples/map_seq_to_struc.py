from util import logging_ini
from align.tmalign import align_2struct

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import StructureAlignment

import os.path

parser = PDBParser(PERMISSIVE=1)

f1 = '../tmp/icode/1IL3.pdb'
f2 = '../tmp/icode/4IMV.pdb'

tmalign = '../../Softwares/TMalignc/TMalign'

alignment = align_2struct(f1, f2, tmalign)

m1 = parser.get_structure(os.path.basename(f1), f1)[0]
m2 = parser.get_structure(os.path.basename(f2), f2)[0]


mapping = StructureAlignment(alignment, m1, m2)

print(mapping)

print(mapping.map12)