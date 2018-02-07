from util import logging_ini
from align.tmalign import (align_2struct, extract_chain_from_sup,
                           remove_sup_files)

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO

from os.path import basename

from bio.pdb import Extractor

parser = PDBParser(PERMISSIVE=True, QUIET=True)

# f1 = '../tmp/alignment/1BGA.A.sup_all_atm_lig'
# structure = parser.get_structure("TESTE", f1)
# model = structure[0]
# model.detach_child('B')
# chain = model['A']
# chain.id = "B"
# io = PDBIO()
# io.set_structure(structure)
# io.save('../tmp/alignment/1BGA.A.aligned.pdb', preserve_atom_numbering=True)

f1 = '../tmp/icode/1IL3.pdb'
f2 = '../tmp/icode/4IMV.pdb'
outputPath = 'obsolete'
tmalign = '../../Softwares/TMalignc/TMalign'
alignment = align_2struct(f1, f2, outputPath, tmalign=tmalign)

supFile = "%s/1IL3.pdb.sup_all_atm_lig" % outputPath
renameChain = ""
outputFile = "obsolete/1IL3.aligned.pdb"
extract_chain_from_sup(supFile, "A", "E", outputFile, QUIET=True)

remove_sup_files(outputPath)

chain = 'G'
f3 = '../tmp/alignment/1BGA.pdb'
m3 = parser.get_structure("algo", f3)[0]
extractor = Extractor(m3)
outF3 = '../tmp/icode/1BGA.chain%s.pdb' % chain
extractor.extract_residues([chain], outF3)
