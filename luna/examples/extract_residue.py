from MyBio.PDB.PDBParser import PDBParser
from MyBio.extractor import Extractor

parser = PDBParser(PERMISSIVE=True)

f1 = '../tmp/icode/1IL3.pdb'
m = parser.get_structure("algo", f1)[0]
extractor = Extractor(m)

ligand = ('H_7DG', 301, ' ')

out = '../tmp/icode/1BGA.7DG.pdb'

extractor.set_entity(m['A'])
extractor.extract_residues([ligand], out)
