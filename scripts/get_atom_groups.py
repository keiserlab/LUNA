from rdkit import Chem
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.rdDepictor import Compute2DCoords

from MyBio.selector import (ResidueSelector, ResidueSelectorByResSeq)
from MyBio.PDB.PDBIO import PDBIO
from MyBio.PDB.PDBParser import PDBParser
from mol.chemical_feature import FeatureExtractor

from io import StringIO

from pybel import readstring

from interaction.contact import *


fdefName = 'data/MinimalFeatures.fdef'
featFactory = ChemicalFeatures.BuildFeatureFactory(fdefName)
fExt = FeatureExtractor(featFactory)

pdbPath = "../tmp/pharm"

pdbId = "3QQK"
pdbFile = "%s/%s.pdb" % (pdbPath, pdbId)

parser = PDBParser(PERMISSIVE=True)

struct = parser.get_structure(pdbId, pdbFile)
chainA = struct[0]["A"]

lig = ("H_X02", 497, " ")

targetMolecules = set([497])

resSel = ResidueSelectorByResSeq(targetMolecules)

io = PDBIO()
io.set_structure(chainA)


targetMolecules = [chainA[2], chainA[8], chainA[20], chainA[82], chainA[lig]]
# targetMolecules = [chainA[lig]]

path = '.'
for res in list(targetMolecules):

    proximal = get_contacts_for_entity(model=struct[0], source=res,
                                       radius=1.99999999999, level='R')

    covResidues = set()
    for pair in proximal:

        if (pair[1] != res):
            covResidues.add(pair[1])

    resSel = ResidueSelector(list(covResidues) + [res])

    fh = StringIO()
    io.save(fh, select=resSel, preserve_atom_numbering=True,
            write_conects=True)
    fh.seek(0)
    pdbBlock = ''.join(fh.readlines())

    obMol = readstring("pdb", pdbBlock)
    obMol.OBMol.AddHydrogens(False, True, 7)

    molBlock = obMol.write('mol')

    rdMol = Chem.MolFromMolBlock(molBlock)

    Compute2DCoords(rdMol)

    atomFeatures = fExt.get_atom_features(rdMol)

    drawer = rdMolDraw2D.MolDraw2DSVG(800, 600)
    drawer.DrawMolecule(rdMol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace('svg:', '')

    output = "%s/%s.svg" % (path, res.get_id()[1])
    with open(output, "w") as fh:
        fh.write(svg)
