from rdkit import Chem
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.AllChem import Compute2DCoords
from rdkit.Chem.Draw import rdMolDraw2D
from math import (cos, sin, radians)
from file.util import *

import glob


def add_atom(mol, pos=None):
    newAtom = Chem.rdchem.Atom(10)
    newAtom.SetNoImplicit(True)

    atmIdx = rwm.AddAtom(newAtom)

    rwm.GetConformer().SetAtomPosition(atmIdx, pos)

    return atmIdx

fdefName = 'data/MinimalFeatures.fdef'
featFactory = ChemicalFeatures.BuildFeatureFactory(fdefName)

path = "../tmp/clustering"
molFiles = glob.glob('%s/ligands_mol/*.mol' % path)
mols = [Chem.MolFromMolFile(m) for m in molFiles]

path = "../tmp/pharm"

for m in mols:
    Compute2DCoords(m)

tmpColors1 = {
    "Acceptor": [102,194,165],
    "Donor": [255,217,47],
    "Aromatic": [231,138,195],
    "Hydrophobic": [166,216,84],
    "PosIonizable": [141,160,203],
    "NegIonizable": [252,141,98]
}

tmpColors2 = {
    "Acceptor": [253,180,98],
    "Donor": [255,255,179],
    "Aromatic": [190,186,218],
    "Hydrophobic": [141,211,199],
    "PosIonizable": [128,177,211],
    "NegIonizable": [251,128,114]
}

tmpColorsSafeColorBlind1 = {
    "Acceptor": [209,229,240],
    "Donor": [103,169,207],
    "Aromatic": [253,219,199],
    "Hydrophobic": [239,138,98],
    "PosIonizable": [33,102,172],
    "NegIonizable": [178,24,43]
}

tmpColorsSafeColorBlind2 = {
    "Acceptor": [252,141,89],
    "Donor": [145,191,219],
    "Aromatic": [224,243,248],
    "Hydrophobic": [254,224,144],
    "PosIonizable": [69,117,180],
    "NegIonizable": [215,48,39]
}

tmpColors = tmpColorsSafeColorBlind2

colors = {}
for t in tmpColors:
    for i, c in enumerate(tmpColors[t]):
        tmpColors[t][i] = c / 255

    colors[t] = tuple(tmpColors[t])


for mol in mols:
    rwm = Chem.RWMol(mol)

    drawer = rdMolDraw2D.MolDraw2DSVG(800, 600)
    opts = drawer.drawOptions()


    highlight = {}
    addAtom = []

    cutoff = 0.2

    atomTypes = {}
    for f in featFactory.GetFeaturesForMol(mol):
        atoms = set(f.GetAtomIds())

        for a in atoms:
            types = atomTypes[a] if (a in atomTypes) else []
            types.append(f.GetFamily())
            atomTypes[a] = types

    for a in atomTypes:
        types = atomTypes[a]
        centroid = list(rwm.GetConformer().GetAtomPosition(a))

        if (len(types) == 1):
            pos = centroid
            atmIdx = add_atom(rwm, centroid)
            highlight[atmIdx] = colors[types[0]]
            opts.atomLabels[atmIdx] = ''
        else:
            sliceRad = radians(360 / len(types))
            for i, atmType in enumerate(types):
                rad = i * sliceRad
                pos = [cutoff * cos(rad), cutoff * sin(rad), 0]
                adjustedPos = [x + y for x, y in zip(centroid, pos)]
                atmIdx = add_atom(rwm, adjustedPos)
                highlight[atmIdx] = colors[atmType]
                opts.atomLabels[atmIdx] = ''

        # atmIdx = add_atom(rwm, centroid)
        # highlight[atmIdx] = (1,1,1)
        # opts.atomLabels[atmIdx] = ''


    atoms = [x for x in highlight]

    opts.flagCloseContactsDist = -1000

    opts.legendFontSize = 20
    opts.padding = 0.2

    opts.useBWAtomPalette()

    entry = get_filename(mol.GetProp("_Name"))

    drawer.SetFontSize(0.5)
    drawer.DrawMolecule(rwm, highlightAtoms=atoms,
                        highlightAtomColors=highlight,
                        highlightBonds=[],
                        legend=entry)

    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace('svg:', '')

    output = "%s/%s.svg" % (path, entry)

    with open(output, "w") as fh:
        fh.write(svg)