from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.AllChem import Compute2DCoords

from math import (cos, sin, radians)


class ColorPallete:

    def __init__(self, color_map=None):
        color_map = color_map or {}
        self.color_map = color_map

    def add_color(self, key, rgb_tuple):
        self.color_map[key] = rgb_tuple

    def get_normalized_color(self, key):
        return tuple(c / 255 for c in self.color_map[key])


def add_atom(mol, pos=None):
    new_atm = Chem.rdchem.Atom(10)
    new_atm.SetNoImplicit(True)

    atm_id = mol.AddAtom(new_atm)

    mol.GetConformer().SetAtomPosition(atm_id, pos)

    return atm_id


def ligand_pharm_figure(rdmol, atm_types, output, colors=None):
    # Make a copy of the molecule
    rwm = Chem.RWMol(rdmol)
    Compute2DCoords(rwm)

    drawer = rdMolDraw2D.MolDraw2DSVG(800, 600)
    opts = drawer.drawOptions()

    highlight = {}
    cutoff = 0.2

    for atm_id in atm_types:
        features = list(atm_types[atm_id])
        centroid = list(rwm.GetConformer().GetAtomPosition(atm_id))

        if (len(features) == 1):
            pos = centroid
            new_atm_id = add_atom(rwm, centroid)
            highlight[new_atm_id] = colors.get_normalized_color(features[0].name)
            opts.atomLabels[new_atm_id] = ''
        else:
            sliceRad = radians(360 / len(features))
            for i, feature in enumerate(features):
                rad = i * sliceRad
                pos = [cutoff * cos(rad), cutoff * sin(rad), 0]
                adj_Pos = [x + y for x, y in zip(centroid, pos)]
                new_atm_id = add_atom(rwm, adj_Pos)
                highlight[new_atm_id] = colors.get_normalized_color(feature.name)
                opts.atomLabels[new_atm_id] = ''

    atoms = [x for x in highlight]
    opts.flagCloseContactsDist = -1000
    opts.legendFontSize = 20
    opts.padding = 0.2
    opts.useBWAtomPalette()

    drawer.SetFontSize(0.5)
    drawer.DrawMolecule(rwm, highlightAtoms=atoms,
                        highlightAtomColors=highlight,
                        highlightBonds=[],
                        legend=rdmol.GetProp("_Name"))
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace('svg:', '')

    with open(output, "w") as fh:
        fh.write(svg)
