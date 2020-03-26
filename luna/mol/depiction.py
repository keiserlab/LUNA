from math import cos, sin, radians

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.AllChem import Compute2DCoords

from luna.mol.wrappers.base import MolWrapper
from luna.util.default_values import ATOM_TYPES_COLOR


class PharmacophoreDepiction:

    def __init__(self, feature_extractor=None, colors=ATOM_TYPES_COLOR, add_legend=True, fig_ext="png",
                 fig_size=(800, 800), font_size=0.5, circle_dist=0.2, circle_radius=0.3,
                 use_bw_atom_palette=True, svg_opts=None):

        self.feature_extractor = feature_extractor
        self.colors = colors
        self.add_legend = add_legend
        self.fig_ext = fig_ext
        self.fig_size = fig_size
        self.font_size = font_size
        self.circle_dist = circle_dist
        self.circle_radius = circle_radius
        self.use_bw_atom_palette = use_bw_atom_palette

        if svg_opts is None:
            svg_opts = {
                "flagCloseContactsDist": -1000,
                "legendFontSize": 20,
                "padding": 0.2
            }
        self.svg_opts = svg_opts or {}

    def _perceive_atm_types(self, rdmol):
        return self.feature_extractor.get_features_by_atoms(rdmol)

    def plot_fig(self, mol, output, atm_types=None, legend=None):
        rdmol = MolWrapper(mol).as_rdkit()

        # Make a copy of the molecule
        rwm = Chem.RWMol(rdmol)
        Compute2DCoords(rwm)

        if atm_types is None:
            atm_types = self._perceive_atm_types(rdmol)

        if self.fig_ext == "png":
            drawer = rdMolDraw2D.MolDraw2DCairo(*self.fig_size)
        elif self.fig_ext == "svg":
            drawer = rdMolDraw2D.MolDraw2DSVG(*self.fig_size)
        else:
            # raise
            pass

        opts = drawer.drawOptions()

        highlight = {}

        for atm_id in atm_types:
            centroid = list(rwm.GetConformer().GetAtomPosition(atm_id))
            valid_features = [f for f in atm_types[atm_id] if f.name in self.colors]

            if valid_features:
                if len(valid_features) == 1:
                    pos = centroid
                    atmIdx = self._add_dummy_atom(rwm, centroid)
                    highlight[atmIdx] = self.colors.get_normalized_color(valid_features[0].name)
                    opts.atomLabels[atmIdx] = ''
                else:
                    sliceRad = radians(360 / len(valid_features))
                    for i, feature in enumerate(valid_features):
                        rad = i * sliceRad
                        pos = [self.circle_dist * cos(rad), self.circle_dist * sin(rad), 0]
                        adj_Pos = [x + y for x, y in zip(centroid, pos)]
                        new_atm_id = self._add_dummy_atom(rwm, adj_Pos)
                        highlight[new_atm_id] = self.colors.get_normalized_color(feature.name)
                        opts.atomLabels[new_atm_id] = ''

        atoms = [x for x in highlight]
        radius = {a: self.circle_radius for a in atoms}

        # if self.svg_opts:
        #     for k, v in self.svg_opts:
        #     opts = **self.svg_opts

        opts.flagCloseContactsDist = -1000
        opts.legendFontSize = 20
        opts.padding = 0.02

        if self.use_bw_atom_palette:
            opts.useBWAtomPalette()
        else:
            opts.useDefaultAtomPalette()

        legend = legend or rdmol.GetProp("_Name")

        drawer.SetFontSize(self.font_size)
        drawer.DrawMolecule(rwm, highlightAtoms=atoms, highlightAtomColors=highlight,
                            highlightBonds=[], highlightAtomRadii=radius, legend=legend)
        drawer.FinishDrawing()

        if self.fig_ext == "png":
            drawer.WriteDrawingText(output)
        elif self.fig_ext == "svg":
            svg = drawer.GetDrawingText().replace('svg:', '')
            with open(output, "w") as fh:
                fh.write(svg)

    def _add_dummy_atom(self, mol, pos=None):
        new_atm = Chem.rdchem.Atom(10)
        new_atm.SetNoImplicit(True)
        atm_id = mol.AddAtom(new_atm)
        mol.GetConformer().SetAtomPosition(atm_id, pos)

        return atm_id
