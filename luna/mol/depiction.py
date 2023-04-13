from math import cos, sin, radians

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.AllChem import Compute2DCoords

from luna.wrappers.base import MolWrapper
from luna.util.default_values import ATOM_TYPES_COLOR


class PharmacophoreDepiction:
    """Draw molecules and depict pharmacophoric properties as colored circles.

    Parameters
    ----------
    feature_extractor : :class:`~luna.mol.features.FeatureExtractor`
        Perceive pharmacophoric properties from molecules.
    colors : :class:`~luna.util.ColorPallete`
        Color scheme for pharmacophoric properties perceived by
        ``feature_extractor``. The default value is \
        :const:`~luna.util.default_values.ATOM_TYPES_COLOR`.
    format : {'png', 'svg'}
        The output file format. The default value is 'png'.
    figsize : tuple of (float, float)
        Width and height in inches. The default value is (800, 800).
    font_size : float
        The font size. The units are, roughly, pixels.
        The default value is 0.5.
    circle_dist : float
        Distance between circles (pharmacophoric properties).
        The default value is 0.2.
    circle_radius :
        Circles' radius size of pharmacophoric properties.
        The default value is 0.3.
    use_bw_atom_palette : bool
        Use a black & white palette for atoms and bonds.

    Examples
    --------

    First, let's read a molecule (glutamine).

    >>> from luna.wrappers.base import MolWrapper
    >>> mol = MolWrapper.from_smiles("N[C@@H](CCC(N)=O)C(O)=O")

    Now, create a feature factory and instantiate a new FeatureExtractor
    object.

    >>> from luna.util.default_values import ATOM_PROP_FILE
    >>> from rdkit.Chem import ChemicalFeatures
    >>> from luna.mol.features import FeatureExtractor
    >>> feature_factory = ChemicalFeatures.BuildFeatureFactory(ATOM_PROP_FILE)
    >>> feature_extractor = FeatureExtractor(feature_factory)

    Instantiate a new PharmacophoreDepiction object with the desired
    configuration. For example, you can provide a color scheme for
    pharmacophoric properties, the image size, and its format.

    >>> from luna.util.default_values import ATOM_TYPES_COLOR
    >>> from luna.mol.depiction import PharmacophoreDepiction
    pd = PharmacophoreDepiction(feature_extractor=feature_extractor,
                                colors=ATOM_TYPES_COLOR,
                                fig_size=(500, 500),
                                format="svg")

    Finally, you can draw the molecule with annotated pharmacophoric
    properties.

    >>> pd.plot_fig(mol, "output.svg")
    """

    def __init__(self,
                 feature_extractor=None,
                 colors=ATOM_TYPES_COLOR,
                 format="png",
                 fig_size=(800, 800),
                 font_size=0.5,
                 circle_dist=0.2,
                 circle_radius=0.3,
                 use_bw_atom_palette=True):

        self.feature_extractor = feature_extractor
        self.colors = colors
        self.format = format
        self.fig_size = fig_size
        self.font_size = font_size
        self.circle_dist = circle_dist
        self.circle_radius = circle_radius
        self.use_bw_atom_palette = use_bw_atom_palette

    def _perceive_atm_types(self, rdmol):
        if self.feature_extractor is not None:
            return self.feature_extractor.get_features_by_atoms(rdmol)
        return {}

    def plot_fig(self,
                 mol_obj,
                 output=None,
                 atm_types=None,
                 legend=None):
        """Draw the molecule ``mol_obj`` and depict its pharmacophoric
        properties.

        Parameters
        ----------
        mol_obj : :class:`~luna.wrappers.base.MolWrapper`, \
                    :class:`rdkit.Chem.rdchem.Mol`, or \
                    :class:`openbabel.pybel.Molecule`
            The molecule.
        output : str
            The output file where the molecule will be drawn.
            If None, returns a drawing object
            (:class:`~rdkit.Chem.Draw.rdMolDraw2D.MolDraw2DCairo` or
            :class:`~rdkit.Chem.Draw.rdMolDraw2D.MolDraw2DSVG`).
        atm_types : dict or None
            A pre-annotated dictionary for mapping atoms and pharmacophoric
            properties. If None, try to perceive properties with
            ``feature_extractor``.
        legend : str
            A title for the figure.

        Returns
        -------
        drawer : None or a drawing object \
            (:class:`~rdkit.Chem.Draw.rdMolDraw2D.MolDraw2DCairo` or \
             :class:`~rdkit.Chem.Draw.rdMolDraw2D.MolDraw2DSVG`)
        """

        rdmol = MolWrapper(mol_obj).as_rdkit()

        # Make a copy of the molecule
        rwm = Chem.RWMol(rdmol)
        Compute2DCoords(rwm)

        if atm_types is None:
            atm_types = self._perceive_atm_types(rdmol)

        if self.format == "png":
            drawer = rdMolDraw2D.MolDraw2DCairo(*self.fig_size)
        elif self.format == "svg":
            drawer = rdMolDraw2D.MolDraw2DSVG(*self.fig_size)
        else:
            # raise
            pass

        opts = drawer.drawOptions()

        highlight = {}
        for atm_id in atm_types:
            centroid = list(rwm.GetConformer().GetAtomPosition(atm_id))
            valid_features = [f for f in atm_types[atm_id]
                              if f.name in self.colors]

            if valid_features:
                if len(valid_features) == 1:
                    pos = centroid
                    atmIdx = self._add_dummy_atom(rwm, centroid)
                    feat_name = valid_features[0].name
                    highlight[atmIdx] = \
                        self.colors.get_normalized_color(feat_name)
                    opts.atomLabels[atmIdx] = ''
                else:
                    sliceRad = radians(360 / len(valid_features))
                    for i, feature in enumerate(valid_features):
                        rad = i * sliceRad
                        pos = [self.circle_dist * cos(rad),
                               self.circle_dist * sin(rad),
                               0]
                        adj_Pos = [x + y for x, y in zip(centroid, pos)]
                        new_atm_id = self._add_dummy_atom(rwm, adj_Pos)
                        highlight[new_atm_id] = \
                            self.colors.get_normalized_color(feature.name)
                        opts.atomLabels[new_atm_id] = ''

        atoms = [x for x in highlight]
        radius = {a: self.circle_radius for a in atoms}

        opts.flagCloseContactsDist = -1000
        opts.legendFontSize = 20
        opts.padding = 0.02

        if self.use_bw_atom_palette:
            opts.useBWAtomPalette()
        else:
            opts.useDefaultAtomPalette()

        legend = legend or ""

        drawer.SetFontSize(self.font_size)
        drawer.DrawMolecule(rwm,
                            highlightAtoms=atoms,
                            highlightAtomColors=highlight,
                            highlightBonds=[],
                            highlightAtomRadii=radius,
                            legend=legend)
        drawer.FinishDrawing()

        if output:
            if self.format == "png":
                drawer.WriteDrawingText(output)
            elif self.format == "svg":
                svg = drawer.GetDrawingText().replace('svg:', '')
                with open(output, "w") as fh:
                    fh.write(svg)
        else:
            return drawer

    def _add_dummy_atom(self, mol, pos=None):
        new_atm = Chem.rdchem.Atom(10)
        atm_id = mol.AddAtom(new_atm)
        mol.GetConformer().SetAtomPosition(atm_id, pos)
        new_atm.SetNoImplicit(True)

        return atm_id
