import re

from luna.wrappers.base import MolWrapper
from luna.util.default_values import MIN_FDEF_FILE
from luna.util.exceptions import FingerprintNotCreated, IllegalArgumentError

from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AtomPairs import Pairs, Torsions
from rdkit.Chem import AllChem
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Pharm2D import Generate
from rdkit.Chem.Pharm2D.SigFactory import SigFactory

import logging

logger = logging.getLogger()


class FingerprintGenerator():
    """Generate molecular fingerprints for the molecule ``mol``.

    Parameters
    ----------
    mol : :class:`~luna.wrappers.base.MolWrapper`, \
            :class:`rdkit.Chem.rdchem.Mol`, or \
            :class:`openbabel.pybel.Molecule`, optional
        The molecule.

    Examples
    --------

    First, create a new `FingerprintGenerator` object.

    >>> from luna.mol.fingerprint import FingerprintGenerator
    >>> fg = FingerprintGenerator()

    Now, let's read a molecule (glutamine) and set it to the
    `FingerprintGenerator` object.

    >>> from luna.wrappers.base import MolWrapper
    >>> fg.mol = MolWrapper.from_smiles("N[C@@H](CCC(N)=O)C(O)=O")

    Finally, you can call any available function to generate the desired
    fingerprint type. In the below example, a count ECFP4 fingerprint of
    size 1,024 is created.

    >>> fp = fg.morgan_fp(radius=2, length=1024, type=2)
    >>> print(fp.GetNonzeroElements())
    {1: 1, 80: 2, 140: 1, 147: 2, 389: 1, 403: 1, 540: 1, 545: 1, 650: 2,\
728: 1, 739: 1, 767: 1, 786: 1, 807: 3, 820: 1, 825: 1, 874: 1, 893: 2, 900: 1}

    You can then continue using the `FingerprintGenerator` object to
    create other fingerprint types. For example, let's create a 2D
    pharmacophore fingerprint.

    >>> fp = fg.pharm2d_fp()
    >>> print(fp.GetNumOnBits())
    90
    """

    def __init__(self, mol_obj=None):
        self.mol = mol_obj

    @property
    def mol(self):
        """:class:`~luna.wrappers.base.MolWrapper`, \
                :class:`rdkit.Chem.rdchem.Mol`, or \
                :class:`openbabel.pybel.Molecule`: The molecule."""
        return self._mol_obj

    @mol.setter
    def mol(self, mol_obj):
        self._mol_obj = mol_obj

        if mol_obj:
            self._rdmol = MolWrapper(mol_obj).as_rdkit()

    def rdk_fp(self):
        """Generate an RDKit topological fingerprint for the molecule ``mol``.

        Raises
        ------
        FingerprintNotCreated
            If the fingerprint could not be created.
        """
        try:
            return FingerprintMols.FingerprintMol(self._rdmol)
        except Exception as e:
            logger.exception(e)
            raise FingerprintNotCreated("Fingerprint could not be created.")

    def maccs_keys_fp(self):
        """Generate a MACCS keys fingerprint for the molecule ``mol``.

        Raises
        ------
        FingerprintNotCreated
            If the fingerprint could not be created.
        """
        try:
            return MACCSkeys.GenMACCSKeys(self._rdmol)
        except Exception as e:
            logger.exception(e)
            raise FingerprintNotCreated("Fingerprint could not be created.")

    def atom_pairs_fp(self):
        """Generate an atom pairs fingerprint for the molecule ``mol``.

        Raises
        ------
        FingerprintNotCreated
            If the fingerprint could not be created.
        """
        try:
            return Pairs.GetAtomPairFingerprint(self._rdmol)
        except Exception as e:
            logger.exception(e)
            raise FingerprintNotCreated("Fingerprint could not be created.")

    def torsion_fp(self):
        """Generate a topological torsion fingerprint for the molecule ``mol``.

        Raises
        ------
        FingerprintNotCreated
            If the fingerprint could not be created.
        """
        try:
            func = Torsions.GetTopologicalTorsionFingerprintAsIntVect
            return func(self._rdmol)
        except Exception as e:
            logger.exception(e)
            raise FingerprintNotCreated("Fingerprint could not be created.")

    def morgan_fp(self, radius=2, length=2048, features=False, type=1):
        """Generate a Morgan fingerprint for the molecule ``mol``.

        Parameters
        ----------
        radius : int
            Define the maximum radius of the circular neighborhoods
            considered for each atom. The default value is 2, which is
            roughly equivalent to ECFP4 and FCFP4.
        length : int
            The length of the fingerprint. The default value is 2,048.
        features : bool
            If True, use pharmacophoric properties (FCFP) instead of
            atomic invariants (ECFP). The default value is False.
        type : {1, 2, 3}
            Define the type of the Morgan fingerprint function to be
            used, where:

            * ``1`` means **GetMorganFingerprintAsBitVect()**. \
                    It returns an explicit bit vector of size ``length`` \
                    (hashed fingerprint), where 0s and 1s represent the \
                    presence or absence of a given feature, respectively.
            * ``2`` means **GetHashedMorganFingerprint()**. \
                    It returns a sparse int vector ``length`` elements \
                    long (hashed fingerprint) containing the occurrence \
                    number of each feature.
            * ``3`` means **GetMorganFingerprint()**. \
                    It returns a sparse int vector 2^32 elements long \
                    containing the occurrence number of each feature.

            The default value is ``2``.

        Raises
        ------
        FingerprintNotCreated
            If the fingerprint could not be created.
        IllegalArgumentError
            If ``type`` is a value other than 1, 2, or 3.
        """
        try:
            if type == 1:
                func = AllChem.GetMorganFingerprintAsBitVect
                return func(self._rdmol, radius=radius,
                            nBits=length, useFeatures=features)
            elif type == 2:
                func = AllChem.GetHashedMorganFingerprint
                return func(self._rdmol, radius=radius,
                            nBits=length, useFeatures=features)
            elif type == 3:
                func = AllChem.GetMorganFingerprint
                return func(self._rdmol, radius=radius,
                            useFeatures=features)
            else:
                raise IllegalArgumentError("Informed type '%s' is invalid. "
                                           "Available options are 1, 2, or 3."
                                           % str(type))
        except Exception as e:
            logger.exception(e)
            raise FingerprintNotCreated("Fingerprint could not be created.")

    def pharm2d_fp(self, sig_factory=None):
        """Generate a 2D pharmacophore fingerprint for the molecule ``mol``.

        Parameters
        ----------
        sig_factory : RDKit :class:`~rdkit.Chem.Pharm2D.SigFactory`, optional
            Factory object for producing signatures. The default signature
            factory is defined as shown below:

            >>> feat_factory = \
ChemicalFeatures.BuildFeatureFactory(MIN_FDEF_FILE)
            >>> sig_factory = SigFactory(feat_factory,
            ...                          minPointCount=2,
            ...                          maxPointCount=3,
            ...                          trianglePruneBins=False)
            >>> sig_factory.SetBins([(0, 2), (2, 5), (5, 8)])
            >>> sig_factory.Init()

        Raises
        ------
        FingerprintNotCreated
            If the fingerprint could not be created.
        """
        try:
            if sig_factory is None:
                sig_factory = _prepare_pharm2d_fp()["sigFactory"]
            return Generate.Gen2DFingerprint(self._rdmol, sig_factory)
        except Exception as e:
            logger.exception(e)
            raise FingerprintNotCreated("The fingerprint could not "
                                        "be created.")


def available_fp_functions():
    """Return a list of all fingerprints available at
    `FingerprintGenerator`."""
    regex = re.compile(".*([a-zA-Z]+)_fp", flags=0)
    funcs = list(filter(regex.match, dir(FingerprintGenerator)))
    return funcs


def _prepare_morgan_fp(fp_opt):
    params = {}
    if "radius" in fp_opt:
        params["radius"] = fp_opt["radius"]
    if "length" in fp_opt:
        params["length"] = fp_opt["length"]
    if "features" in fp_opt:
        params["features"] = fp_opt["features"]
    if "type" in fp_opt:
        params["type"] = fp_opt["type"]

    return params


def _prepare_pharm2d_fp(fp_opt=None):
    fp_opt = fp_opt or {}

    params = {}
    if "sigFactory" in fp_opt:
        sig_factory = fp_opt["sigFactory"]
    else:
        feat_factory = ChemicalFeatures.BuildFeatureFactory(MIN_FDEF_FILE)
        sig_factory = SigFactory(feat_factory, minPointCount=2,
                                 maxPointCount=3, trianglePruneBins=False)
        sig_factory.SetBins([(0, 2), (2, 5), (5, 8)])
        sig_factory.Init()

    params["sigFactory"] = sig_factory
    return params


def _prepare_fp_params(fp_function, fp_opt):
    if fp_function == "pharm2d_fp" and type(fp_opt) is dict:
        return _prepare_pharm2d_fp(fp_opt)
    elif fp_function == "morgan_fp" and type(fp_opt) is dict:
        return _prepare_morgan_fp(fp_opt)
    else:
        return {}


def generate_fp_for_mols(mols, fp_function=None, fp_opt=None, critical=False):
    """Generate molecular fingerprints for a sequence of molecules.

    Parameters
    ----------
    mols : iterable of :class:`~luna.wrappers.base.MolWrapper`, \
                :class:`rdkit.Chem.rdchem.Mol`, or \
                :class:`openbabel.pybel.Molecule`
        A sequence of molecules.
    fp_function : str
        The fingerprint function to use. The default value is 'pharm2d_fp'.

        To check out the list of available functions, call the function
        :py:meth:`available_fp_functions`.
    fp_opt : dict, optional
        A set of parameters to pass to ``fp_function``.
    critical : bool
        If True, raises any exceptions caught during the generation of
        fingerprints.Otherwise, ignores all exceptions (the default).
        The error messages are always printed to the logging output.

    Returns
    -------
     : list of dict
        A list of dictionaries where each item contains the molecule name
        and its fingerprint.

        The dict is defined as follows:

            * ``mol`` (str): the molecule name;
            * ``fp`` (RDKit \
:class:`~rdkit.DataStructs.cDataStructs.ExplicitBitVect` \
or :class:`~rdkit.DataStructs.cDataStructs.SparseBitVect`): the fingerprint;

    Raises
    ------
    IllegalArgumentError
        If ``fp_function`` is not a function available in
        `FingerprintGenerator`.

    Examples
    --------

    First, let's define a set of molecules.

    >>> from luna.wrappers.base import MolWrapper
    >>> mols = [MolWrapper.from_smiles("N[C@@H](CCC(N)=O)C(O)=O"),
    ...         MolWrapper.from_smiles("C[C@@H](C(=O)O)N"),
    ...         MolWrapper.from_smiles("C1=CC(=CC=C1CC(C(=O)O)N)O")]

    Now, you can generate fingerprints for these molecules using the function
    :meth:`generate_fp_for_mols`. For example, let's create count ECFP4
    fingerprints of size 1,024 for the above molecules.

    >>> from luna.mol.fingerprint import generate_fp_for_mols
    >>> fps = generate_fp_for_mols(mols,
    ...                            fp_function="morgan_fp",
    ...                            fp_opt={"length": 1024})

    Then, you can loop through the results as shown below:

    >>> for d in fps:
    >>>     print(f"{d['mol'].ljust(25)} - \
{len(d['fp'].GetNonzeroElements())}")
    N[C@@H](CCC(N)=O)C(O)=O   - 19
    C[C@@H](C(=O)O)N          - 12
    C1=CC(=CC=C1CC(C(=O)O)N)O - 24
    """

    funcs = available_fp_functions()
    if fp_function not in funcs:
        raise IllegalArgumentError("The fingerprint function is "
                                   "not available.")

    if fp_function is None:
        fp_function = "pharm2d_fp"
        logger.debug("No fingerprint function was defined. So, the default "
                     "fingerprint type will be used: 2D Pharmacophore "
                     "fingerprint.")

    logger.debug("Generating molecular fingerprints for %d molecules."
                 % len(mols))

    params = _prepare_fp_params(fp_function, fp_opt)
    fpg = FingerprintGenerator()

    fp_mols = []
    for idx, mol in enumerate(mols):
        try:
            fpg.mol = mol
            fp = getattr(fpg, fp_function)(**params)
            fp_mols.append({"fp": fp, "mol": mol.GetProp("_Name")})
        except Exception as e:
            logger.error("Molecule at position %d failed. Name: %s"
                         % (idx, mol.GetProp("_Name")))
            logger.exception(e)

            if critical:
                raise

    logger.debug("%d molecular fingerprint(s) created." % len(fp_mols))

    return fp_mols
