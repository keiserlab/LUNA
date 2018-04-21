from util.exceptions import FingerprintNotCreated

import logging
logger = logging.getLogger(__name__)


class FingerprintGenerator():

    def __init__(self, mol=None):
        self.molecule = mol

    def set_molecule(self, mol):
        self.molecule = mol

    def topological_fp(self):
        from rdkit.Chem.Fingerprints import FingerprintMols

        try:
            return FingerprintMols.FingerprintMol(self.molecule)
        except Exception as e:
            logger.exception(e)
            raise FingerprintNotCreated("Fingerprint could not be created.")

    def maccs_keys_fp(self):
        from rdkit.Chem import MACCSkeys

        try:
            return MACCSkeys.GenMACCSKeys(self.molecule)
        except Exception as e:
            logger.exception(e)
            raise FingerprintNotCreated("Fingerprint could not be created.")

    def atom_pairs_fp(self):
        from rdkit.Chem.AtomPairs import Pairs

        try:
            return Pairs.GetAtomPairFingerprint(self.molecule)
        except Exception as e:
            logger.exception(e)
            raise FingerprintNotCreated("Fingerprint could not be created.")

    def torsion_fp(self):
        from rdkit.Chem.AtomPairs import Pairs

        try:
            return Pairs.GetAtomPairFingerprint(self.molecule)
        except Exception as e:
            logger.exception(e)
            raise FingerprintNotCreated("Fingerprint could not be created.")

    def morgan_fp(self, radius=2, nBits=2048, features=True, type=1):
        """ Create a Morgan fingerprint.

            @param type: defines the type of Morgan fingerprint function
                         to be used.

                         Options:
                            1: GetMorganFingerprintAsBitVect()
                            2: GetHashedMorganFingerprint()
                            3: GetMorganFingerprint
            @type int
        """
        from rdkit.Chem import AllChem
        from util.exceptions import IllegalArgumentError

        try:
            if (type == 1):
                return AllChem.GetMorganFingerprintAsBitVect(self.molecule,
                                                             radius=radius,
                                                             nBits=nBits,
                                                             useFeatures=features)
            elif (type == 2):
                return AllChem.GetHashedMorganFingerprint(self.molecule,
                                                          radius=radius,
                                                          nBits=nBits,
                                                          useFeatures=features)
            elif (type == 3):
                return AllChem.GetMorganFingerprint(self.molecule,
                                                    radius=radius,
                                                    useFeatures=features)
            else:
                raise IllegalArgumentError("Informed type '%s' is invalid. "
                                           "Available options are 1, 2, or 3."
                                           % str(type))
        except Exception as e:
            logger.exception(e)
            raise FingerprintNotCreated("Fingerprint could not be created.")

    def pharm2d_fp(self, sigFactory=None):

        from rdkit.Chem.Pharm2D import Generate
        from util.exceptions import IllegalArgumentError

        try:
            if (sigFactory is None):
                raise IllegalArgumentError("SigFactory object is obligatory.")

            return Generate.Gen2DFingerprint(self.molecule, sigFactory)
        except Exception as e:
            logger.exception(e)
            raise FingerprintNotCreated("Fingerprint could not be created.")


def available_fp_functions():
    import re

    regex = re.compile(".*([a-zA-Z]+)_fp", flags=0)
    funcs = list(filter(regex.match, dir(FingerprintGenerator)))

    return funcs


def prepare_morgan_fp(fpOpt):
    params = {}
    if ("radius" in fpOpt):
        params["radius"] = fpOpt["radius"]
    if ("nBits" in fpOpt):
        params["nBits"] = fpOpt["nBits"]
    if ("features" in fpOpt):
        params["features"] = fpOpt["features"]
    if ("type" in fpOpt):
        params["type"] = fpOpt["type"]

    return params


def prepare_pharm2d_fp(fpOpt):
    from rdkit.Chem import ChemicalFeatures
    from rdkit.Chem.Pharm2D.SigFactory import SigFactory

    params = {}
    if ("sigFactory" in fpOpt):
        sigFactory = fpOpt["sigFactory"]
    else:
        fdefName = 'data/MinimalFeatures.fdef'
        featFactory = ChemicalFeatures.BuildFeatureFactory(fdefName)
        sigFactory = SigFactory(featFactory, minPointCount=2, maxPointCount=3,
                                trianglePruneBins=False)
        sigFactory.SetBins([(0, 2), (2, 5), (5, 8)])
        sigFactory.Init()

    params["sigFactory"] = sigFactory
    return params


def prepare_fp_params(fpFunction, fpOpt):

    if (fpFunction == "pharm2d_fp" and type(fpOpt) is dict):
        params = prepare_pharm2d_fp(fpOpt)
    elif (fpFunction == "morgan_fp" and type(fpOpt) is dict):
        params = prepare_morgan_fp(fpOpt)
    else:
        params = {}

    return params


def generate_fp_for_mols(mols, fpFunction=None, fpOpt=None, critical=False):
    from util.exceptions import IllegalArgumentError
    import logging
    logger = logging.getLogger(__name__)

    funcs = available_fp_functions()
    if (fpFunction not in funcs):
        raise IllegalArgumentError("Fingerprint function not available.")

    if (fpFunction is None):
        fpFunction = "pharm2d_fp"

        logger.info("No fingerprint function was defined.")
        logger.info("The default fingerprint type will be used: 2D "
                    "Pharmacophore fingerprint.")

    logger.info("Trying to generate fingerprints for %d molecules."
                % len(mols))

    params = prepare_fp_params(fpFunction, fpOpt)
    fpg = FingerprintGenerator()

    fpMols = []
    for idx, mol in enumerate(mols):
        try:
            fpg.set_molecule(mol)
            fp = getattr(fpg, fpFunction)(**params)
            fpMols.append({"fp": fp, "mol": mol.GetProp("_Name")})
        except Exception as e:
            logger.info("Molecule at position %d failed. Name: %s" %
                        (idx, mol.GetProp("_Name")))
            logger.exception(e)

            if (critical):
                raise

    logger.info("%d fingerprint(s) created." % len(fpMols))

    return fpMols
