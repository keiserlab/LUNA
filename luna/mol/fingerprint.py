import re

from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem import AllChem
from rdkit.Chem.Pharm2D import Generate
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Pharm2D.SigFactory import SigFactory

from luna.util.exceptions import FingerprintNotCreated, IllegalArgumentError

import logging

logger = logging.getLogger()


class FingerprintGenerator():

    def __init__(self, mol=None):
        self.molecule = mol

    def set_molecule(self, mol):
        self.molecule = mol

    def topological_fp(self):
        try:
            return FingerprintMols.FingerprintMol(self.molecule)
        except Exception as e:
            logger.exception(e)
            raise FingerprintNotCreated("Fingerprint could not be created.")

    def maccs_keys_fp(self):
        try:
            return MACCSkeys.GenMACCSKeys(self.molecule)
        except Exception as e:
            logger.exception(e)
            raise FingerprintNotCreated("Fingerprint could not be created.")

    def atom_pairs_fp(self):
        try:
            return Pairs.GetAtomPairFingerprint(self.molecule)
        except Exception as e:
            logger.exception(e)
            raise FingerprintNotCreated("Fingerprint could not be created.")

    def torsion_fp(self):
        try:
            return Pairs.GetAtomPairFingerprint(self.molecule)
        except Exception as e:
            logger.exception(e)
            raise FingerprintNotCreated("Fingerprint could not be created.")

    def morgan_fp(self, radius=2, nBits=2048, features=False, type=1):
        """ Create a Morgan fingerprint.

            @param type: defines the type of Morgan fingerprint function
                         to be used.

                         Options:
                            1: GetMorganFingerprintAsBitVect()
                            2: GetHashedMorganFingerprint()
                            3: GetMorganFingerprint
            @type int
        """
        try:
            if (type == 1):
                return AllChem.GetMorganFingerprintAsBitVect(self.molecule, radius=radius, nBits=nBits, useFeatures=features)
            elif (type == 2):
                return AllChem.GetHashedMorganFingerprint(self.molecule, radius=radius, nBits=nBits, useFeatures=features)
            elif (type == 3):
                return AllChem.GetMorganFingerprint(self.molecule, radius=radius, useFeatures=features)
            else:
                raise IllegalArgumentError("Informed type '%s' is invalid. Available options are 1, 2, or 3." % str(type))
        except Exception as e:
            logger.exception(e)
            raise FingerprintNotCreated("Fingerprint could not be created.")

    def pharm2d_fp(self, sigFactory=None):
        try:
            if (sigFactory is None):
                raise IllegalArgumentError("SigFactory object is obligatory.")

            return Generate.Gen2DFingerprint(self.molecule, sigFactory)
        except Exception as e:
            logger.exception(e)
            raise FingerprintNotCreated("Fingerprint could not be created.")


def available_fp_functions():
    regex = re.compile(".*([a-zA-Z]+)_fp", flags=0)
    funcs = list(filter(regex.match, dir(FingerprintGenerator)))
    return funcs


def prepare_morgan_fp(fp_opt):
    params = {}
    if ("radius" in fp_opt):
        params["radius"] = fp_opt["radius"]
    if ("nBits" in fp_opt):
        params["nBits"] = fp_opt["nBits"]
    if ("features" in fp_opt):
        params["features"] = fp_opt["features"]
    if ("type" in fp_opt):
        params["type"] = fp_opt["type"]

    return params


def prepare_pharm2d_fp(fp_opt):
    params = {}
    if ("sigFactory" in fp_opt):
        sigFactory = fp_opt["sigFactory"]
    else:
        fdefName = 'data/MinimalFeatures.fdef'
        featFactory = ChemicalFeatures.BuildFeatureFactory(fdefName)
        sigFactory = SigFactory(featFactory, minPointCount=2, maxPointCount=3, trianglePruneBins=False)
        sigFactory.SetBins([(0, 2), (2, 5), (5, 8)])
        sigFactory.Init()

    params["sigFactory"] = sigFactory
    return params


def prepare_fp_params(fp_function, fp_opt):
    if (fp_function == "pharm2d_fp" and type(fp_opt) is dict):
        return prepare_pharm2d_fp(fp_opt)
    elif (fp_function == "morgan_fp" and type(fp_opt) is dict):
        return prepare_morgan_fp(fp_opt)
    else:
        return {}


def generate_fp_for_mols(mols, fp_function=None, fp_opt=None, critical=False):
    funcs = available_fp_functions()
    if fp_function not in funcs:
        raise IllegalArgumentError("Fingerprint function not available.")

    if fp_function is None:
        fp_function = "pharm2d_fp"
        logger.info("No fingerprint function was defined.")
        logger.info("The default fingerprint type will be used: 2D Pharmacophore fingerprint.")

    logger.info("Trying to generate fingerprints for %d molecules." % len(mols))

    params = prepare_fp_params(fp_function, fp_opt)
    fpg = FingerprintGenerator()

    fp_mols = []
    for idx, mol in enumerate(mols):
        try:
            fpg.set_molecule(mol)
            fp = getattr(fpg, fp_function)(**params)
            fp_mols.append({"fp": fp, "mol": mol.GetProp("_Name")})
        except Exception as e:
            logger.info("Molecule at position %d failed. Name: %s" % (idx, mol.GetProp("_Name")))
            logger.exception(e)

            if critical:
                raise

    logger.info("%d fingerprint(s) created." % len(fp_mols))

    return fp_mols
