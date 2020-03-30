from collections import defaultdict

# Open Babel
from openbabel.openbabel import OBMol, OBSmartsPattern
from openbabel.pybel import Molecule as PybelMol
# RDKit
from rdkit.Chem import Mol as RDMol

from luna.util import stringcase as case
from luna.util.exceptions import MoleculeObjectTypeError
from luna.mol.wrappers.base import MolWrapper

import logging
logger = logging.getLogger()


class ChemicalFeature():

    def __init__(self, name):
        self.name = name

    def format_name(self, case_func="sentencecase"):
        func = getattr(case, case_func)
        return func(self.name)

    # Special methods
    def __repr__(self):
        return "<Feature=%s>" % self.name

    def __eq__(self, other):
        """Overrides the default implementation"""
        if isinstance(self, other.__class__):
            return self.name == other.name
        return False

    def __ne__(self, other):
        """Overrides the default implementation"""
        return not self.__eq__(other)

    def __lt__(self, other):
        return self.name < other.name

    def __hash__(self):
        return hash(self.name)


class MolChemicalFeature:

    def __init__(self, family, atom_ids):
        self.family = family
        self.atom_ids = atom_ids

    def GetAtomIds(self):
        return self.atom_ids

    def GetFamily(self):
        return self.family


class FeatureExtractor:

    def __init__(self, feature_factory):
        self.feature_factory = feature_factory

    def set_feature_factory(self, feature_factory):
        self.feature_factory = feature_factory

    def get_features_by_atoms(self, mol_obj, atm_map=None):
        if isinstance(mol_obj, RDMol):
            perceived_features = self.feature_factory.GetFeaturesForMol(mol_obj)
        elif isinstance(mol_obj, OBMol):
            perceived_features = self._get_features_from_obmol(mol_obj)
        elif isinstance(mol_obj, PybelMol):
            perceived_features = self._get_features_from_obmol(mol_obj.OBMol)
        else:
            logger.exception("Objects of type '%s' are not currently accepted." % mol_obj.__class__)
            raise MoleculeObjectTypeError("Objects of type '%s' are not currently accepted." % mol_obj.__class__)

        atm_features = defaultdict(set)
        for f in perceived_features:
            for atm_idx in f.GetAtomIds():
                tmp_atm_idx = atm_idx
                if atm_map is not None:
                    if atm_idx in atm_map:
                        tmp_atm_idx = atm_map[atm_idx]
                    else:
                        logger.warning("There is no corresponding mapping to the index '%d'. It will be ignored." % atm_idx)

                feature = ChemicalFeature(f.GetFamily())
                atm_features[tmp_atm_idx].add(feature)

        return atm_features

    def get_features_by_groups(self, mol_obj, atm_map=None):

        if isinstance(mol_obj, MolWrapper):
            mol_obj = mol_obj.unwrap()

        if isinstance(mol_obj, RDMol):
            perceived_features = self.feature_factory.GetFeaturesForMol(mol_obj)
        elif isinstance(mol_obj, OBMol):
            perceived_features = self._get_features_from_obmol(mol_obj)
        elif isinstance(mol_obj, PybelMol):
            perceived_features = self._get_features_from_obmol(mol_obj.OBMol)
        else:
            logger.exception("Objects of type '%s' are not currently accepted." % mol_obj.__class__)
            raise MoleculeObjectTypeError("Objects of type '%s' are not currently accepted." % mol_obj.__class__)

        grp_features = {}
        for f in perceived_features:
            atm_ids = sorted(list(f.GetAtomIds()))

            if atm_map is not None:
                tmp_atm_ids = []
                for atm_id in atm_ids:
                    if atm_id in atm_map:
                        tmp_atm_ids.append(atm_map[atm_id])
                    else:
                        logger.warning("There is no corresponding mapping to the index '%d'. It will be ignored." % atm_id)
                atm_ids = tmp_atm_ids

            key = ','.join([str(x) for x in atm_ids])
            if key in grp_features:
                grp_obj = grp_features[key]
            else:
                grp_obj = {"atm_ids": atm_ids, "features": []}

            features = set(grp_obj["features"])
            feature = ChemicalFeature(f.GetFamily())
            features.add(feature)
            grp_obj["features"] = list(features)
            grp_features[key] = grp_obj

        return grp_features

    def _get_features_from_obmol(self, ob_mol):
        grp_features = defaultdict(set)
        for key, smarts in self.feature_factory.GetFeatureDefs().items():
            grp_type = key.split(".")[0]

            ob_smart = OBSmartsPattern()
            ob_smart.Init(str(smarts))
            ob_smart.Match(ob_mol)

            matches = [x for x in ob_smart.GetMapList()]
            if matches:
                for match in matches:
                    cur_ids = set(match)
                    exists = False
                    remove_ids = []

                    for ids in grp_features[grp_type]:
                        ids = set(ids)
                        # If there is any other group of the same type that already contains the current atoms.
                        if cur_ids.issubset(ids):
                            exists = True
                            # It just need to find one bigger group.
                            break
                        # If the current group contains atoms from already added group atoms.
                        elif ids.issubset(cur_ids):
                            exists = True
                            # Find all smaller groups to remove them.
                            remove_ids.append(tuple(ids))

                    if exists is False:
                        grp_features[grp_type].add(tuple(cur_ids))
                    elif remove_ids:
                        for ids in remove_ids:
                            grp_features[grp_type].remove(ids)
                        grp_features[grp_type].add(tuple(cur_ids))

        return [MolChemicalFeature(family, atom_ids) for family in grp_features for atom_ids in grp_features[family]]
