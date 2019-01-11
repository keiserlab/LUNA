import util.stringcase as case

import logging
logger = logging.getLogger(__name__)


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

    def __hash__(self):
        return hash(self.name)


class FeatureExtractor():

    def __init__(self, feature_factory=None):
        self.feature_factory = feature_factory

    def set_feature_factory(self, feature_factory):
        self.feature_factory = feature_factory

    def get_features_by_atoms(self, rdMol, atm_map=None):
        perceived_features = self.feature_factory.GetFeaturesForMol(rdMol)

        atm_features = {}
        for f in perceived_features:
            for atm_idx in f.GetAtomIds():
                tmp_atm_idx = atm_idx
                if (atm_map is not None):
                    if (atm_idx in atm_map):
                        tmp_atm_idx = atm_map[atm_idx]
                    else:
                        logger.warning("Does not exist a corresponding mapping"
                                       " to the index '%d'" % atm_idx)

                if (tmp_atm_idx in atm_features):
                    features = set(atm_features[tmp_atm_idx])
                else:
                    features = set()

                feature = ChemicalFeature(f.GetFamily())
                features.add(feature)
                atm_features[tmp_atm_idx] = features

        return atm_features

    def get_features_by_groups(self, rdMol, atm_map=None):
        perceived_features = self.feature_factory.GetFeaturesForMol(rdMol)

        grp_features = {}
        for f in perceived_features:
            atm_ids = sorted(list(f.GetAtomIds()))

            if atm_map is not None:
                tmp_atm_ids = []
                for i in range(0, len(atm_ids)):
                    if (atm_ids[i] in atm_map):
                        tmp_atm_ids.append(atm_map[atm_ids[i]])
                    else:
                        logger.warning("Does not exist a corresponding mapping to the index '%d'. It will be ignored."
                                       % atm_ids[i])
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
