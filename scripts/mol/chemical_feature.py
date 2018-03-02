import logging
logger = logging.getLogger(__name__)


class ChemicalFeature():

    def __init__(self, name):
        self.name = name

    # Special methods

    def __repr__(self):
        return "<Feature=%s>" % self.name

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        """Overrides the default implementation"""
        if isinstance(self, other.__class__):
            return self.name == other.name
        return False

    def __ne__(self, other):
        """Overrides the default implementation (unnecessary in Python 3)"""
        return not self.__eq__(other)


class FeatureExtractor():

    def __init__(self, featureFactory=None):
        self.featureFactory = featureFactory

    def set_feature_factory(self, featureFactory):
        self.featureFactory = featureFactory

    def get_features_by_atoms(self, rdMol, atomMap=None):
        perceivedFeatures = self.featureFactory.GetFeaturesForMol(rdMol)

        atomFeatures = {}
        for f in perceivedFeatures:
            for atomIdx in f.GetAtomIds():
                tmpAtomIdx = atomIdx
                if (atomMap is not None):
                    if (atomIdx in atomMap):
                        tmpAtomIdx = atomMap[atomIdx]
                    else:
                        logger.warning("Does not exist a corresponding mapping"
                                       " to the index '%d'" % atomIdx)

                if (tmpAtomIdx in atomFeatures):
                    features = set(atomFeatures[tmpAtomIdx])
                else:
                    features = set()

                feature = ChemicalFeature(f.GetFamily())
                features.add(feature)
                atomFeatures[tmpAtomIdx] = features

        return atomFeatures

    def get_features_by_groups(self, rdMol, atomMap=None):
        perceivedFeatures = self.featureFactory.GetFeaturesForMol(rdMol)

        groupFeatures = {}
        for f in perceivedFeatures:
            atomIds = sorted(list(f.GetAtomIds()))

            if (atomMap is not None):
                tmpAtomIds = []
                for i in range(0, len(atomIds)):
                    if (atomIds[i] in atomMap):
                        tmpAtomIds.append(atomMap[atomIds[i]])
                    else:
                        logger.warning("Does not exist a corresponding mapping"
                                       " to the index '%d'. "
                                       "It will be ignored."
                                       % atomIds[i])
                atomIds = tmpAtomIds

            key = ','.join([str(x) for x in atomIds])
            if key in groupFeatures:
                groupObj = groupFeatures[key]
            else:
                groupObj = {"atomIds": atomIds, "features": []}

            features = set(groupObj["features"])
            feature = ChemicalFeature(f.GetFamily())
            features.add(feature)
            groupObj["features"] = list(features)
            groupFeatures[key] = groupObj

        return groupFeatures
