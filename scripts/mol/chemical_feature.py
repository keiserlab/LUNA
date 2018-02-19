

class FeatureExtractor():

    def __init__(self, featureFactory=None):
        self.featureFactory = featureFactory

    def set_feature_factory(self, featureFactory):
        self.featureFactory = featureFactory

    def get_atom_features(self, mol):
        features = self.featureFactory.GetFeaturesForMol(mol)

        atomFeatures = {}
        for f in features:
            for atmIdx in f.GetAtomIds():

                atomObj = mol.GetAtomWithIdx(atmIdx)

                if (atomObj.GetPDBResidueInfo() is None):
                    key = atmIdx
                else:
                    key = atomObj.GetPDBResidueInfo().GetSerialNumber()

                if (key in atomFeatures):
                    types = atomFeatures[key]
                else:
                    types = []

                types.append(f.GetFamily())
                atomFeatures[key] = types

        return atomFeatures
