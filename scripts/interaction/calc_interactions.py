import interaction.util as iu
from mol.chemical_feature import ChemicalFeature

from operator import (le, ge)

from itertools import (combinations, product)
from collections import defaultdict

import logging
logger = logging.getLogger(__name__)


class InteractionType():

    def __init__(self, comp1, comp2, type, params=None):
        self.comp1 = comp1
        self.comp2 = comp2
        self.type = type

        if (params is None):
            params = {}
        self._params = params
        self._expand_dict()

    @property
    def params(self):
        return self._params

    def _expand_dict(self):
        for key in self._params:
            self.__dict__[key] = self._params[key]


class InteractionConf():

    def __init__(self, conf):
        self._conf = conf
        self._expand_dict()

    @property
    def conf(self):
        return self._conf

    @property
    def keys(self):
        return [k for k in self._conf]

    def __getattr__(self, key):
        if key in self._conf:
            return self._conf[key]
        else:
            logger.info("Key '%s' does not exist." % key)
            return None

    def _expand_dict(self):
        for key in self._conf:
            self.__dict__[key] = self._conf[key]

    def add(self, key, val):
        if key not in self._conf:
            self._conf[key] = val
            self.__dict__[key] = val
        else:
            logger.info("Key '%s' already exists." % key)

    def alter(self, key, val):
        if key in self._conf:
            self.conf[key] = val
        else:
            logger.info("Key '%s' does not exist." % key)


class DefaultInteractionConf(InteractionConf):

    def __init__(self):

        conf = {}

        # Hydrogen bond
        conf["max_da_dist_hb_inter"] = 3.9
        conf["max_ha_dist_hb_inter"] = 2.5
        conf["min_dha_ang_hb_inter"] = 120

        # Ionic interactions
        conf["max_dist_repuls_inter"] = 6
        conf["max_dist_attract_inter"] = 6

        # Aromatic stacking
        conf["max_cc_dist_pi_pi_inter"] = 6
        conf["min_dihed_ang_pi_pi_inter"] = 30
        conf["max_disp_ang_pi_pi_inter"] = 20

        # Hydrophobic interaction
        conf["max_dist_hydrop_inter"] = 4.5

        # Cation-pi interaction
        conf["max_dist_cation_pi_inter"] = 6

        # Halogen bond.
        # Interaction model: C-X ---- A-R,
        # Where C is a carbon, X a halogen, A an acceptor and
        # R is an atom bonded to A.
        # Distance X-A when A is an single atom.
        conf["max_xa_dist_xbond_inter"] = 4
        # Distance X-A when A is an aromatic ring, so C comes from Centroid.
        conf["max_xc_dist_xbond_inter"] = 4.5
        # Ref: Halogen bonds in biological molecules [Auffinger, 2004]
        conf["min_cxa_ang_xbond_inter"] = 120
        conf["min_xar_ang_xbond_inter"] = 80
        conf["max_disp_ang_xbond_inter"] = 60

        # Covalent interactions
        conf["cov_dist_tolerance"] = 0.45

        conf["boundary_cutoff"] = 7

        super().__init__(conf)


class InteractionCalculator():

    def __init__(self, interConf=None, funcs=None):
        if (interConf is None):
            interConf = {}

        self.interConf = interConf

        if (funcs is None):
            funcs = self._default_functions()
        self._funcs = funcs

    @property
    def funcs(self):
        return self._funcs

    def _default_functions(self):
        return {("Hydrophobic", "Hydrophobic"): self.calc_hydrop,
                ("Donor", "Acceptor"): self.calc_hbond,
                ("Aromatic", "Aromatic"): self.calc_pi_pi,
                ("Negative", "Positive"): self.calc_attractive,
                ("Negative", "Negative"): self.calc_repulsive,
                ("Positive", "Positive"): self.calc_repulsive,
                ("Hydrophobe", "Hydrophobe"): self.calc_hydrop,
                ("NegIonizable", "PosIonizable"): self.calc_attractive,
                ("NegIonizable", "NegIonizable"): self.calc_repulsive,
                ("PosIonizable", "PosIonizable"): self.calc_repulsive,
                ("PosIonizable", "Aromatic"): self.calc_cation_pi,
                ("HalogenDonor", "HalogenAcceptor"): self.calc_xbond,
                ("HalogenDonor", "Aromatic"): self.calc_xbond_pi}

    def calc_inter(self, group1, group2, feat1, feat2):
        from util.exceptions import IllegalArgumentError

        func = self.get_function(feat1, feat2)
        if (func is None):
            raise IllegalArgumentError("Does not exist a corresponding"
                                       "function to the features: "
                                       "'%s', '%s'" % (feat1, feat2))

        iTypeList = func((group1, group2, feat1, feat2))
        if isinstance(iTypeList, list) is False:
            iTypeList = []

        return iTypeList

    def calc_cation_pi(self, params):
        interactions = []

        group1, group2 = params[0:2]
        ccDist = iu.calc_euclidean_distance(group1.centroid,
                                            group2.centroid)

        if self.is_within_boundary(ccDist, "boundary_cutoff", le):
            if (self.is_within_boundary(ccDist,
                                        "max_dist_cation_pi_inter", le)):

                params = {"cc_dist": ccDist}
                iType = InteractionType(group1, group2,
                                        "Cation-pi interaction", params)

                interactions.append(iType)
        return interactions

    def calc_pi_pi(self, params):
        interactions = []

        ring1, ring2 = params[0:2]
        ccDist = iu.calc_euclidean_distance(ring1.centroid, ring2.centroid)

        if self.is_within_boundary(ccDist, "boundary_cutoff", le):
            if (self.is_within_boundary(ccDist,
                                        "max_cc_dist_pi_pi_inter", le)):

                dihedralAngle = iu.to_quad1(iu.calc_angle(ring1.normal,
                                                          ring2.normal))
                vectorCC = ring2.centroid - ring1.centroid
                dispAngle1 = iu.to_quad1(iu.calc_angle(ring1.normal, vectorCC))

                # If the angle criterias was not defined, a specific
                # Pi-stacking definition is not possible as it depends
                # on angle criterias. Therefore, a more general classification
                # is used instead, i.e., all interactions will be Pi-stacking.
                if ("min_dihed_ang_pi_pi_inter" not in self.interConf.conf and
                        "max_disp_ang_pi_pi_inter" not in self.interConf.conf):
                    interType = "Pi-stacking"
                else:
                    if (self.is_within_boundary(dihedralAngle,
                                                "min_dihed_ang_pi_pi_inter", ge)):
                        interType = "Edge-to-face pi-stacking"
                    else:
                        if (self.is_within_boundary(dispAngle1,
                                                    "max_disp_ang_pi_pi_inter", le)):
                            interType = "Face-to-face pi-stacking"
                        else:
                            interType = "Parallel-displaced pi-stacking"

                params = {"cc_dist": ccDist,
                          "dihed_angle": dihedralAngle,
                          "disp_angle": dispAngle1}
                iType = InteractionType(ring1, ring2,
                                        interType, params)

                interactions.append(iType)
        return interactions

    def calc_hydrop(self, params):
        interactions = []

        group1, group2 = params[0:2]
        ccDist = iu.calc_euclidean_distance(group1.centroid,
                                            group2.centroid)

        if self.is_within_boundary(ccDist, "boundary_cutoff", le):
            if (self.is_within_boundary(ccDist,
                                        "max_dist_hydrop_inter", le)):

                params = {"cc_dist": ccDist}
                iType = InteractionType(group1, group2,
                                        "Hydrophobic", params)

                interactions.append(iType)
        return interactions

    def calc_xbond_pi(self, params):
        interactions = []
        group1, group2, feat1, feat2 = params

        if (feat1 == "Aromatic" and feat2 == "HalogenDonor"):
            donorGroup = group2
            ringGroup = group1
        else:
            donorGroup = group1
            ringGroup = group2

        # There are always just one donor/acceptor atom.
        donorAtom = donorGroup.atoms[0]
        # Recover only carbon coordinates
        xCarbonCoords = [x for x in donorAtom.nbCoords.coords
                         if x.atomicNumber == 6]

        # Interaction model: C-X ---- A-R
        # Defining the XA distance, in which A is the ring center
        xaDist = iu.calc_euclidean_distance(donorGroup.centroid,
                                            ringGroup.centroid)

        if self.is_within_boundary(xaDist, "boundary_cutoff", le):
            if (self.is_within_boundary(xaDist,
                                        "max_xc_dist_xbond_inter", le)):

                axVect = donorGroup.centroid - ringGroup.centroid
                dispAngle = iu.to_quad1(iu.calc_angle(ringGroup.normal,
                                                      axVect))

                if (self.is_within_boundary(dispAngle,
                                            "max_disp_ang_xbond_inter", le)):
                    # Interaction model: C-X ---- A-R
                    # XA vector is always the same
                    xaVect = ringGroup.centroid - donorGroup.centroid
                    # Defining angle CXA, in which A is the ring center
                    for carbonCoord in xCarbonCoords:
                        xcVect = carbonCoord.coord - donorGroup.centroid
                        cxaAngle = iu.calc_angle(xcVect, xaVect)

                        if (self.is_within_boundary(cxaAngle,
                                                    "min_cxa_ang_xbond_inter", ge)):
                            params = {}
                            params["xa_dist"] = xaDist
                            params["cxa_angle"] = cxaAngle
                            params["disp_angle"] = dispAngle

                            iType = InteractionType(donorGroup, ringGroup,
                                                    "Halogen bond", params)

                            interactions.append(iType)
        return interactions

    def calc_xbond(self, params):

        interactions = []
        group1, group2, feat1, feat2 = params

        if (feat1 == "HalogenAcceptor" and feat2 == "HalogenDonor"):
            donorGroup = group2
            acceptorGroup = group1
        else:
            donorGroup = group1
            acceptorGroup = group2

        # There are always just one donor/acceptor atom.
        donorAtom = donorGroup.atoms[0]
        acceptorAtom = acceptorGroup.atoms[0]

        # Recover only carbon coordinates
        xCarbonCoords = [x for x in donorAtom.nbCoords.coords
                         if x.atomicNumber == 6]

        # Interaction model: C-X ---- A-R
        # Distance XY.
        xaDist = iu.calc_euclidean_distance(donorGroup.centroid,
                                            acceptorGroup.centroid)

        if self.is_within_boundary(xaDist, "boundary_cutoff", le):
            if (self.is_within_boundary(xaDist,
                                        "max_xa_dist_xbond_inter", le)):

                # Interaction model: C-X ---- A-R
                validCXAAngles = []
                # XA vector is always the same
                xaVect = acceptorGroup.centroid - donorGroup.centroid
                # Defining angle CXY
                for carbonCoord in xCarbonCoords:
                    xcVect = carbonCoord.coord - donorGroup.centroid
                    cxaAngle = iu.calc_angle(xcVect, xaVect)
                    if (self.is_within_boundary(cxaAngle,
                                                "min_cxa_ang_xbond_inter", ge)):
                        validCXAAngles.append(cxaAngle)

                # Interaction model: C-X ---- A-R
                # Defining angle RYX
                validXARAngles = []
                # AX vector is always the same
                axVect = donorGroup.centroid - acceptorGroup.centroid
                for nbAtomCoord in acceptorAtom.nbCoords.coords:
                    arVect = nbAtomCoord.coord - acceptorGroup.centroid
                    xarAngle = iu.calc_angle(axVect, arVect)

                    if (self.is_within_boundary(xarAngle,
                                                "min_xar_ang_xbond_inter", ge)):
                        validXARAngles.append(xarAngle)

                for anglePair in product(validCXAAngles, validXARAngles):
                    params = {}
                    params["xa_dist"] = xaDist
                    params["cxa_angle"] = anglePair[0]
                    params["xar_angle"] = anglePair[1]

                    iType = InteractionType(donorGroup, acceptorGroup,
                                            "Halogen bond", params)

                    interactions.append(iType)
        return interactions

    def calc_hbond(self, params):
        interactions = []

        group1, group2, feat1, feat2 = params

        if (feat1 == "Acceptor" and feat2 == "Donor"):
            donorGroup = group2
            acceptorGroup = group1
        else:
            donorGroup = group1
            acceptorGroup = group2

        donorAtom = donorGroup.atoms[0]

        daDist = iu.calc_euclidean_distance(donorGroup.centroid,
                                            acceptorGroup.centroid)

        if self.is_within_boundary(daDist, "boundary_cutoff", le):
            if self.is_within_boundary(daDist, "max_da_dist_hb_inter", le):

                if donorAtom.get_parent().get_id()[0] == "W":
                    haDist = daDist - 1

                    if (self.is_within_boundary(haDist,
                                                "max_ha_dist_hb_inter", le)):
                        params = {}
                        params["da_dist"] = daDist
                        params["ha_dist"] = -1
                        params["dha_angle"] = -1

                        iType = InteractionType(donorGroup, acceptorGroup,
                                                "Hydrogen bond", params)

                        interactions.append(iType)
                else:
                    # Recover only hydrogen coordinates
                    hydrogCoords = [x for x in donorAtom.nbCoords.coords
                                    if x.atomicNumber == 1]

                    for hydrogCoord in hydrogCoords:
                            haDist = iu.calc_euclidean_distance(hydrogCoord.coord,
                                                                acceptorGroup.centroid)

                            dhVect = donorGroup.centroid - hydrogCoord.coord
                            haVect = acceptorGroup.centroid - hydrogCoord.coord
                            dhaAngle = iu.calc_angle(haVect, dhVect)

                            if (self.is_within_boundary(haDist,
                                                        "max_ha_dist_hb_inter", le) and
                                self.is_within_boundary(dhaAngle,
                                                        "min_dha_ang_hb_inter", ge)):

                                params = {}
                                params["da_dist"] = daDist
                                params["ha_dist"] = haDist
                                params["dha_angle"] = dhaAngle

                                iType = InteractionType(donorGroup, acceptorGroup,
                                                        "Hydrogen bond", params)

                                interactions.append(iType)
        return interactions

    def calc_attractive(self, params):
        interactions = []

        group1, group2 = params[0:2]
        ccDist = iu.calc_euclidean_distance(group1.centroid,
                                            group2.centroid)

        if self.is_within_boundary(ccDist, "boundary_cutoff", le):
            if (self.is_within_boundary(ccDist,
                                        "max_dist_attract_inter", le)):

                params = {}
                params["cc_dist"] = ccDist

                iType = InteractionType(group1, group2,
                                        "Attractive", params)

                interactions.append(iType)
        return interactions

    def calc_repulsive(self, params):
        interactions = []

        group1, group2 = params[0:2]
        ccDist = iu.calc_euclidean_distance(group1.centroid,
                                            group2.centroid)

        if self.is_within_boundary(ccDist, "boundary_cutoff", le):
            if (self.is_within_boundary(ccDist,
                                        "max_dist_repuls_inter", le)):

                params = {}
                params["cc_dist"] = ccDist

                iType = InteractionType(group1, group2,
                                        "Repulsive", params)

                interactions.append(iType)
        return interactions

    def is_within_boundary(self, value, key, func):
        if key not in self.interConf.conf:
            return True
        else:
            return func(value, self.interConf.conf[key])

    def is_feature_pair_valid(self, feat1, feat2):
        if isinstance(feat1, ChemicalFeature):
            feat1 = feat1.name
        if isinstance(feat2, ChemicalFeature):
            feat2 = feat2.name

        funcs = self._funcs
        return (True if ((feat1, feat2) in funcs or
                         (feat2, feat1) in funcs) else False)

    def get_function(self, feat1, feat2):
        if isinstance(feat1, ChemicalFeature):
            feat1 = feat1.name
        if isinstance(feat2, ChemicalFeature):
            feat2 = feat2.name

        funcs = self._funcs
        if (feat1, feat2) in funcs:
            return funcs[(feat1, feat2)]
        elif (feat2, feat1) in funcs:
            return funcs[(feat2, feat1)]
        else:
            return None

    def set_function(self, pair, func):
        self._funcs[pair] = func


class DependentInteractions:

    def __init__(self,
                 IGNORE_SAME_ENTITY=True,
                 IGNORE_PROT_PROT=True,
                 IGNORE_PROT_LIG=False,
                 IGNORE_LIG_LIG=True,
                 IGNORE_H2O_PAIRS=True):

        self.IGNORE_SAME_ENTITY = IGNORE_SAME_ENTITY
        self.IGNORE_PROT_PROT = IGNORE_PROT_PROT
        self.IGNORE_PROT_LIG = IGNORE_PROT_LIG
        self.IGNORE_LIG_LIG = IGNORE_LIG_LIG
        self.IGNORE_H2O_PAIRS = IGNORE_H2O_PAIRS

    def get_dependent_interactions(self, interactions):

        dependentInteractions = set()

        h2oPairs = defaultdict(dict)

        for interaction in interactions:
            if interaction.type == "Hydrogen bond":
                compParent1 = interaction.comp1.atoms[0].get_parent()
                compParent2 = interaction.comp2.atoms[0].get_parent()

                if compParent1.get_id()[0] == "W":
                    h2oKey = interaction.comp1
                    secondKey = interaction.comp2

                    if secondKey not in h2oPairs[h2oKey]:
                        h2oPairs[h2oKey][secondKey] = interaction

                if compParent2.get_id()[0] == "W":
                    h2oKey = interaction.comp2
                    secondKey = interaction.comp1

                    if secondKey not in h2oPairs[h2oKey]:
                        h2oPairs[h2oKey][secondKey] = interaction

        for h2oKey in h2oPairs:

            compPairs = combinations(h2oPairs[h2oKey].keys(), 2)

            for (comp1, comp2) in compPairs:
                compParent1 = comp1.atoms[0].get_parent()
                compParent2 = comp2.atoms[0].get_parent()

                isSameEntity = compParent1 == compParent2
                if (self.IGNORE_SAME_ENTITY and isSameEntity):
                    continue

                isProtProt = (compParent1.get_id()[0] == " " and
                              compParent2.get_id()[0] == " ")
                if (self.IGNORE_PROT_PROT and isProtProt):
                    continue

                isProtLig = ((compParent1.get_id()[0] == " " and
                              compParent2.get_id()[0].startswith("H_")) or
                             (compParent1.get_id()[0].startswith("H_") and
                              compParent2.get_id()[0] == " "))
                if (self.IGNORE_PROT_LIG and isProtLig):
                    continue

                isLigLig = (compParent1.get_id()[0].startswith("H_") and
                            compParent2.get_id()[0].startswith("H_"))
                if (self.IGNORE_LIG_LIG and isLigLig):
                    continue

                isH2OPair = (compParent1.get_id()[0] == "W" or
                             compParent2.get_id()[0] == "W")
                if (self.IGNORE_H2O_PAIRS and isH2OPair):
                    continue

                params = {"dependOf": [h2oPairs[h2oKey][comp1],
                                       h2oPairs[h2oKey][comp2]]}

                iType = InteractionType(comp1, comp2,
                                        "Water-bridged hydrogen bond", params)

                dependentInteractions.add(iType)

                # print(comp1.chemicalFeatures, comp2.chemicalFeatures)
                # print(h2oKey.atoms[0].get_parent())
                # print(compParent1, compParent2)
                # print(comp1, comp2)
                # print()

        print(len(dependentInteractions))
        exit()




def calc_all_interactions(targetCompGroups, nearbyCompGroups,
                          conf=DefaultInteractionConf()):

    iCalc = InteractionCalculator(conf)

    interactions = []

    for (targetCompGroup, nearbyCompGroup) in product(targetCompGroups,
                                                      nearbyCompGroups):

        for (targetAtmsGrp, nbAtmsGrp) in product(targetCompGroup.atomGroups,
                                                  nearbyCompGroup.atomGroups):

            featurePairs = list(product(targetAtmsGrp.chemicalFeatures,
                                        nbAtmsGrp.chemicalFeatures))

            featurePairs = filter(lambda x: iCalc.is_feature_pair_valid(*x),
                                  featurePairs)

            for featPair in featurePairs:
                featTuple = (featPair[0].name, featPair[1].name)

                calcInterParams = (targetAtmsGrp, nbAtmsGrp) + featTuple
                iTypeList = iCalc.calc_inter(*calcInterParams)
                interactions.extend(iTypeList)

    diCalc = DependentInteractions()
    diCalc.get_dependent_interactions(interactions)

    exit()

    print(len(interactions))
    from data.summary import count_interaction_types
    count = count_interaction_types(interactions)
    print(count)

    logger.info("Number of possible interactions: %d" % len(interactions))

    exit()

def filter_interactions():
    pass