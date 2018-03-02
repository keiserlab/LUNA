import interaction.util as iu
import interaction.plane_regression as plane

from mol.chemical_feature import ChemicalFeature

from operator import (le, ge)

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
        conf["min_dha_ang_hb_inter"] = 120
        conf["max_da_dist_hb_inter"] = 3.9
        conf["max_ha_dist_hb_inter"] = 2.5

        # Ionic interactions
        conf["max_dist_repuls_inter"] = 6
        conf["max_dist_attract_inter"] = 6

        # Aromatic stacking
        conf["max_cc_dist_pi_pi_inter"] = 6
        conf["min_dihed_ang_pi_pi_inter"] = 30
        conf["max_disp_ang_pi_pi_inter"] = 20

        # Hydrophobic interaction
        conf["max_hydrop_inter"] = 4.5

        # Cation-pi interaction
        conf["max_cation_pi_inter"] = 6

        # Halogen bond.
        # Interaction model: C-X ---- A-R,
        # Where C is a carbon, X a halogen, A an acceptor and R is an atom bonded to A.
        # Distance X-A when A is an single atom.
        conf["max_xa_dist_xbond_inter"] = 4
        # Distance X-A when A is an aromatic ring, so C comes from Centroid.
        conf["max_xc_dist_xbond_inter"] = 4.5
        # Ref: Halogen bonds in biological molecules [Auffinger, 2004]
        conf["min_cxa_angle_xbond_inter"] = 120
        conf["min_xar_angle_xbond_inter"] = 80
        conf["max_disp_angle_xbond_inter"] = 60

        # Covalent interactions
        conf["cov_dist_tolerance"] = 0.45

        conf["boundary_cutoff"] = 7

        super().__init__(conf)


def get_feature_tuple(pair):
    return (pair[0].name, pair[1].name)


def calc_all_interactions(targetCompGroups, nearbyCompGroups):
    from itertools import product

    iConf = DefaultInteractionConf()
    conf = {"boundary_cutoff": iConf.boundary_cutoff}
    # iConf = InteractionConf(conf)

    iCalc = InteractionCalculator(iConf)

    interactions = []
    for compGroups in nearbyCompGroups:
        for groupPair in product(targetCompGroups.atomGroups,
                                 compGroups.atomGroups):

            ligAtomGroups, compAtomGroups = groupPair

            featurePairs = list(product(ligAtomGroups.chemicalFeatures,
                                        compAtomGroups.chemicalFeatures))

            featurePairs = filter(lambda x: iCalc.is_feature_pair_valid(*x),
                                  featurePairs)

            for featPair in featurePairs:

                featTuple = get_feature_tuple(featPair)

                iTypeList = iCalc.calc_inter(*groupPair, *featTuple)

                if (len(iTypeList) != 0):
                    print("############################")
                    print(ligAtomGroups.atoms[0].get_parent())
                    print(compAtomGroups.atoms[0].get_parent())
                    print(".....")
                    print(ligAtomGroups.atoms)
                    print(compAtomGroups.atoms)
                    print(featTuple)

                    for iType in iTypeList:
                        print(iType.type)
                        print(iType.params)

                    interactions.extend(iTypeList)
                    print()
                    print()

    logger.info("Number of possible interactions: %d" % len(interactions))


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
                                        "max_cation_pi_inter", le)):

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
                                        "max_hydrop_inter", le)):

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
                                            "max_disp_angle_xbond_inter", le)):
                    # Interaction model: C-X ---- A-R
                    # XA vector is always the same
                    xaVect = ringGroup.centroid - donorGroup.centroid
                    # Defining angle CXA, in which A is the ring center
                    for carbonCoord in xCarbonCoords:
                        xcVect = carbonCoord.coord - donorGroup.centroid
                        cxaAngle = iu.calc_angle(xcVect, xaVect)

                        if (self.is_within_boundary(cxaAngle,
                                                    "min_cxa_angle_xbond_inter", ge)):
                            params = {}
                            params["xa_dist"] = xaDist
                            params["cxa_angle"] = cxaAngle
                            params["disp_angle"] = dispAngle

                            iType = InteractionType(donorGroup, ringGroup,
                                                    "Halogen bond", params)

                            interactions.append(iType)
        return interactions

    def calc_xbond(self, params):
        from itertools import product

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
                                                "min_cxa_angle_xbond_inter", ge)):
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
                                                "min_xar_angle_xbond_inter", ge)):
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

        # Recover only hydrogen coordinates
        hydrogCoords = [x for x in donorAtom.nbCoords.coords
                        if x.atomicNumber == 1]

        for hydrogCoord in hydrogCoords:
            daDist = iu.calc_euclidean_distance(donorGroup.centroid,
                                                acceptorGroup.centroid)

            if self.is_within_boundary(daDist, "boundary_cutoff", le):

                haDist = iu.calc_euclidean_distance(hydrogCoord.coord,
                                                    acceptorGroup.centroid)

                # dhVect = hydrogCoord.coord - donorGroup.centroid
                # haVect = acceptorGroup.centroid - hydrogCoord.coord
                # dhaAngle = iu.calc_angle(dhVect, haVect)

                dhVect = donorGroup.centroid - hydrogCoord.coord
                haVect = acceptorGroup.centroid - hydrogCoord.coord
                dhaAngle = iu.calc_angle(haVect, dhVect)

                if (self.is_within_boundary(daDist,
                                            "max_da_dist_hb_inter", le) and
                    self.is_within_boundary(haDist,
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
        if (key not in self.interConf.conf):
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
