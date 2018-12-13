import interaction.util as iu
from mol.chemical_feature import ChemicalFeature

from operator import (le, ge)

from itertools import (chain, combinations, product)
from collections import defaultdict

from util.exceptions import IllegalArgumentError

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

    def __init__(self, inter_conf=DefaultInteractionConf(), funcs=None):
        self.inter_conf = inter_conf

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
                ("HalogenDonor", "Aromatic"): self.calc_xbond_pi,

                ("NegativelyIonizable", "PositivelyIonizable"): self.calc_attractive,
                ("NegativelyIonizable", "NegativelyIonizable"): self.calc_repulsive,
                ("PositivelyIonizable", "PositivelyIonizable"): self.calc_repulsive,
                ("PositivelyIonizable", "Aromatic"): self.calc_cation_pi,
                }

                # TODO: Incluir:
                # Acceptor - XDonor
                # Acceptor - YDonor
                # Anion - pi system
                # Weak donor - acceptor
                # Weak donor - weak acceptor
                # Hydrogen bond with pi system
                # weak donor - pi system


    def calc_inter(self, group1, group2, feat1, feat2):
        func = self.get_function(feat1, feat2)
        if (func is None):
            raise IllegalArgumentError("Does not exist a corresponding"
                                       "function to the features: "
                                       "'%s', '%s'" % (feat1, feat2))

        inter_type_list = func((group1, group2, feat1, feat2))
        if isinstance(inter_type_list, list) is False:
            inter_type_list = []

        return inter_type_list

    def calc_cation_pi(self, params):
        interactions = []

        group1, group2 = params[0:2]
        ccDist = iu.calc_euclidean_distance(group1.centroid,
                                            group2.centroid)

        if self.is_within_boundary(ccDist, "boundary_cutoff", le):
            if (self.is_within_boundary(ccDist,
                                        "max_dist_cation_pi_inter", le)):

                params = {"dist_cation_pi_inter": ccDist}
                iType = InteractionType(group1, group2,
                                        "Cation-pi", params)

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

                # If the angle criteria were not defined, a specific
                # Pi-stacking definition is not possible as it depends
                # on angle criteria. Therefore, a more general classification
                # is used instead, i.e., all interactions will be Pi-stacking.
                if ("min_dihed_ang_pi_pi_inter" not in self.inter_conf.conf and
                        "max_disp_ang_pi_pi_inter" not in self.inter_conf.conf):
                    interType = "Pi-stacking"
                else:
                    if (self.is_within_boundary(dihedralAngle,
                                                "min_dihed_ang_pi_pi_inter", ge)):
                        interType = "Edge-to-face pi-stacking"
                    elif (self.is_within_boundary(dispAngle1,
                                                  "max_disp_ang_pi_pi_inter", le)):
                        interType = "Face-to-face pi-stacking"
                    else:
                        interType = "Parallel-displaced pi-stacking"

                params = {"cc_dist_pi_pi_inter": ccDist,
                          "dihed_ang_pi_pi_inter": dihedralAngle,
                          "disp_ang_pi_pi_inter": dispAngle1}

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

                params = {"dist_hydrop_inter": ccDist}
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

                            params = {"xc_dist_xbond_inter": xaDist,
                                      "disp_ang_xbond_inter": dispAngle,
                                      "cxa_ang_xbond_inter": cxaAngle}

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
        # TODO: adicionar S, P aqui
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

                # TODO: Porque eu estou percorrendo o vetor de Carbonos ligados ao halogênio?
                # Isso não faz sentido. Irá existir somente 1 carbono por vez.
                # TODO: Tentar entender isso.

                # Defining the angle CXA
                for carbonCoord in xCarbonCoords:
                    xcVect = carbonCoord.coord - donorGroup.centroid
                    cxaAngle = iu.calc_angle(xcVect, xaVect)
                    if (self.is_within_boundary(cxaAngle,
                                                "min_cxa_ang_xbond_inter", ge)):
                        validCXAAngles.append(cxaAngle)

                # Interaction model: C-X ---- A-R
                # Defining angle RAX
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
                    params = {"xa_dist_xbond_inter": xaDist,
                              "cxa_ang_xbond_inter": anglePair[0],
                              "xar_ang_xbond_inter": anglePair[1]}

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
                        params = {"da_dist_hb_inter": daDist,
                                  "ha_dist_hb_inter": -1,
                                  "dha_ang_hb_inter": -1}
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

                                params = {"da_dist_hb_inter": daDist,
                                          "ha_dist_hb_inter": haDist,
                                          "dha_ang_hb_inter": dhaAngle}

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

                params = {"dist_attract_inter": ccDist}
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

                params = {"dist_repuls_inter": ccDist}
                iType = InteractionType(group1, group2,
                                        "Repulsive", params)

                interactions.append(iType)
        return interactions

    def is_within_boundary(self, value, key, func):
        if key not in self.inter_conf.conf:
            return True
        else:
            return func(value, self.inter_conf.conf[key])

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

        hbondSet = set()
        attracSet = set()
        h2oPairs = defaultdict(dict)
        dependentInteractions = set()

        # TODO: break code into two functions: Water bridged and SB

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

                hbondSet.add(interaction)
            elif interaction.type == "Attractive":
                attracSet.add(interaction)

        for h2oKey in h2oPairs:
            compPairs = combinations(h2oPairs[h2oKey].keys(), 2)

            for (comp1, comp2) in compPairs:
                compParent1 = comp1.compound
                compParent2 = comp2.compound

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

        # It will try to match Hydrogen bonds and Attractive interactions
        # involving the same chemical groups to attribute salt bridges.
        sbGroups = set()
        for (hbond, attractive) in product(hbondSet, attracSet):

            condA = (attractive.comp1.has_atom(hbond.comp1.atoms[0]) and
                     attractive.comp2.has_atom(hbond.comp2.atoms[0]))

            condB = (attractive.comp1.has_atom(hbond.comp2.atoms[0]) and
                     attractive.comp2.has_atom(hbond.comp1.atoms[0]))

            # If an acceptor atom belongs to a negative group, and the donor
            # to a positive group (and vice-versa), it means that the
            # interaction occurs between the same meioties.
            # However, just one condition should occur. For, example,
            # it is not possible that an acceptor atom belongs to a negative
            # and positive group at the same time.
            if condA ^ condB:

                key1 = (attractive.comp1, attractive.comp2)
                key2 = (attractive.comp2, attractive.comp1)

                if key1 in sbGroups or key2 in sbGroups:
                    continue

                sbGroups.add(key1)
                params = {"dependOf": [hbond, attractive]}
                iType = InteractionType(attractive.comp1, attractive.comp2,
                                        "Salt bridge", params)
                dependentInteractions.add(iType)

        return dependentInteractions


class InteractionFilter:

    def __init__(self, inter_conf=DefaultInteractionConf(),
                 funcs=None, use_inter_criteria=True,
                 use_ignore_tests=True,
                 ignore_same_entity=True,
                 ignore_prot_prot=True,
                 ignore_prot_lig=False,
                 ignore_lig_lig=True,
                 ignore_h2o_pairs=True,
                 positive_list=None):

        self.inter_conf = inter_conf

        if funcs is None:
            funcs = self._default_functions()
        self._funcs = funcs

        self.use_inter_criteria = use_inter_criteria
        self.use_ignore_tests = use_ignore_tests
        self.ignore_same_entity = ignore_same_entity
        self.ignore_prot_prot = ignore_prot_prot
        self.ignore_prot_lig = ignore_prot_lig
        self.ignore_lig_lig = ignore_lig_lig
        self.ignore_h2o_pairs = ignore_h2o_pairs

        if positive_list is None:
            positive_list = []
        self.positive_list = positive_list

    @property
    def funcs(self):
        return self._funcs

    def _default_functions(self):
        return {"Hydrogen bond": self.filter_hbond,
                "Attractive": self.filter_attractive,
                "Cation-pi": self.filter_cation_pi,
                "Edge-to-face pi-stacking": self.filter_pi_pi,
                "Face-to-face pi-stacking": self.filter_pi_pi,
                "Parallel-displaced pi-stacking": self.filter_pi_pi,
                "Pi-stacking": self.filter_pi_pi,
                "Hydrophobic": self.filter_hydrop,
                "Halogen bond": self.filter_xbond,
                "Repulsive": self.filter_repulsive}

    def filter(self, interaction):

        if interaction in self.positive_list:
            return True

        compParent1 = interaction.comp1.compound
        compParent2 = interaction.comp2.compound

        if self.use_ignore_tests:
            isSameEntity = compParent1 == compParent2
            if self.ignore_same_entity and isSameEntity:
                return False

            isProtProt = (compParent1.get_id()[0] == " " and
                          compParent2.get_id()[0] == " ")
            if self.ignore_prot_prot and isProtProt:
                return False

            isProtLig = ((compParent1.get_id()[0] == " " and
                          compParent2.get_id()[0].startswith("H_")) or
                         (compParent1.get_id()[0].startswith("H_") and
                          compParent2.get_id()[0] == " "))
            if self.ignore_prot_lig and isProtLig:
                return False

            isLigLig = (compParent1.get_id()[0].startswith("H_") and
                        compParent2.get_id()[0].startswith("H_"))
            if self.ignore_lig_lig and isLigLig:
                return False

            isH2OPair = (compParent1.get_id()[0] == "W" or
                         compParent2.get_id()[0] == "W")

            if self.ignore_h2o_pairs and isH2OPair:
                return False

        if (self.use_inter_criteria is False or
                interaction.type not in self.funcs):
            return True
        else:
            func = self.funcs[interaction.type]
            return func(interaction)

    def filter_hbond(self, interaction):
        isValid = False

        if self.is_within_boundary(interaction.da_dist_hb_inter,
                                   "max_da_dist_hb_inter", le):

            if (interaction.ha_dist_hb_inter == -1 and
                    interaction.dha_ang_hb_inter == -1):
                # If the angle DHA and the distance HA is equal -1, it is
                # expected that the molecules to be water molecules. "
                if (interaction.comp1.compound.get_id()[0] == "W" or
                        interaction.comp2.compound.get_id()[0] == "W"):

                    ha_dist = interaction.da_dist_hb_inter - 1
                    if (self.is_within_boundary(ha_dist,
                                                "max_ha_dist_hb_inter", le)):
                        isValid = True
            else:
                if (self.is_within_boundary(interaction.ha_dist_hb_inter,
                                            "max_ha_dist_hb_inter", le) and
                    self.is_within_boundary(interaction.dha_ang_hb_inter,
                                            "min_dha_ang_hb_inter", ge)):
                        isValid = True
        return isValid

    def filter_attractive(self, interaction):
        return self.is_within_boundary(interaction.dist_attract_inter,
                                       "max_dist_attract_inter", le)

    def filter_cation_pi(self, interaction):
        return self.is_within_boundary(interaction.dist_cation_pi_inter,
                                       "max_dist_cation_pi_inter", le)

    def filter_pi_pi(self, interaction):
        isValid = False
        if self.is_within_boundary(interaction.cc_dist_pi_pi_inter,
                                   "max_cc_dist_pi_pi_inter", le):

            # If the angle criteria were not defined, it is not possible
            # to test if the interaction fits the requirements.
            # So, it is better to not filter out the interaction.
            if ("min_dihed_ang_pi_pi_inter" not in self.inter_conf.conf and
                    "max_disp_ang_pi_pi_inter" not in self.inter_conf.conf):
                isValid = True
            else:
                if self.is_within_boundary(interaction.dihed_ang_pi_pi_inter,
                                           "min_dihed_ang_pi_pi_inter", ge):
                    # It overwrites a general Pi-stacking interaction.
                    interaction.type = "Edge-to-face pi-stacking"
                    isValid = True
                elif self.is_within_boundary(interaction.disp_ang_pi_pi_inter,
                                             "max_disp_ang_pi_pi_inter", le):
                    # It overwrites a general Pi-stacking interaction.
                    interaction.type = "Face-to-face pi-stacking"
                    isValid = True
                else:
                    # It overwrites a general Pi-stacking interaction.
                    interaction.type = "Parallel-displaced pi-stacking"
                    isValid = True

        return isValid

    def filter_hydrop(self, interaction):
        return self.is_within_boundary(interaction.dist_hydrop_inter,
                                       "max_dist_hydrop_inter", le)

    def filter_xbond(self, interaction):
        # Halogen bond involving a PI system (aromatic ring)
        if (interaction.comp1 == "Aromatic" or
                interaction.comp2 == "Aromatic"):
            return (self.is_within_boundary(interaction.xc_dist_xbond_inter,
                                            "max_xc_dist_xbond_inter", le) and
                    self.is_within_boundary(interaction.disp_ang_xbond_inter,
                                            "max_disp_ang_xbond_inter", le) and
                    self.is_within_boundary(interaction.cxa_ang_xbond_inter,
                                            "min_cxa_ang_xbond_inter", ge))
        # Classical halogen bond, i.e., does not envolve a PI system
        else:
            return (self.is_within_boundary(interaction.xa_dist_xbond_inter,
                                            "max_xa_dist_xbond_inter", le) and
                    self.is_within_boundary(interaction.cxa_ang_xbond_inter,
                                            "min_cxa_ang_xbond_inter", ge) and
                    self.is_within_boundary(interaction.xar_ang_xbond_inter,
                                            "min_xar_ang_xbond_inter", ge))

    def filter_repulsive(self, interaction):
        return self.is_within_boundary(interaction.dist_repuls_inter,
                                       "max_dist_repuls_inter", le)

    def is_within_boundary(self, value, key, func):
        if key not in self.inter_conf.conf:
            return True
        else:
            return func(value, self.inter_conf.conf[key])


def calc_all_interactions(trgt_comp_grps, nb_comp_grps,
                          conf=DefaultInteractionConf()):

    iCalc = InteractionCalculator(conf)

    interactions = []
    for (trgt_comp_group, nb_comp_grp) in product(trgt_comp_grps,
                                                  nb_comp_grps):
        for (trgt_atms_grp, nb_atms_grp) in product(trgt_comp_group.atomGroups,
                                                    nb_comp_grp.atomGroups):

            featurePairs = list(product(trgt_atms_grp.chemicalFeatures,
                                        nb_atms_grp.chemicalFeatures))

            featurePairs = filter(lambda x: iCalc.is_feature_pair_valid(*x),
                                  featurePairs)

            for featPair in featurePairs:
                featTuple = (featPair[0].name, featPair[1].name)

                calc_inter_params = (trgt_atms_grp, nb_atms_grp) + featTuple
                inter_type_list = iCalc.calc_inter(*calc_inter_params)
                interactions.extend(inter_type_list)

    diCalc = DependentInteractions()
    dependent_inter = diCalc.get_dependent_interactions(interactions)
    interactions.extend(dependent_inter)

    keep_inter = list(chain.from_iterable(i.dependOf for i in dependent_inter))
    iFilter = InteractionFilter(positive_list=keep_inter,
                                use_inter_criteria=False)
    filtered_inter = [i for i in interactions if iFilter.filter(i)]

    logger.info("Number of possible interactions: %d" % len(filtered_inter))

    return filtered_inter


def apply_interaction_criteria(interactions, conf=DefaultInteractionConf(),
                               inter_funcs=None):

    # First, it does not use Ignore tests because the interest here is only to
    # apply a filtering based on the defined interaction criteria
    iFilter = InteractionFilter(use_ignore_tests=False,
                                inter_conf=conf,
                                funcs=inter_funcs)
    filtered_inter = [i for i in interactions if iFilter.filter(i)]

    aux_set = set(filtered_inter)
    keep_inter = set()
    remove_inter = set()
    for i in filtered_inter:
        if "dependOf" in i.params:
            remove = False
            for di in i.dependOf:
                if di not in aux_set:
                    remove = True
                    break

            if remove:
                remove_inter.update(set(i.dependOf))
                aux_set.remove(i)
            else:
                keep_inter.update(set(i.dependOf))

    for i in remove_inter:
        if i in aux_set and i not in keep_inter:
            aux_set.remove(i)

    filtered_inter = list(aux_set)

    return filtered_inter
