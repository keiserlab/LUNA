import mol.interaction.math as iu

from mol.interaction.conf import DefaultInteractionConf
from mol.interaction.filter import InteractionFilter
from mol.interaction.type import InteractionType
from mol.features import ChemicalFeature

from operator import (le, ge)
from itertools import (chain, combinations, product)
from collections import defaultdict
from util.exceptions import IllegalArgumentError

import logging
logger = logging.getLogger(__name__)


class InteractionCalculator:

    def __init__(self, inter_conf=DefaultInteractionConf(), inter_filter=None,
                 inter_funcs=None, add_proximal=False, add_cov_inter=False,
                 add_clash=False, add_orphan_h2o_pair=False):

        self.inter_conf = inter_conf

        self.add_proximal = add_proximal
        self.add_cov_inter = add_cov_inter
        self.add_clash = add_clash
        self.add_orphan_h2o_pair = add_orphan_h2o_pair

        if inter_filter is None:
            inter_filter = InteractionFilter(inter_conf=inter_conf)
        self.inter_filter = inter_filter

        if inter_funcs is None:
            inter_funcs = self._default_functions()

        self._inter_funcs = inter_funcs

    @property
    def funcs(self):
        return self._inter_funcs

    def calc_interactions(self, trgt_comp_grps, nb_comp_grps=None):
        all_interactions = []

        # If nb_comp_grps was not informed, it uses the trgt_comp_grps as the neighbors.
        # In this case, the interactions will be target x target.
        nb_comp_grps = nb_comp_grps or trgt_comp_grps

        for (trgt_comp_group, nb_comp_grp) in product(trgt_comp_grps, nb_comp_grps):
            for (trgt_atms_grp, nb_atms_grp) in product(trgt_comp_group.atm_grps, nb_comp_grp.atm_grps):

                if isinstance(self.inter_filter, InteractionFilter):
                    if not self.inter_filter.is_valid_pair(trgt_atms_grp, nb_atms_grp):
                        continue

                if self.add_proximal:
                    all_interactions.extend(self.calc_proximal((trgt_atms_grp, nb_atms_grp)))
                if self.add_cov_inter:
                    pass
                if self.add_clash:
                    pass

                featurePairs = list(product(trgt_atms_grp.features, nb_atms_grp.features))

                featurePairs = filter(lambda x: self.is_feature_pair_valid(*x), featurePairs)

                for featPair in featurePairs:
                    calc_inter_params = (trgt_atms_grp, nb_atms_grp) + featPair
                    interactions = self.resolve_interactions(*calc_inter_params)
                    all_interactions.extend(interactions)

        dependent_interactions = self.find_dependent_interactions(all_interactions)
        all_interactions.extend(dependent_interactions)

        # Get only unique interactions.
        all_interactions = set(all_interactions)

        if not self.add_orphan_h2o_pair:
            all_interactions = self.remove_orphan_h2o_pairs(all_interactions)

        logger.info("Number of potential interactions: %d" % len(all_interactions))

        return list(all_interactions)

    def resolve_interactions(self, group1, group2, feat1, feat2):
        func = self.get_function(feat1.name, feat2.name)
        if (func is None):
            raise IllegalArgumentError("It does not exist a corresponding function to the features: '%s' and '%s'."
                                       % (feat1, feat2))

        interactions = func((group1, group2, feat1, feat2))

        if isinstance(interactions, list) is False:
            interactions = []
        else:
            group1.interactions.extend(interactions)
            group2.interactions.extend(interactions)

        return interactions

    def find_dependent_interactions(self, interactions):
        hbond_set = set()
        attractive_set = set()
        h2o_pairs = defaultdict(dict)
        dependent_interactions = set()

        # Save all hydrogen bonds involving waters and attractive interactions.
        for inter in interactions:
            if inter.type == "Hydrogen bond":
                comp1 = inter.atm_grp1.compound
                comp2 = inter.atm_grp2.compound

                if comp1.is_water():
                    h2o_key = inter.atm_grp1
                    secondKey = inter.atm_grp2

                    if secondKey not in h2o_pairs[h2o_key]:
                        h2o_pairs[h2o_key][secondKey] = inter

                if comp2.is_water():
                    h2o_key = inter.atm_grp2
                    secondKey = inter.atm_grp1

                    if secondKey not in h2o_pairs[h2o_key]:
                        h2o_pairs[h2o_key][secondKey] = inter

                hbond_set.add(inter)
            elif inter.type == "Attractive":
                attractive_set.add(inter)

        for h2o_key in h2o_pairs:
            pairs = combinations(h2o_pairs[h2o_key].keys(), 2)

            for (atm_grp1, atm_grp2) in pairs:
                comp1 = atm_grp1.compound
                comp2 = atm_grp2.compound

                if isinstance(self.inter_filter, InteractionFilter):
                    if not self.inter_filter.is_valid_pair(atm_grp1, atm_grp2):
                        continue

                params = {"depend_of": [h2o_pairs[h2o_key][atm_grp1], h2o_pairs[h2o_key][atm_grp2]]}

                inter = InteractionType(atm_grp1, atm_grp2, "Water-bridged hydrogen bond", params)

                dependent_interactions.add(inter)

        # It will try to match Hydrogen bonds and Attractive interactions
        # involving the same chemical groups to attribute salt bridges.
        sb_groups = set()
        for (hbond, attractive) in product(hbond_set, attractive_set):

            condA = (attractive.atm_grp1.has_atom(hbond.atm_grp1.atoms[0]) and
                     attractive.atm_grp2.has_atom(hbond.atm_grp2.atoms[0]))

            condB = (attractive.atm_grp1.has_atom(hbond.atm_grp2.atoms[0]) and
                     attractive.atm_grp2.has_atom(hbond.atm_grp1.atoms[0]))

            # If an acceptor atom belongs to a negative group, and the donor
            # to a positive group (and vice-versa), it means that the
            # interaction occurs between the same meioties.
            # However, just one condition should occur. For example, it is not
            # possible that an acceptor atom belongs to a negative
            # and positive group at the same time.
            if condA ^ condB:
                key1 = (attractive.atm_grp1, attractive.atm_grp2)
                key2 = (attractive.atm_grp2, attractive.atm_grp1)

                if key1 in sb_groups or key2 in sb_groups:
                    continue

                if isinstance(self.inter_filter, InteractionFilter):
                    if not self.inter_filter.is_valid_pair(attractive.atm_grp1, attractive.atm_grp2):
                        continue

                sb_groups.add(key1)
                params = {"depend_of": [hbond, attractive]}

                inter = InteractionType(attractive.atm_grp1, attractive.atm_grp2, "Salt bridge", params)
                dependent_interactions.add(inter)

        return dependent_interactions

    def remove_orphan_h2o_pairs(self, interactions):
        interactions = set(interactions)
        valid_hbs = set(chain.from_iterable(i.required_interactions for i in interactions))

        orphan_hbs = set([i for i in interactions
                          if ((i.atm_grp1.compound.is_water() or i.atm_grp2.compound.is_water()) and
                              i.type == "Hydrogen bond" and i not in valid_hbs)])

        return interactions - orphan_hbs

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
                ("PositivelyIonizable", "Aromatic"): self.calc_cation_pi}

        # TODO: Incluir:
        # Acceptor - XDonor
        # Acceptor - YDonor
        # Anion - pi system
        # Weak donor - acceptor
        # Weak donor - weak acceptor
        # Hydrogen bond with pi system
        # weak donor - pi system

    def calc_cation_pi(self, params):
        interactions = []

        group1, group2, feat1, feat2 = params
        ccDist = iu.euclidean_distance(group1.centroid, group2.centroid)

        if self.is_within_boundary(ccDist, "boundary_cutoff", le):
            if (self.is_within_boundary(ccDist, "max_dist_cation_pi_inter", le)):
                params = {"dist_cation_pi_inter": ccDist}
                inter = InteractionType(group1, group2, "Cation-pi", params)

                interactions.append(inter)
        return interactions

    def calc_pi_pi(self, params):
        interactions = []

        ring1, ring2, feat1, feat2 = params
        ccDist = iu.euclidean_distance(ring1.centroid, ring2.centroid)

        if self.is_within_boundary(ccDist, "boundary_cutoff", le):
            if (self.is_within_boundary(ccDist, "max_cc_dist_pi_pi_inter", le)):

                dihedralAngle = iu.to_quad1(iu.angle(ring1.normal, ring2.normal))
                vectorCC = ring2.centroid - ring1.centroid
                dispAngle1 = iu.to_quad1(iu.angle(ring1.normal, vectorCC))

                # If the angle criteria were not defined, a specific
                # Pi-stacking definition is not possible as it depends
                # on angle criteria. Therefore, a more general classification
                # is used instead, i.e., all interactions will be Pi-stacking.
                if ("min_dihed_ang_pi_pi_inter" not in self.inter_conf.conf and
                        "max_disp_ang_pi_pi_inter" not in self.inter_conf.conf):
                    interType = "Pi-stacking"
                else:
                    if (self.is_within_boundary(dihedralAngle, "min_dihed_ang_pi_pi_inter", ge)):
                        interType = "Edge-to-face pi-stacking"
                    elif (self.is_within_boundary(dispAngle1, "max_disp_ang_pi_pi_inter", le)):
                        interType = "Face-to-face pi-stacking"
                    else:
                        interType = "Parallel-displaced pi-stacking"

                params = {"cc_dist_pi_pi_inter": ccDist,
                          "dihed_ang_pi_pi_inter": dihedralAngle,
                          "disp_ang_pi_pi_inter": dispAngle1}

                inter = InteractionType(ring1, ring2, interType, params)

                interactions.append(inter)
        return interactions

    def calc_hydrop(self, params):
        interactions = []

        group1, group2, feat1, feat2 = params
        ccDist = iu.euclidean_distance(group1.centroid, group2.centroid)

        if self.is_within_boundary(ccDist, "boundary_cutoff", le):
            if (self.is_within_boundary(ccDist, "max_dist_hydrop_inter", le)):

                params = {"dist_hydrop_inter": ccDist}
                inter = InteractionType(group1, group2, "Hydrophobic", params)

                interactions.append(inter)
        return interactions

    def calc_xbond_pi(self, params):
        interactions = []
        group1, group2, feat1, feat2 = params

        if (feat1.name == "Aromatic" and feat2.name == "HalogenDonor"):
            donorGroup = group2
            ringGroup = group1
        else:
            donorGroup = group1
            ringGroup = group2

        # There are always just one donor/acceptor atom.
        donorAtom = donorGroup.atoms[0]
        # Recover only carbon coordinates
        xCarbonCoords = [x for x in donorAtom.nb_coords.coords if x.atomic_num == 6]

        # Interaction model: C-X ---- A-R
        # Defining the XA distance, in which A is the ring center
        xaDist = iu.euclidean_distance(donorGroup.centroid, ringGroup.centroid)

        if self.is_within_boundary(xaDist, "boundary_cutoff", le):
            if (self.is_within_boundary(xaDist, "max_xc_dist_xbond_inter", le)):

                axVect = donorGroup.centroid - ringGroup.centroid
                dispAngle = iu.to_quad1(iu.angle(ringGroup.normal, axVect))

                if (self.is_within_boundary(dispAngle, "max_disp_ang_xbond_inter", le)):
                    # Interaction model: C-X ---- A-R
                    # XA vector is always the same
                    xaVect = ringGroup.centroid - donorGroup.centroid
                    # Defining angle CXA, in which A is the ring center
                    for carbonCoord in xCarbonCoords:
                        xcVect = carbonCoord.coord - donorGroup.centroid
                        cxaAngle = iu.angle(xcVect, xaVect)

                        if (self.is_within_boundary(cxaAngle, "min_cxa_ang_xbond_inter", ge)):

                            params = {"xc_dist_xbond_inter": xaDist,
                                      "disp_ang_xbond_inter": dispAngle,
                                      "cxa_ang_xbond_inter": cxaAngle}

                            inter = InteractionType(group1, group2, "Halogen bond", params)

                            interactions.append(inter)
        return interactions

    def calc_xbond(self, params):

        interactions = []
        group1, group2, feat1, feat2 = params

        if (feat1.name == "HalogenAcceptor" and feat2.name == "HalogenDonor"):
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
        xCarbonCoords = [x for x in donorAtom.nb_coords.coords if x.atomic_num == 6]

        # Interaction model: C-X ---- A-R
        # Distance XY.
        xaDist = iu.euclidean_distance(donorGroup.centroid, acceptorGroup.centroid)

        if self.is_within_boundary(xaDist, "boundary_cutoff", le):
            if (self.is_within_boundary(xaDist, "max_xa_dist_xbond_inter", le)):

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
                    cxaAngle = iu.angle(xcVect, xaVect)
                    if (self.is_within_boundary(cxaAngle, "min_cxa_ang_xbond_inter", ge)):
                        validCXAAngles.append(cxaAngle)

                # Interaction model: C-X ---- A-R
                # Defining angle RAX
                validXARAngles = []
                # AX vector is always the same
                axVect = donorGroup.centroid - acceptorGroup.centroid
                for nbAtomCoord in acceptorAtom.nb_coords.coords:
                    arVect = nbAtomCoord.coord - acceptorGroup.centroid
                    xarAngle = iu.angle(axVect, arVect)

                    if (self.is_within_boundary(xarAngle, "min_xar_ang_xbond_inter", ge)):
                        validXARAngles.append(xarAngle)

                for anglePair in product(validCXAAngles, validXARAngles):
                    params = {"xa_dist_xbond_inter": xaDist,
                              "cxa_ang_xbond_inter": anglePair[0],
                              "xar_ang_xbond_inter": anglePair[1]}

                    inter = InteractionType(group1, group2, "Halogen bond", params)

                    interactions.append(inter)
        return interactions

    def calc_hbond(self, params):
        interactions = []

        group1, group2, feat1, feat2 = params

        if (feat1.name == "Acceptor" and feat2.name == "Donor"):
            donorGroup = group2
            acceptorGroup = group1
        else:
            donorGroup = group1
            acceptorGroup = group2

        donorAtom = donorGroup.atoms[0]

        daDist = iu.euclidean_distance(donorGroup.centroid, acceptorGroup.centroid)

        if self.is_within_boundary(daDist, "boundary_cutoff", le):
            if self.is_within_boundary(daDist, "max_da_dist_hb_inter", le):

                if donorAtom.get_parent().get_id()[0] == "W":
                    haDist = daDist - 1
                    if (self.is_within_boundary(haDist, "max_ha_dist_hb_inter", le)):
                        params = {"da_dist_hb_inter": daDist,
                                  "ha_dist_hb_inter": -1,
                                  "dha_ang_hb_inter": -1}

                        inter = InteractionType(group1, group2, "Hydrogen bond", params)

                        interactions.append(inter)
                else:
                    # Recover only hydrogen coordinates
                    hydrogCoords = [x for x in donorAtom.nb_coords.coords if x.atomic_num == 1]

                    for hydrogCoord in hydrogCoords:
                        haDist = iu.euclidean_distance(hydrogCoord.coord, acceptorGroup.centroid)

                        dhVect = donorGroup.centroid - hydrogCoord.coord
                        haVect = acceptorGroup.centroid - hydrogCoord.coord
                        dhaAngle = iu.angle(haVect, dhVect)

                        if (self.is_within_boundary(haDist, "max_ha_dist_hb_inter", le) and
                                self.is_within_boundary(dhaAngle, "min_dha_ang_hb_inter", ge)):

                            params = {"da_dist_hb_inter": daDist,
                                      "ha_dist_hb_inter": haDist,
                                      "dha_ang_hb_inter": dhaAngle}

                            inter = InteractionType(group1, group2, "Hydrogen bond", params)

                            interactions.append(inter)
        return interactions

    def calc_attractive(self, params):
        interactions = []

        group1, group2, feat1, feat2 = params
        ccDist = iu.euclidean_distance(group1.centroid, group2.centroid)

        if self.is_within_boundary(ccDist, "boundary_cutoff", le):
            if (self.is_within_boundary(ccDist, "max_dist_attract_inter", le)):

                params = {"dist_attract_inter": ccDist}

                inter = InteractionType(group1, group2, "Attractive", params)

                interactions.append(inter)
        return interactions

    def calc_repulsive(self, params):
        interactions = []

        group1, group2, feat1, feat2 = params
        ccDist = iu.euclidean_distance(group1.centroid, group2.centroid)

        if self.is_within_boundary(ccDist, "boundary_cutoff", le):
            if (self.is_within_boundary(ccDist, "max_dist_repuls_inter", le)):

                params = {"dist_repuls_inter": ccDist}
                inter = InteractionType(group1, group2, "Repulsive", params)

                interactions.append(inter)
        return interactions

    def calc_proximal(self, params):
        interactions = []

        group1, group2 = params

        ccDist = iu.euclidean_distance(group1.centroid, group2.centroid)
        if self.is_within_boundary(ccDist, "boundary_cutoff", le):
            params = {"dist_proximal": ccDist}
            inter = InteractionType(group1, group2, "Proximal", params)
            interactions.append(inter)
            group1.interactions.append(inter)
            group2.interactions.append(inter)

        return interactions

    def is_within_boundary(self, value, key, func):
        if key not in self.inter_conf.conf:
            return True
        else:
            return func(value, self.inter_conf.get_value(key))

    def is_feature_pair_valid(self, feat1, feat2):
        if isinstance(feat1, ChemicalFeature):
            feat1 = feat1.name
        if isinstance(feat2, ChemicalFeature):
            feat2 = feat2.name

        funcs = self.funcs
        return (True if ((feat1, feat2) in funcs or
                         (feat2, feat1) in funcs) else False)

    def get_function(self, feat1, feat2):
        if isinstance(feat1, ChemicalFeature):
            feat1 = feat1.name
        if isinstance(feat2, ChemicalFeature):
            feat2 = feat2.name

        funcs = self.funcs
        if (feat1, feat2) in funcs:
            return funcs[(feat1, feat2)]
        elif (feat2, feat1) in funcs:
            return funcs[(feat2, feat1)]
        else:
            return None

    def set_function(self, pair, func):
        self.funcs[pair] = func


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
        if "depend_of" in i.params:
            remove = False
            for di in i.depend_of:
                if di not in aux_set:
                    remove = True
                    break

            if remove:
                remove_inter.update(set(i.depend_of))
                aux_set.remove(i)
            else:
                keep_inter.update(set(i.depend_of))

    for i in remove_inter:
        if i in aux_set and i not in keep_inter:
            aux_set.remove(i)

    filtered_inter = list(aux_set)

    return filtered_inter
