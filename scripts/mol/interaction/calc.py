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

                params = {"depends_on": [h2o_pairs[h2o_key][atm_grp1], h2o_pairs[h2o_key][atm_grp2]]}

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
                params = {"depends_on": [hbond, attractive]}

                inter = InteractionType(attractive.atm_grp1, attractive.atm_grp2, "Salt bridge", params)
                dependent_interactions.add(inter)

        return dependent_interactions

    def remove_orphan_h2o_pairs(self, interactions):
        interactions = set(interactions)
        valid_hbs = set(chain.from_iterable(i.required_interactions for i in interactions))

        orphan_hbs = set([i for i in interactions
                          if ((i.atm_grp1.compound.is_water() or i.atm_grp2.compound.is_water()) and
                              i.type == "Hydrogen bond" and i not in valid_hbs)])

        interactions -= orphan_hbs

        # Clear the references of each interaction from the AtomGroups objects.
        for i in orphan_hbs:
            i.clear_refs()

        return interactions

    def _default_functions(self):
        return {
                # ("Hydrophobic", "Hydrophobic"): self.calc_hydrop,
                # ("Hydrophobe", "Hydrophobe"): self.calc_hydrop,
                # ("Aromatic", "Aromatic"): self.calc_pi_pi,

                # ("Donor", "Acceptor"): self.calc_hbond,
                ("WeakDonor", "Acceptor"): self.calc_weak_hbond,

                # ("HalogenDonor", "HalogenAcceptor"): self.calc_xbond,
                # ("HalogenDonor", "Aromatic"): self.calc_xbond_pi,

                # ("Negative", "Positive"): self.calc_attractive,
                # ("Negative", "Negative"): self.calc_repulsive,
                # ("Positive", "Positive"): self.calc_repulsive,
                # ("NegIonizable", "PosIonizable"): self.calc_attractive,
                # ("NegIonizable", "NegIonizable"): self.calc_repulsive,
                # ("PosIonizable", "PosIonizable"): self.calc_repulsive,
                # ("PosIonizable", "Aromatic"): self.calc_cation_pi,
                # ("NegativelyIonizable", "PositivelyIonizable"): self.calc_attractive,
                # ("NegativelyIonizable", "NegativelyIonizable"): self.calc_repulsive,
                # ("PositivelyIonizable", "PositivelyIonizable"): self.calc_repulsive,
                # ("PositivelyIonizable", "Aromatic"): self.calc_cation_pi
                }

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
        cc_dist = iu.euclidean_distance(group1.centroid, group2.centroid)

        if self.is_within_boundary(cc_dist, "boundary_cutoff", le):
            if (self.is_within_boundary(cc_dist, "max_dist_cation_pi_inter", le)):
                params = {"dist_cation_pi_inter": cc_dist}
                inter = InteractionType(group1, group2, "Cation-pi", params)

                interactions.append(inter)
        return interactions

    def calc_pi_pi(self, params):
        interactions = []

        ring1, ring2, feat1, feat2 = params
        cc_dist = iu.euclidean_distance(ring1.centroid, ring2.centroid)

        if self.is_within_boundary(cc_dist, "boundary_cutoff", le):
            if (self.is_within_boundary(cc_dist, "max_cc_dist_pi_pi_inter", le)):

                dihedral_angle = iu.to_quad1(iu.angle(ring1.normal, ring2.normal))
                vector_cc = ring2.centroid - ring1.centroid
                disp_angle1 = iu.to_quad1(iu.angle(ring1.normal, vector_cc))

                # If the angle criteria were not defined, a specific
                # Pi-stacking definition is not possible as it depends
                # on angle criteria. Therefore, a more general classification
                # is used instead, i.e., all interactions will be Pi-stacking.
                if ("min_dihed_ang_pi_pi_inter" not in self.inter_conf.conf and
                        "max_disp_ang_pi_pi_inter" not in self.inter_conf.conf):
                    inter_type = "Pi-stacking"
                else:
                    if (self.is_within_boundary(dihedral_angle, "min_dihed_ang_pi_pi_inter", ge)):
                        inter_type = "Edge-to-face pi-stacking"
                    elif (self.is_within_boundary(disp_angle1, "max_disp_ang_pi_pi_inter", le)):
                        inter_type = "Face-to-face pi-stacking"
                    else:
                        inter_type = "Parallel-displaced pi-stacking"

                params = {"cc_dist_pi_pi_inter": cc_dist,
                          "dihed_ang_pi_pi_inter": dihedral_angle,
                          "disp_ang_pi_pi_inter": disp_angle1}

                inter = InteractionType(ring1, ring2, inter_type, params)

                interactions.append(inter)
        return interactions

    def calc_hydrop(self, params):
        interactions = []

        group1, group2, feat1, feat2 = params
        cc_dist = iu.euclidean_distance(group1.centroid, group2.centroid)

        if self.is_within_boundary(cc_dist, "boundary_cutoff", le):
            if (self.is_within_boundary(cc_dist, "max_dist_hydrop_inter", le)):

                params = {"dist_hydrop_inter": cc_dist}
                inter = InteractionType(group1, group2, "Hydrophobic", params)

                interactions.append(inter)
        return interactions

    def calc_xbond_pi(self, params):
        interactions = []
        group1, group2, feat1, feat2 = params

        if (feat1.name == "Aromatic" and feat2.name == "HalogenDonor"):
            donor_grp = group2
            ring_grp = group1
        else:
            donor_grp = group1
            ring_grp = group2

        # There are always just one donor/acceptor atom.
        donor_atm = donor_grp.atoms[0]
        # Recover only carbon coordinates
        x_carbon_coords = [x for x in donor_atm.nb_coords if x.atomic_num == 6]

        # Interaction model: C-X ---- A-R
        # Defining the XA distance, in which A is the ring center
        xa_dist = iu.euclidean_distance(donor_grp.centroid, ring_grp.centroid)

        if self.is_within_boundary(xa_dist, "boundary_cutoff", le):
            if (self.is_within_boundary(xa_dist, "max_xc_dist_xbond_inter", le)):

                ax_vect = donor_grp.centroid - ring_grp.centroid
                disp_angle = iu.to_quad1(iu.angle(ring_grp.normal, ax_vect))

                if (self.is_within_boundary(disp_angle, "max_disp_ang_xbond_inter", le)):
                    # Interaction model: C-X ---- A-R
                    # XA vector is always the same
                    xa_vect = ring_grp.centroid - donor_grp.centroid
                    # Defining angle CXA, in which A is the ring center
                    for coord in x_carbon_coords:
                        xc_vect = coord.vector - donor_grp.centroid
                        cxa_angle = iu.angle(xc_vect, xa_vect)

                        if (self.is_within_boundary(cxa_angle, "min_cxa_ang_xbond_inter", ge)):
                            params = {"xc_dist_xbond_inter": xa_dist,
                                      "disp_ang_xbond_inter": disp_angle,
                                      "cxa_ang_xbond_inter": cxa_angle}

                            inter = InteractionType(group1, group2, "Halogen bond", params)

                            interactions.append(inter)
        return interactions

    def calc_xbond(self, params):

        interactions = []
        group1, group2, feat1, feat2 = params

        if (feat1.name == "HalogenAcceptor" and feat2.name == "HalogenDonor"):
            donor_grp = group2
            acceptor_grp = group1
        else:
            donor_grp = group1
            acceptor_grp = group2

        # There are always just one donor/acceptor atom.
        donor_atm = donor_grp.atoms[0]
        acceptor_atm = acceptor_grp.atoms[0]

        # TODO: Add inside the if
        # Recover only carbon coordinates
        # TODO: adicionar S, P aqui
        x_carbon_coords = [x for x in donor_atm.nb_coords if x.atomic_num == 6]

        # Interaction model: C-X ---- A-R
        # Distance XA.
        xa_dist = iu.euclidean_distance(donor_grp.centroid, acceptor_grp.centroid)

        if self.is_within_boundary(xa_dist, "boundary_cutoff", le):
            if (self.is_within_boundary(xa_dist, "max_xa_dist_xbond_inter", le)):

                # Interaction model: C-X ---- A-R
                valid_CXA_angles = []
                # XA vector is always the same
                xa_vect = acceptor_grp.centroid - donor_grp.centroid

                # Defining the angle CXA
                # It may happen that X is covalently bound to more than one group.
                # In such cases the halogen may also form more than one halogen bond.
                # Ref: Cavallo, G. et al. The Halogen Bond. (2016).
                for coord in x_carbon_coords:
                    xc_vect = coord.vector - donor_grp.centroid
                    cxa_angle = iu.angle(xc_vect, xa_vect)
                    if (self.is_within_boundary(cxa_angle, "min_cxa_ang_xbond_inter", ge)):
                        valid_CXA_angles.append(cxa_angle)

                # Interaction model: C-X ---- A-R
                # Defining angle RAX
                valid_XAR_angles = []
                # AX vector is always the same.
                # Obs: the donor_grp is the Halogen.
                ax_vect = donor_grp.centroid - acceptor_grp.centroid
                for coord in acceptor_atm.nb_coords:
                    ar_vect = coord.vector - acceptor_grp.centroid
                    xar_angle = iu.angle(ax_vect, ar_vect)

                    if (self.is_within_boundary(xar_angle, "min_xar_ang_xbond_inter", ge)):
                        valid_XAR_angles.append(xar_angle)

                for pair in product(valid_CXA_angles, valid_XAR_angles):
                    params = {"xa_dist_xbond_inter": xa_dist,
                              "cxa_ang_xbond_inter": pair[0],
                              "xar_ang_xbond_inter": pair[1]}

                    inter = InteractionType(group1, group2, "Halogen bond", params)
                    interactions.append(inter)
        return interactions

    def calc_hbond(self, params):
        interactions = []

        group1, group2, feat1, feat2 = params

        if (feat1.name == "Acceptor" and feat2.name == "Donor"):
            donor_grp = group2
            acceptor_grp = group1
        else:
            donor_grp = group1
            acceptor_grp = group2

        donor_atm = donor_grp.atoms[0]

        da_dist = iu.euclidean_distance(donor_grp.centroid, acceptor_grp.centroid)

        if self.is_within_boundary(da_dist, "boundary_cutoff", le):
            if self.is_within_boundary(da_dist, "max_da_dist_hb_inter", le):

                if donor_atm.get_parent().is_water():
                    ha_dist = da_dist - 1
                    if (self.is_within_boundary(ha_dist, "max_ha_dist_hb_inter", le)):
                        params = {"da_dist_hb_inter": da_dist,
                                  "ha_dist_hb_inter": -1,
                                  "dha_ang_hb_inter": -1}

                        inter = InteractionType(group1, group2, "Hydrogen bond", params)
                        interactions.append(inter)
                else:
                    # Recover only hydrogen coordinates
                    hydrog_coords = [x for x in donor_atm.nb_coords if x.atomic_num == 1]

                    for coord in hydrog_coords:
                        ha_dist = iu.euclidean_distance(coord.vector, acceptor_grp.centroid)

                        dh_vect = donor_grp.centroid - coord.vector
                        ha_vect = acceptor_grp.centroid - coord.vector
                        dha_angle = iu.angle(ha_vect, dh_vect)

                        if (self.is_within_boundary(ha_dist, "max_ha_dist_hb_inter", le) and
                                self.is_within_boundary(dha_angle, "min_dha_ang_hb_inter", ge)):

                            params = {"da_dist_hb_inter": da_dist,
                                      "ha_dist_hb_inter": ha_dist,
                                      "dha_ang_hb_inter": dha_angle}

                            inter = InteractionType(group1, group2, "Hydrogen bond", params)
                            interactions.append(inter)
        return interactions

    def calc_weak_hbond(self, params):
        interactions = []

        group1, group2, feat1, feat2 = params

        if (feat1.name == "Acceptor" and feat2.name == "WeakDonor"):
            donor_grp = group2
            acceptor_grp = group1
        else:
            donor_grp = group1
            acceptor_grp = group2

        donor_atm = donor_grp.atoms[0]
        acceptor_atm = acceptor_grp.atoms[0]

        da_dist = iu.euclidean_distance(donor_grp.centroid, acceptor_grp.centroid)

        if self.is_within_boundary(da_dist, "boundary_cutoff", le):
            if self.is_within_boundary(da_dist, "max_da_dist_whb_inter", le):

                # Interaction model: D-H ---- A-R.
                # Recover only hydrogen coordinates
                hydrog_coords = [x for x in donor_atm.nb_coords if x.atomic_num == 1]

                # It may happen that D is covalently bound to more than one hydrogen atom.
                # In such cases, it's necessary to check the distances and angles for each atom.
                for h_coord in hydrog_coords:
                    ha_dist = iu.euclidean_distance(h_coord.vector, acceptor_grp.centroid)

                    hd_vect = donor_grp.centroid - h_coord.vector
                    ha_vect = acceptor_grp.centroid - h_coord.vector

                    # TODO: Test if it there is any difference in how the values are passed
                    dha_angle = iu.angle(hd_vect, ha_vect)

                    if (self.is_within_boundary(ha_dist, "max_ha_dist_whb_inter", le) and
                            self.is_within_boundary(dha_angle, "min_dha_ang_whb_inter", ge)):

                        if acceptor_grp.compound.is_water():
                            params = {"da_dist_whb_inter": da_dist,
                                      "ha_dist_whb_inter": ha_dist,
                                      "dha_ang_whb_inter": dha_angle,
                                      "har_ang_whb_inter": -1,
                                      "dar_ang_whb_inter": -1}

                            inter = InteractionType(group1, group2, "Weak hydrogen bond", params)
                            interactions.append(inter)
                        else:
                            # Interaction model: D-H ---- A-R.
                            # R coordinates, in which R is a heavy atom.
                            r_coords = [x for x in acceptor_atm.nb_coords if x.atomic_num != 1]

                            # Interaction model: D-H ---- A-R
                            # AH vector is always the same.
                            ah_vect = h_coord.vector - acceptor_grp.centroid
                            # AD vector is always the same.
                            ad_vect = donor_grp.centroid - acceptor_grp.centroid

                            lower_xar_angle = None
                            lower_dar_angle = None
                            for r_coord in r_coords:
                                ar_vect = r_coord.vector - acceptor_grp.centroid
                                har_angle = iu.angle(ah_vect, ar_vect)

                                # Update the XAR and DAR angles with the lowest values.
                                if lower_xar_angle is None or har_angle < lower_xar_angle:
                                    lower_xar_angle = har_angle
                                    # It is not necessary to check if it is lower, because HAR is already checked.
                                    # It is sufficient to check only HAR.
                                    lower_dar_angle = iu.angle(ad_vect, ar_vect)

                            if lower_xar_angle is None or lower_dar_angle is None:
                                lower_xar_angle = -1
                                lower_dar_angle = -1

                            if (self.is_within_boundary(lower_xar_angle, "min_har_ang_whb_inter", ge) and
                                    self.is_within_boundary(lower_dar_angle, "min_dar_ang_whb_inter", ge)):

                                params = {"da_dist_whb_inter": da_dist,
                                          "ha_dist_whb_inter": ha_dist,
                                          "dha_ang_whb_inter": dha_angle,
                                          "har_ang_whb_inter": lower_xar_angle,
                                          "dar_ang_whb_inter": lower_dar_angle}

                                inter = InteractionType(group1, group2, "Weak hydrogen bond", params)
                                interactions.append(inter)

    def calc_attractive(self, params):
        group1, group2, feat1, feat2 = params

        interactions = []
        cc_dist = iu.euclidean_distance(group1.centroid, group2.centroid)
        if self.is_within_boundary(cc_dist, "boundary_cutoff", le):
            if (self.is_within_boundary(cc_dist, "max_dist_attract_inter", le)):

                params = {"dist_attract_inter": cc_dist}
                inter = InteractionType(group1, group2, "Attractive", params)

                interactions.append(inter)
        return interactions

    def calc_repulsive(self, params):
        group1, group2, feat1, feat2 = params

        interactions = []
        cc_dist = iu.euclidean_distance(group1.centroid, group2.centroid)
        if self.is_within_boundary(cc_dist, "boundary_cutoff", le):
            if (self.is_within_boundary(cc_dist, "max_dist_repuls_inter", le)):

                params = {"dist_repuls_inter": cc_dist}
                inter = InteractionType(group1, group2, "Repulsive", params)

                interactions.append(inter)
        return interactions

    def calc_proximal(self, params):
        group1, group2 = params

        interactions = []
        cc_dist = iu.euclidean_distance(group1.centroid, group2.centroid)
        if self.is_within_boundary(cc_dist, "max_dist_proximal", le):
            params = {"dist_proximal": cc_dist}
            inter = InteractionType(group1, group2, "Proximal", params)
            interactions.append(inter)

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
        if "depends_on" in i.params:
            remove = False
            for di in i.depends_on:
                if di not in aux_set:
                    remove = True
                    break

            if remove:
                remove_inter.update(set(i.depends_on))
                aux_set.remove(i)
            else:
                keep_inter.update(set(i.depends_on))

    for i in remove_inter:
        if i in aux_set and i not in keep_inter:
            aux_set.remove(i)

    filtered_inter = list(aux_set)

    return filtered_inter
