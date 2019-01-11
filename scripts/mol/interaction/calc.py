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
                 add_clash=False, add_orphan_h2o_pair=False, strict_donor_rules=False):

        self.inter_conf = inter_conf

        self.add_proximal = add_proximal
        self.add_cov_inter = add_cov_inter
        self.add_clash = add_clash
        self.add_orphan_h2o_pair = add_orphan_h2o_pair
        self.strict_donor_rules = strict_donor_rules

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
                comp1 = inter.atm_grp1.compounds[0]
                comp2 = inter.atm_grp2.compounds[0]

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
                comp1 = atm_grp1.compounds[0]
                comp2 = atm_grp2.compounds[0]

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
                          if ((i.atm_grp1.is_water() or i.atm_grp2.is_water()) and
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
                # ("WeakDonor", "Acceptor"): self.calc_weak_hbond,
                # ("Donor", "Aromatic"): self.calc_hbond_pi,
                # ("WeakDonor", "Aromatic"): self.calc_hbond_pi,

                # ("HalogenDonor", "Acceptor"): self.calc_xbond,
                # ("HalogenDonor", "Aromatic"): self.calc_xbond_pi,

                ("Amide", "Aromatic"): self.calc_amide_pi,


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
        # Chalcogen bond
        # Anion - pi system

        # Weak donor - weak acceptor

        # Dipole-dipole interaction
        # Covalent interaction
        # Clash
        # Amide-pi bond
        # Agostic and Hydrogen-Bonding X–H· · · M

        # Sulfur interactions - Bissantz 2010; Taylor 2016

        # Metalic complex

        # REF: https://onlinelibrary.wiley.com/doi/epdf/10.1002/anie.200390319
        # aromatic between hbond arrays

    def calc_cation_pi(self, params):
        interactions = []

        group1, group2, feat1, feat2 = params
        cc_dist = iu.euclidean_distance(group1.centroid, group2.centroid)

        if (self.is_within_boundary(cc_dist, "boundary_cutoff", le) and
                self.is_within_boundary(cc_dist, "max_dist_cation_pi_inter", le)):
                params = {"dist_cation_pi_inter": cc_dist}
                inter = InteractionType(group1, group2, "Cation-pi", params)

                interactions.append(inter)
        return interactions

    def calc_pi_pi(self, params):
        interactions = []

        ring1, ring2, feat1, feat2 = params
        cc_dist = iu.euclidean_distance(ring1.centroid, ring2.centroid)

        if (self.is_within_boundary(cc_dist, "boundary_cutoff", le) and
                self.is_within_boundary(cc_dist, "max_cc_dist_pi_pi_inter", le)):

                dihedral_angle = iu.to_quad1(iu.angle(ring1.normal, ring2.normal))
                vector_cc = ring2.centroid - ring1.centroid
                disp_angle = iu.to_quad1(iu.angle(ring1.normal, vector_cc))

                # If the angle criteria were not defined, a specific Pi-stacking definition is not
                # possible as it depends on angle criteria. Therefore, a more general classification
                # is used instead, i.e., all interactions will be Pi-stacking.
                if ("min_dihed_ang_pi_pi_inter" not in self.inter_conf.conf and
                        "max_disp_ang_pi_pi_inter" not in self.inter_conf.conf):
                    inter_type = "Pi-stacking"
                else:
                    if (self.is_within_boundary(dihedral_angle, "min_dihed_ang_pi_pi_inter", ge)):
                        inter_type = "Edge-to-face pi-stacking"
                    elif (self.is_within_boundary(disp_angle, "max_disp_ang_pi_pi_inter", le)):
                        inter_type = "Face-to-face pi-stacking"
                    else:
                        inter_type = "Parallel-displaced pi-stacking"

                params = {"cc_dist_pi_pi_inter": cc_dist,
                          "dihed_ang_pi_pi_inter": dihedral_angle,
                          "disp_ang_pi_pi_inter": disp_angle}

                inter = InteractionType(ring1, ring2, inter_type, params)

                interactions.append(inter)
        return interactions

    def calc_amide_pi(self, params):
        interactions = []

        group1, group2, feat1, feat2 = params

        if (feat1.name == "Aromatic" and feat2.name == "Amide"):
            ring_grp = group1
            amide_grp = group2
        elif (feat2.name == "Aromatic" and feat1.name == "Amide"):
            ring_grp = group2
            amide_grp = group1
        else:
            logger.warning("Amide-aromatic interactions requires an aromatic and and amide group.")
            logger.warning("However, the informed groups have the features '%s' and '%s'" %
                           (group1.features, group2.features))
            return []

        # Distance between the amide and ring centroids.
        cc_dist = iu.euclidean_distance(ring_grp.centroid, amide_grp.centroid)

        print("#############################")
        print(ring_grp)
        print(ring_grp.compounds)
        print("Amide:")
        print(amide_grp)
        print(amide_grp.compounds)
        print()
        print("CC dist: ", cc_dist)
        print()

        if (self.is_within_boundary(cc_dist, "boundary_cutoff", le) and
                self.is_within_boundary(cc_dist, "max_cc_dist_amide_pi_inter", le)):

            dihedral_angle = iu.to_quad1(iu.angle(ring_grp.normal, amide_grp.normal))
            vector_cc = amide_grp.centroid - ring_grp.centroid
            disp_angle = iu.to_quad1(iu.angle(ring_grp.normal, vector_cc))

            print("Dihedral: ", dihedral_angle)
            print("Displacement: ", disp_angle)
            print()

            if (self.is_within_boundary(dihedral_angle, "max_dihed_ang_amide_pi_inter", le) and
                    self.is_within_boundary(disp_angle, "max_disp_ang_pi_pi_inter", le)):

                print("Valid interaction")
                print()

                params = {"cc_dist_amide_pi_inter": cc_dist,
                          "dihed_ang_amide_pi_inter": dihedral_angle,
                          "disp_ang_amide_pi_inter": disp_angle}

                inter = InteractionType(group1, group2, "Amide-aromatic stacking", params)

                interactions.append(inter)

        # print()

        return interactions

    def calc_hydrop(self, params):
        interactions = []

        group1, group2, feat1, feat2 = params
        cc_dist = iu.euclidean_distance(group1.centroid, group2.centroid)

        if (self.is_within_boundary(cc_dist, "boundary_cutoff", le) and
                self.is_within_boundary(cc_dist, "max_dist_hydrop_inter", le)):

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

        # Interaction model: C-X ---- A
        # Defining the XA distance, in which A is the ring center
        xa_dist = iu.euclidean_distance(donor_grp.centroid, ring_grp.centroid)

        if (self.is_within_boundary(xa_dist, "boundary_cutoff", le) and
                self.is_within_boundary(xa_dist, "max_xc_dist_xbond_inter", le)):

            ax_vect = donor_grp.centroid - ring_grp.centroid
            disp_angle = iu.to_quad1(iu.angle(ring_grp.normal, ax_vect))

            if (self.is_within_boundary(disp_angle, "max_disp_ang_xbond_inter", le)):
                # Interaction model: C-X ---- A
                # XA vector is always the same
                xa_vect = ring_grp.centroid - donor_grp.centroid

                # Defining angle CXA, in which A is the ring center
                # It may happen that X is covalently bound to more than one group.
                # In such cases the halogen may also form more than one halogen bond.
                # Ref: Cavallo, G. et al. The Halogen Bond. (2016).
                carbon_coords = [x for x in donor_atm.nb_coords if x.atomic_num == 6]
                for coord in carbon_coords:
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

        if len(group1.atoms) != 1 or len(group2.atoms) != 1:
            logger.warning("One or more invalid atom groups were informed: '%s' and '%s" % (group1, group2))
            logger.warning("In halogen bonds, halogen donor and acceptor groups should always contain only one atom.")
            return []

        if (feat1.name == "Acceptor" and feat2.name == "HalogenDonor"):
            donor_grp = group2
            acceptor_grp = group1
        else:
            donor_grp = group1
            acceptor_grp = group2

        # There are always just one donor/acceptor atom.
        donor_atm = donor_grp.atoms[0]
        acceptor_atm = acceptor_grp.atoms[0]

        # Interaction model: C-X ---- A-R
        # Distance XA.
        xa_dist = iu.euclidean_distance(donor_grp.centroid, acceptor_grp.centroid)

        if (self.is_within_boundary(xa_dist, "boundary_cutoff", le) and
                self.is_within_boundary(xa_dist, "max_xa_dist_xbond_inter", le)):

            # Interaction model: C-X ---- A-R
            # XA vector is always the same
            xa_vect = acceptor_grp.centroid - donor_grp.centroid

            # Defining the angle CXA
            # It may happen that X is covalently bound to more than one group.
            # In such cases the halogen may also form more than one halogen bond.
            # Ref: Cavallo, G. et al. The Halogen Bond. (2016).
            carbon_coords = [x for x in donor_atm.nb_coords if x.atomic_num == 6]

            # Interaction model: C-X ---- A-R.
            # R coordinates, in which R is a heavy atom.
            r_coords = [x for x in acceptor_atm.nb_coords if x.atomic_num != 1]

            for coord in carbon_coords:
                xc_vect = coord.vector - donor_grp.centroid
                cxa_angle = iu.angle(xc_vect, xa_vect)

                if (self.is_within_boundary(cxa_angle, "min_cxa_ang_xbond_inter", ge)):
                    # If no heavy atom is bonded to the acceptor, it means that only hydrogens
                    # may be bonded to it. Then, we do not need to calculate the angles because
                    # hydrogens are too dynamic what means that an acceptor could be ionized or not
                    # at a specific moment in time. That is the reason for why the strict rule is
                    # only applied to donor atoms.
                    if len(r_coords) == 0:
                        params = {"xa_dist_xbond_inter": xa_dist,
                                  "cxa_ang_xbond_inter": cxa_angle,
                                  "xar_ang_xbond_inter": -1}

                        inter = InteractionType(group1, group2, "Halogen bond", params)
                        interactions.append(inter)
                    else:
                        # AX vector is always the same.
                        # Obs: the donor_grp is the Halogen (X).
                        ax_vect = donor_grp.centroid - acceptor_grp.centroid

                        lowest_xar_angle = None
                        for r_coord in r_coords:
                            ar_vect = r_coord.vector - acceptor_grp.centroid
                            xar_angle = iu.angle(ax_vect, ar_vect)

                            # Update the XAR angle with the lowest value.
                            if lowest_xar_angle is None or xar_angle < lowest_xar_angle:
                                lowest_xar_angle = xar_angle

                        # The angle will be None when any R (heavy atom) atom was found.
                        # In this case, the criteria must always fail.
                        if lowest_xar_angle is None:
                            lowest_xar_angle = -1

                        if self.is_within_boundary(lowest_xar_angle, "min_xar_ang_xbond_inter", ge):
                            params = {"xa_dist_xbond_inter": xa_dist,
                                      "cxa_ang_xbond_inter": cxa_angle,
                                      "xar_ang_xbond_inter": lowest_xar_angle}

                            inter = InteractionType(group1, group2, "Halogen bond", params)
                            interactions.append(inter)

        return interactions

    def calc_hbond(self, params):
        interactions = []
        group1, group2, feat1, feat2 = params

        if len(group1.atoms) != 1 or len(group2.atoms) != 1:
            logger.warning("One or more invalid atom groups were informed: '%s' and '%s" % (group1, group2))
            logger.warning("In hydrogen bonds, donor and acceptor groups should always contain only one atom.")
            return []

        if (feat1.name == "Acceptor" and feat2.name == "Donor"):
            donor_grp = group2
            acceptor_grp = group1
        else:
            donor_grp = group1
            acceptor_grp = group2

        donor_atm = donor_grp.atoms[0]
        acceptor_atm = acceptor_grp.atoms[0]

        # Interaction model: D-H ---- A-R.
        # DA distance
        da_dist = iu.euclidean_distance(donor_grp.centroid, acceptor_grp.centroid)

        if (self.is_within_boundary(da_dist, "boundary_cutoff", le) and
                self.is_within_boundary(da_dist, "max_da_dist_hb_inter", le)):

            # Interaction model: D-H ---- A-R.
            # Recover only hydrogen coordinates bonded to the donor.
            hydrog_coords = [x for x in donor_atm.nb_coords if x.atomic_num == 1]

            # Interaction model: D-H ---- A-R.
            # R coordinates, in which R is a heavy atom.
            r_coords = [x for x in acceptor_atm.nb_coords if x.atomic_num != 1]

            # Firstly, it checks if it is not necessary to apply a strict hbond rule, i.e.,
            # hydrogens must exist and all geometrical criteria should be evaluated.
            # Then it checks if no hydrogen is bonded to the donor, or if the donor has hydrogens
            # and only hydrogens as neighbours (water, solvents, ammonia, SH2). In the latter case, the
            # hydrogens can be positioned in many different set of ways, and each run of a tool like
            # OpenBabel would vary the hydrogen bond list when one applies this algorithm.
            if (self.strict_donor_rules is False and
                    (len(hydrog_coords) == 0 or
                        len(hydrog_coords) == len(donor_atm.nb_coords))):

                # When the position of the hydrogen cannot be defined, it assumes the hydrogen to be located 1A
                # away from the donor in a line formed by the donor and the acceptor.
                ha_dist = da_dist - 1
                if self.is_within_boundary(ha_dist, "max_ha_dist_hb_inter", le):

                    # If no heavy atom is bonded to the acceptor, it means that only hydrogens
                    # may be bonded to it. Then, we do not need to calculate the angles because
                    # hydrogens are too dynamic, what means that an acceptor could be ionized or not
                    # at a specific moment in time and its hydrogens may be positioned in different ways.
                    # That is the reason for why the strict rule is only applied to donor atoms.
                    if len(r_coords) == 0:
                        params = {"da_dist_hb_inter": da_dist,
                                  "ha_dist_hb_inter": -1,
                                  "dha_ang_hb_inter": -1,
                                  "har_ang_hb_inter": -1,
                                  "dar_ang_hb_inter": -1}

                        inter = InteractionType(group1, group2, "Hydrogen bond", params)
                        interactions.append(inter)
                    else:
                        # AD vector is always the same.
                        ad_vect = donor_grp.centroid - acceptor_grp.centroid

                        lowest_dar_angle = None
                        for r_coord in r_coords:
                            ar_vect = r_coord.vector - acceptor_grp.centroid
                            dar_angle = iu.angle(ad_vect, ar_vect)

                            # Update the DAR angle with the lowest value.
                            if lowest_dar_angle is None or dar_angle < lowest_dar_angle:
                                lowest_dar_angle = dar_angle

                        # The angle will be None when any R (heavy atom) atom was found.
                        # In this case, the criteria must always fail.
                        if lowest_dar_angle is None:
                            lowest_dar_angle = -1

                        if self.is_within_boundary(lowest_dar_angle, "min_dar_ang_hb_inter", ge):
                            params = {"da_dist_hb_inter": da_dist,
                                      "ha_dist_hb_inter": -1,
                                      "dha_ang_hb_inter": -1,
                                      "har_ang_hb_inter": -1,
                                      "dar_ang_hb_inter": lowest_dar_angle}

                            inter = InteractionType(group1, group2, "Hydrogen bond", params)
                            interactions.append(inter)
            else:
                # It may happen that D is covalently bound to more than one hydrogen atom.
                # In this case, it is necessary to check the distances and angles for each atom.
                # It will produce a hydrogen bond for each valid hydrogen.
                for h_coord in hydrog_coords:
                    ha_dist = iu.euclidean_distance(h_coord.vector, acceptor_grp.centroid)

                    hd_vect = donor_grp.centroid - h_coord.vector
                    ha_vect = acceptor_grp.centroid - h_coord.vector
                    dha_angle = iu.angle(hd_vect, ha_vect)

                    if (self.is_within_boundary(ha_dist, "max_ha_dist_hb_inter", le) and
                            self.is_within_boundary(dha_angle, "min_dha_ang_hb_inter", ge)):

                        # If no heavy atom is bonded to the acceptor, it means that only hydrogens
                        # may be bonded to it. Then, we do not need to calculate the angles because
                        # hydrogens are too dynamic what means that an acceptor could be ionized or not
                        # at a specific moment in time. That is the reason for why the strict rule is
                        # only applied to donor atoms.
                        if len(r_coords) == 0:
                            params = {"da_dist_hb_inter": da_dist,
                                      "ha_dist_hb_inter": ha_dist,
                                      "dha_ang_hb_inter": dha_angle,
                                      "har_ang_hb_inter": -1,
                                      "dar_ang_hb_inter": -1}

                            inter = InteractionType(group1, group2, "Hydrogen bond", params)
                            interactions.append(inter)
                        else:
                            # Interaction model: D-H ---- A-R
                            # AH vector is always the same.
                            ah_vect = h_coord.vector - acceptor_grp.centroid
                            # AD vector is always the same.
                            ad_vect = donor_grp.centroid - acceptor_grp.centroid

                            # Interaction model: D-H ---- A-R
                            # Check the angles formed at the acceptor.
                            # When A is covalently bonded to more than one R atom, it is necessary to
                            # evaluate all possible angles D-A-R and H-A-R. In this case, all angles should
                            # satisfy the angle criteria. To do so, we could analyze only the lowest D-A-R and
                            # H-A-R angles. It guarantees that all angles will satisfy the criteria.
                            # OBS: it may happen that each one of the angles would belong to a different R atom.
                            lowest_har_angle = None
                            lowest_dar_angle = None
                            for r_coord in r_coords:
                                ar_vect = r_coord.vector - acceptor_grp.centroid
                                har_angle = iu.angle(ah_vect, ar_vect)
                                dar_angle = iu.angle(ad_vect, ar_vect)

                                # Update the HAR angle with the lowest value.
                                if lowest_har_angle is None or har_angle < lowest_har_angle:
                                    lowest_har_angle = har_angle
                                # Update the DAR angle with the lowest value.
                                if lowest_dar_angle is None or dar_angle < lowest_dar_angle:
                                    lowest_dar_angle = dar_angle

                            # The angles will be None when any R (heavy atom) atom was found.
                            # In this case, the criteria must always fail.
                            if lowest_har_angle is None or lowest_dar_angle is None:
                                lowest_har_angle = -1
                                lowest_dar_angle = -1

                            if (self.is_within_boundary(lowest_har_angle, "min_har_ang_hb_inter", ge) and
                                    self.is_within_boundary(lowest_dar_angle, "min_dar_ang_hb_inter", ge)):

                                # Only the lowest D-A-R and H-A-R angles are provided.
                                params = {"da_dist_hb_inter": da_dist,
                                          "ha_dist_hb_inter": ha_dist,
                                          "dha_ang_hb_inter": dha_angle,
                                          "har_ang_hb_inter": lowest_har_angle,
                                          "dar_ang_hb_inter": lowest_dar_angle}

                                inter = InteractionType(group1, group2, "Hydrogen bond", params)
                                interactions.append(inter)

        return interactions

    def calc_weak_hbond(self, params):
        interactions = []

        group1, group2, feat1, feat2 = params

        if len(group1.atoms) != 1 or len(group2.atoms) != 1:
            logger.warning("One or more invalid atom groups were informed: '%s' and '%s" % (group1, group2))
            logger.warning("In weak hydrogen bonds, weak donor and acceptor groups should always contain only one atom.")
            return []

        if (feat1.name == "Acceptor" and feat2.name == "WeakDonor"):
            donor_grp = group2
            acceptor_grp = group1
        else:
            donor_grp = group1
            acceptor_grp = group2

        donor_atm = donor_grp.atoms[0]
        acceptor_atm = acceptor_grp.atoms[0]

        da_dist = iu.euclidean_distance(donor_grp.centroid, acceptor_grp.centroid)
        if (self.is_within_boundary(da_dist, "boundary_cutoff", le) and
                self.is_within_boundary(da_dist, "max_da_dist_whb_inter", le)):

            # Interaction model: D-H ---- A-R.
            # Recover only hydrogen coordinates bonded to the donor.
            hydrog_coords = [x for x in donor_atm.nb_coords if x.atomic_num == 1]

            # Interaction model: D-H ---- A-R.
            # R coordinates, in which R is a heavy atom.
            r_coords = [x for x in acceptor_atm.nb_coords if x.atomic_num != 1]

            # It may happen that D is covalently bound to more than one hydrogen atom.
            # In such cases, it's necessary to check the distances and angles for each atom.
            for h_coord in hydrog_coords:
                ha_dist = iu.euclidean_distance(h_coord.vector, acceptor_grp.centroid)

                hd_vect = donor_grp.centroid - h_coord.vector
                ha_vect = acceptor_grp.centroid - h_coord.vector
                dha_angle = iu.angle(hd_vect, ha_vect)

                if (self.is_within_boundary(ha_dist, "max_ha_dist_whb_inter", le) and
                        self.is_within_boundary(dha_angle, "min_dha_ang_whb_inter", ge)):

                    # If no heavy atom is bonded to the acceptor, it means that only hydrogens
                    # may be bonded to it. Then, we do not need to calculate the angles because
                    # hydrogens are too dynamic what means that an acceptor could be ionized or not
                    # at a specific moment in time. That is the reason for why the strict rule is
                    # only applied to donor atoms.
                    if len(r_coords) == 0:
                        params = {"da_dist_whb_inter": da_dist,
                                  "ha_dist_whb_inter": ha_dist,
                                  "dha_ang_whb_inter": dha_angle,
                                  "har_ang_whb_inter": -1,
                                  "dar_ang_whb_inter": -1}

                        inter = InteractionType(group1, group2, "Weak hydrogen bond", params)
                        interactions.append(inter)
                    else:
                        # Interaction model: D-H ---- A-R
                        # AH vector is always the same.
                        ah_vect = h_coord.vector - acceptor_grp.centroid
                        # AD vector is always the same.
                        ad_vect = donor_grp.centroid - acceptor_grp.centroid

                        # Interaction model: D-H ---- A-R
                        # Check the angles formed at the acceptor.
                        # When A is covalently bonded to more than one R atom, it is necessary to
                        # evaluate all possible angles D-A-R and H-A-R. In this case, all angles should
                        # satisfy the angle criteria. To do so, we could analyze only the lowest D-A-R and
                        # H-A-R angles. It guarantees that all angles will satisfy the criteria.
                        # OBS: it may happen that each one of the angles would belong to a different R atom.
                        lowest_har_angle = None
                        lowest_dar_angle = None
                        for r_coord in r_coords:
                            ar_vect = r_coord.vector - acceptor_grp.centroid
                            har_angle = iu.angle(ah_vect, ar_vect)
                            dar_angle = iu.angle(ad_vect, ar_vect)

                            # Update the HAR angle with the lowest value.
                            if lowest_har_angle is None or har_angle < lowest_har_angle:
                                lowest_har_angle = har_angle
                            # Update the DAR angle with the lowest value.
                            if lowest_dar_angle is None or dar_angle < lowest_dar_angle:
                                lowest_dar_angle = dar_angle

                        # The angles will be None when any R (heavy atom) atom was found.
                        # In this case, the criteria must always fail.
                        if lowest_har_angle is None or lowest_dar_angle is None:
                            lowest_har_angle = -1
                            lowest_dar_angle = -1

                        if (self.is_within_boundary(lowest_har_angle, "min_har_ang_whb_inter", ge) and
                                self.is_within_boundary(lowest_dar_angle, "min_dar_ang_whb_inter", ge)):

                            # Only the lowest D-A-R and H-A-R angles are provided.
                            params = {"da_dist_whb_inter": da_dist,
                                      "ha_dist_whb_inter": ha_dist,
                                      "dha_ang_whb_inter": dha_angle,
                                      "har_ang_whb_inter": lowest_har_angle,
                                      "dar_ang_whb_inter": lowest_dar_angle}

                            inter = InteractionType(group1, group2, "Weak hydrogen bond", params)
                            interactions.append(inter)
        return interactions

    def calc_hbond_pi(self, params):
        interactions = []
        group1, group2, feat1, feat2 = params

        if (feat1.name == "Aromatic" and (feat2.name == "Donor" or feat2.name == "WeakDonor")):
            ring_grp = group1
            donor_grp = group2
        elif (feat2.name == "Aromatic" and (feat1.name == "Donor" or feat1.name == "WeakDonor")):
            ring_grp = group2
            donor_grp = group1
        else:
            logger.warning("Hydrogen bond involving pi-systems requires an aromatic and donor (weak donor) groups.")
            logger.warning("However, the informed groups have the features '%s' and '%s'" %
                           (group1.features, group2.features))
            return []

        if len(donor_grp.atoms) != 1:
            logger.warning("Invalid donor (weak donor) group was informed: '%s'" % donor_grp)
            logger.warning("In hydrogen bonds involving pi-systems, donor (weak donor) groups should always contain only one atom.")
            return []

        # There are always just one donor/acceptor atom.
        donor_atm = donor_grp.atoms[0]

        # Interaction model: D-H ---- A, in which A is the ring center.
        da_dist = iu.euclidean_distance(donor_grp.centroid, ring_grp.centroid)
        if (self.is_within_boundary(da_dist, "boundary_cutoff", le) and
                self.is_within_boundary(da_dist, "max_dc_dist_whb_inter", le)):

            # Interaction model: D-H ---- A, in which A is the ring center.
            # Recover only hydrogen coordinates bonded to the donor.
            hydrog_coords = [x for x in donor_atm.nb_coords if x.atomic_num == 1]

            # Firstly, it checks if it is not necessary to apply a strict hbond rule, i.e.,
            # hydrogens must exist and all geometrical criteria should be evaluated.
            # Then it checks if no hydrogen is bonded to the donor, or if the donor has hydrogens
            # and only hydrogens as neighbours (water, solvents, ammonia, SH2). In the latter case, the
            # hydrogens can be positioned in many different set of ways, and each run of a tool like
            # OpenBabel would vary the hydrogen bond list when one applies this algorithm.
            if (self.strict_donor_rules is False and
                    (len(hydrog_coords) == 0 or
                        len(hydrog_coords) == len(donor_atm.nb_coords))):

                # When the position of the hydrogen cannot be defined, it assumes the hydrogen to be located 1A
                # away from the donor in a line formed by the donor and the acceptor.
                ha_dist = da_dist - 1
                if self.is_within_boundary(ha_dist, "max_hc_dist_whb_inter", le):

                    # Interaction model: D-H ---- A, in which A is the ring center.
                    # Calculate the displacement angle formed between the ring normal and the vector Donor-Centroid.
                    ad_vect = donor_grp.centroid - ring_grp.centroid
                    disp_angle = iu.to_quad1(iu.angle(ring_grp.normal, ad_vect))

                    if (self.is_within_boundary(disp_angle, "max_disp_ang_whb_inter", le)):
                        params = {"dc_dist_whb_inter": da_dist,
                                  "hc_dist_whb_inter": -1,
                                  "dhc_ang_whb_inter": -1,
                                  "disp_ang_whb_inter": disp_angle}

                        inter = InteractionType(group1, group2, "Weak hydrogen bond", params)
                        interactions.append(inter)
            else:
                # It may happen that D is covalently bound to more than one hydrogen atom.
                # In this case, it is necessary to check the distances and angles for each atom.
                # It will produce a hydrogen bond for each valid hydrogen.
                for h_coord in hydrog_coords:
                    ha_dist = iu.euclidean_distance(h_coord.vector, ring_grp.centroid)

                    hd_vect = donor_grp.centroid - h_coord.vector
                    ha_vect = ring_grp.centroid - h_coord.vector
                    dha_angle = iu.angle(hd_vect, ha_vect)

                    if (self.is_within_boundary(ha_dist, "max_hc_dist_whb_inter", le) and
                            self.is_within_boundary(dha_angle, "min_dhc_ang_whb_inter", ge)):

                        # Interaction model: D-H ---- A, in which A is the ring center.
                        # Calculate the displacement angle formed between the ring normal and the vector Donor-Centroid.
                        ad_vect = donor_grp.centroid - ring_grp.centroid
                        disp_angle = iu.to_quad1(iu.angle(ring_grp.normal, ad_vect))

                        if (self.is_within_boundary(disp_angle, "max_disp_ang_whb_inter", le)):
                            params = {"dc_dist_whb_inter": da_dist,
                                      "hc_dist_whb_inter": ha_dist,
                                      "dhc_ang_whb_inter": dha_angle,
                                      "disp_ang_whb_inter": disp_angle}

                            inter = InteractionType(group1, group2, "Weak hydrogen bond", params)
                            interactions.append(inter)

        return interactions

    def calc_attractive(self, params):
        group1, group2, feat1, feat2 = params

        interactions = []
        cc_dist = iu.euclidean_distance(group1.centroid, group2.centroid)
        if (self.is_within_boundary(cc_dist, "boundary_cutoff", le) and
                self.is_within_boundary(cc_dist, "max_dist_attract_inter", le)):

            params = {"dist_attract_inter": cc_dist}
            inter = InteractionType(group1, group2, "Attractive", params)

            interactions.append(inter)
        return interactions

    def calc_repulsive(self, params):
        group1, group2, feat1, feat2 = params

        interactions = []
        cc_dist = iu.euclidean_distance(group1.centroid, group2.centroid)
        if (self.is_within_boundary(cc_dist, "boundary_cutoff", le) and
                self.is_within_boundary(cc_dist, "max_dist_repuls_inter", le)):

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
