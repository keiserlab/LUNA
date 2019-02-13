from util.exceptions import ShellCenterNotFound

from util.default_values import (CHEMICAL_FEATURES_IDS, INTERACTIONS_IDS)

from Bio.KDTree import KDTree

from itertools import (chain, product)
from collections import defaultdict

from mol.interaction.fp.fingerprint import (DEFAULT_SHELL_NBITS, Fingerprint, CountFingerprint)

import numpy as np
import mmh3
import logging

logger = logging.getLogger()


class ShellSearch:

    def __init__(self, atm_grps, bucket_size=10):
        self.atm_grps = atm_grps

        # get the coordinates
        coord_list = [ga.centroid for ga in self.atm_grps]

        # to Nx3 array of type float
        self.coords = np.array(coord_list).astype("f")
        assert(bucket_size > 1)
        assert(self.coords.shape[1] == 3)
        self.kdt = KDTree(3, bucket_size)
        self.kdt.set_coords(self.coords)

    def search(self, center, radius):

        self.kdt.search(center, radius)
        indices = self.kdt.get_indices()
        n_grps_list = []
        atm_grps = self.atm_grps
        for i in indices:
            a = atm_grps[i]
            n_grps_list.append(a)

        return n_grps_list


class ShellManager:

    def __init__(self, num_levels, radius_step, num_bits, shells=None, full_control=True, verbose=False):
        self.num_levels = num_levels
        self.radius_steps = radius_step
        self.num_bits = num_bits

        self.shells = shells or []
        self.verbose = verbose

        self._init_controllers()

    @property
    def num_shells(self):
        return len(self.shells)

    @property
    def num_unique_shells(self):
        return len(self.unique_shells)

    @property
    def unique_shells(self):
        return self.get_valid_shells()

    def find_similar_shell(self, shell):
        found_shell = None
        for added_shell in self.shells:
            # For each group of similar shells, it will exist only one valid shell.
            if added_shell.is_valid:
                # If only one of the shells has interactions, then they are different.
                # But, if none of them have interactions, check if the spheres have equal identifiers.
                if not added_shell.interactions and not shell.interactions:
                    if added_shell.identifier == shell.identifier:
                        found_shell = added_shell
                        break
                # If both shells have interactions, check if the interactions are equal.
                elif added_shell.interactions and shell.interactions:
                    if added_shell.interactions == shell.interactions:
                        found_shell = added_shell
                        break

        return found_shell

    def add_shell(self, shell):
        found_shell = self.find_similar_shell(shell)

        if found_shell:
            if found_shell.level < shell.level:
                shell.is_valid = False
            elif found_shell.level == shell.level and found_shell.identifier <= shell.identifier:
                shell.is_valid = False
            else:
                found_shell.is_valid = False

        self.shells.append(shell)
        self.levels[shell.level].append(shell)
        self.centers[shell.central_atm_grp][shell.level] = shell

    def get_valid_shells(self):
        return [s for s in self.shells if s.is_valid]

    def get_shells_by_identifier(self, identifier, unique_shells=False):
        if unique_shells:
            return [s for s in self.shells if s.identifier == identifier and s.is_valid]
        else:
            return [s for s in self.shells if s.identifier == identifier]

    def get_shells_by_level(self, level, unique_shells=False):
        shells = []

        if level in self.levels:
            if unique_shells:
                shells = [s for s in self.levels[level] if s.is_valid]
            else:
                shells = self.levels[level]
        elif self.verbose:
            logger.warning("The informed level '%d' does not exist." % level)

        return shells

    def get_shells_by_center(self, center, unique_shells=False):
        shells = []

        if center in self.centers:
            if unique_shells:
                shells = [s for s in self.centers[center] if s.is_valid]
            else:
                shells = self.centers[center]
        elif self.verbose:
            logger.warning("The informed center '%s' does not exist." % center)

        return shells

    def get_shell_by_center_and_level(self, center, level, unique_shells=False):
        shell = self.centers.get(center, {}).get(level)

        if shell is None and self.verbose:
            logger.warning("The informed center '%s' does not exist in the level '%d'." % (center, level))

        if unique_shells:
            shell = shell if shell.is_valid else None
            if shell is None and self.verbose:
                logger.warning("The shell found with center '%s' and level '%d' is not unique." % (center, level))

        return shell

    def get_previous_shell(self, center, curr_level, unique_shells=False):
        shell = None

        while curr_level != 0 and shell is None:
            level = curr_level - 1
            shell = self.get_shell_by_center_and_level(center, level, unique_shells)
            curr_level = level

        if shell is None and self.verbose:
            logger.warning("No previous shell centered on '%s' departing from the level '%d' was found."
                           % (center, level))

        return shell

    def get_last_shell(self, center, unique_shells=False):
        shell = None

        shells = self.get_shells_by_center(center, unique_shells)
        if shells:
            shell = sorted(shells, key=int)[-1]
        elif self.verbose:
            logger.warning("No shell centered on '%s' was found." % center)

        return shell

    def get_identifiers(self, level=None, unique_shells=False):
        if level is not None:
            if unique_shells:
                identifiers = [s.identifier for s in get_shells_by_level(level) if s.is_valid]
            else:
                identifiers = [s.identifier for s in get_shells_by_level(level)]
        else:
            if unique_shells:
                identifiers = [s.identifier for s in self.shells if s.is_valid]
            else:
                identifiers = [s.identifier for s in self.shells]

        return sorted(identifiers)

    def to_fingerprint(self, unique_shells=False, fold_to_size=None, count_fp=False):
        indices = self.get_identifiers(unique_shells=unique_shells)
        props = {"num_levels": self.num_levels,
                 "num_bits": self.num_bits,
                 "radius_steps": self.radius_steps}

        if count_fp:
            fp = CountFingerprint(indices, fp_length=self.num_bits, props=props)
        else:
            fp = Fingerprint(indices, fp_length=self.num_bits, props=props)

        if fold_to_size:
            return fp.fold(fold_to_size)
        else:
            return fp

    def _init_controllers(self):
        levels = defaultdict(list)
        centers = defaultdict(lambda: defaultdict(lambda: None))

        for shell in self.shells:
            levels[shell.level].append(shell)
            centers[shell.central_atm_grp][shell.level] = shell

        self.levels = levels
        self.centers = centers


class Shell:

    def __init__(self, central_atm_grp, level, radius, neighborhood=None, inter_tuples=None,
                 np_dtype=np.int64, seed=0, manager=None, is_valid=True, feature_mapper=None):

        self.central_atm_grp = central_atm_grp
        self.level = level
        self.radius = radius
        self.is_valid = is_valid

        if not neighborhood:
            # If inter_tuples was not defined initialize the neighborhood as an empty list.
            # Otherwise, define the neighborhood as the atom groups interacting with the central atom group.
            neighborhood = [] if not inter_tuples else [x[1] for x in self._inter_tuples]
        self._neighborhood = set(neighborhood)

        # Always guarantee that the central atom defines the neighborhood.
        self._neighborhood.add(central_atm_grp)

        if not inter_tuples:
            inter_tuples = []
        self._inter_tuples = set(inter_tuples)

        self._interactions = set([x[0] for x in self._inter_tuples])

        # TODO: FIX IT.
        if not feature_mapper:
            default_dict = {**CHEMICAL_FEATURES_IDS, **INTERACTIONS_IDS}
            # TODO: Built a class for managing the feature maps.
            # TODO: Attribute a unique id for each new feature.
            # feature_mapper = defaultdict(lambda: -1, default_dict)
            feature_mapper = default_dict
        self.feature_mapper = feature_mapper

        self.np_dtype = np_dtype
        self.seed = seed

        self._manager = manager
        self._identifier = self.hash_shell()

    @property
    def neighborhood(self):
        return self._neighborhood

    @property
    def interactions(self):
        return self._interactions

    @property
    def inter_tuples(self):
        return self._inter_tuples

    @property
    def manager(self):
        return self._manager

    @property
    def identifier(self):
        return self._identifier

    @property
    def previous_shell(self):
        shell = self._manager.get_previous_shell(self.central_atm_grp, self.level)
        if shell is None:
            logger.exception("No previous shell centered in '%s' was found." % self.central_atm_grp)
            raise ShellCenterNotFound("No previous shell centered in '%s' was found." % self.central_atm_grp)

        return shell

    def hash_shell(self):
        if self.level == 0:
            data = [self.feature_mapper[cf.format_name()] for cf in self.central_atm_grp.features]

            if len(self.central_atm_grp.atoms) == 1:
                features = set()
                for a in self.central_atm_grp.atoms:
                    for ga in a.atm_grps:
                        if ga != self.central_atm_grp:
                            features.update(ga.features)

                data += [self.feature_mapper[cf.format_name()] for cf in features]
            data.sort()
        else:
            cent_prev_id = self.previous_shell.identifier

            # Initialization of a new feature vector
            data = [(self.level, cent_prev_id)]

            inter_tuples = []
            for (inter, nb_atm_grp) in self._inter_tuples:
                prev_nb_shell = self._manager.get_previous_shell(nb_atm_grp, self.level)
                if prev_nb_shell is None:
                    logger.exception("No previous shell centered in %s was found." % nb_atm_grp)
                    raise ShellCenterNotFound("No previous shell centered in %s was found." % nb_atm_grp)

                # 1st elem: interaction type.
                # 2nd elem: previous identifier of the neighbor atom group;
                inter_tuples.append((self.feature_mapper[inter.type], prev_nb_shell.identifier))

            # Sort the tuples to avoid dependence on the order in which tuples are added.
            sorted_list = sorted(inter_tuples)

            # Join the interaction information to the feature vector.
            data += sorted_list

        np_array = np.array(data, self.np_dtype)

        # TODO: Let user to define hash function
        hashed_shell = mmh3.hash(np_array, self.seed, signed=False)

        return hashed_shell

    def __repr__(self):
        return ("<Shell: level=%d, radius=%d, center=%s, interactions=%d>"
                % (self.level, self.radius, self.central_atm_grp, len(self.interactions)))


class ShellGenerator:

    def __init__(self, num_levels, radius_step, bucket_size=10, seed=0, np_dtype=np.int64,
                 num_bits=DEFAULT_SHELL_NBITS):

        self.num_levels = num_levels
        self.radius_step = radius_step

        self.bucket_size = bucket_size
        self.seed = seed
        self.np_dtype = np_dtype
        self.num_bits = num_bits

    def create_shells(self, neighborhood):
        sm = ShellManager(self.num_levels, self.radius_step, self.num_bits)
        ss = ShellSearch(neighborhood, self.bucket_size)

        neighborhood = set(neighborhood)
        skip_atm_grps = set()
        for level in range(self.num_levels):
            radius = self.radius_step * level

            for atm_grp in neighborhood:
                # Ignore centroids that already reached the limit of possible substructures.
                if atm_grp in skip_atm_grps:
                    continue

                shell = None
                if radius > 0:
                    prev_shell = sm.get_previous_shell(atm_grp, level)
                    if not prev_shell:
                        logger.exception("No previous shell centered in %s was found." % nb_atm_grp)
                        raise ShellCenterNotFound("There are no shells initialized to the atom group '%s'." % atm_grp)

                    prev_atm_grps = prev_shell.neighborhood
                    prev_interactions = prev_shell.interactions

                    nb_atm_grps = set(ss.search(atm_grp.centroid, radius))

                    inter_tuples = set()
                    # For each atom group from the previous shell
                    for prev_atm_grp in prev_atm_grps:
                        for inter in prev_atm_grp.interactions:
                            if inter.get_partner(prev_atm_grp) in nb_atm_grps:

                                new_tuple = (inter, inter.get_partner(prev_atm_grp))

                                # Ignore interactions that already exists in the previous shell only if the previous
                                # source atom group does not correspond to the current one. It avoids duplications
                                # on the list of interactions. For example, without this control, an interaction I
                                # between atom A1 and A2 would appear twice in the list: (I, A1) and (I, A2).
                                # Thus, it keeps only the first interaction that appears while increasing the shell.
                                if inter in prev_interactions and new_tuple not in prev_shell.inter_tuples:
                                    continue

                                inter_tuples.add(new_tuple)

                    # It adds a new shell when there are interactions inside the shell.
                    if inter_tuples:
                        shell_nb = set([x[1] for x in inter_tuples])
                        shell_nb.add(atm_grp)

                        shell = Shell(atm_grp, level, radius, neighborhood=shell_nb, inter_tuples=inter_tuples,
                                      manager=sm, seed=self.seed, np_dtype=self.np_dtype)
                else:
                    shell = Shell(atm_grp, level, radius, manager=sm, seed=self.seed, np_dtype=self.np_dtype)

                if shell:
                    sm.add_shell(shell)
                    last_shell = shell
                else:
                    last_shell = sm.get_last_shell(atm_grp)

                # Evaluate if the limit of possible substructures for the current centroid (atom group) was reached.
                if last_shell:
                    # The limit will be reached when the last shell already contains all interactions
                    # established by the atom groups inside the shell. In this case, expanding the radius
                    # will not result in any new shell because a shell is only created when the atoms inside
                    # the last shell establish interactions with the atom groups found after increasing the radius.
                    all_interactions = tuple(chain.from_iterable([g.interactions for g in last_shell.neighborhood]))
                    # It considers only interactions whose atom groups exist in the neigborhood.
                    valid_interactions = set([i for i in all_interactions
                                             if i.atm_grp1 in neighborhood and i.atm_grp2 in neighborhood])

                    current_interactions = last_shell.interactions

                    if valid_interactions == current_interactions:
                        # If the limit was reached for this centroid, in the next level it can be ignored.
                        skip_atm_grps.add(atm_grp)

            # If all atom groups reached the limit of possible substructures, just leave the loop.
            if len(skip_atm_grps) == len(neighborhood):
                logger.warning("The list of shells cannot be expanded anymore. The maximum number "
                               "of substructures were reached.")
                break

        logger.info("Shells creation finished.")
        logger.info("The last level executed was: %d." % level)
        logger.info("The number of levels defined was: %d." % self.num_levels)
        logger.info("Total number of shells created: %d" % sm.num_shells)
        logger.info("Total number of unique shells created: %d" % sm.num_unique_shells)

        return sm
