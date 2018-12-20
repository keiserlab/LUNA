from util.exceptions import (BitsValueError, ShellCenterNotFound, PymolSessionNotInitialized)
from util.default_values import (CHEMICAL_FEATURES_IDS, INTERACTIONS_IDS)

from mol.interaction.fp.shell_viewer import PymolShellViewer

from Bio.KDTree import KDTree

from rdkit.DataStructs.cDataStructs import (ExplicitBitVect, SparseBitVect)

from itertools import (chain, product)
from collections import defaultdict

from scipy.sparse import (issparse, csr_matrix)

import numpy as np
import mmh3
import logging

# TODO: Remove
from util.file import get_unique_filename

logger = logging.getLogger(__name__)


DEFAULT_NBITS = 2**32
DEFAULT_FP_LENGTH = 1024
DEFAULT_FP_DTYPE = np.bool_


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

        if not shells:
            shells = []
        self.shells = shells

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

    def get_shells_by_level(self, level):
        shells = []

        if level in self.levels:
            shells = self.levels[level]
        elif self.verbose:
            logger.warning("The informed level '%d' does not exist." % level)

        return shells

    def get_shells_by_center(self, center):
        shells = []

        if center in self.centers:
            shells = self.centers[center]
        elif self.verbose:
            logger.warning("The informed center '%s' does not exist." % center)

        return shells

    def get_shell_by_center_and_level(self, center, level):
        shell = self.centers.get(center, {}).get(level)

        if not shell and self.verbose:
            logger.warning("The informed center '%s' does not exist in the level '%d'." % (center, level))

        return shell

    def get_previous_shell(self, center, curr_level):
        shell = None

        while curr_level != 0 and shell is None:
            level = curr_level - 1
            shell = self.get_shell_by_center_and_level(center, level)
            curr_level = level

        if shell is None and self.verbose:
            logger.warning("No previous shell centered on '%s' departing from the level '%d' was found."
                           % (center, level))

        return shell

    def get_last_shell(self, center):
        shell = None

        shells = self.get_shells_by_center(center)
        if shells:
            shell = sorted(shells, key=int)[-1]
        elif self.verbose:
            logger.warning("No shell centered on '%s' was found." % center)

    def get_identifiers(self, unique_shells=True, level=None):
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

            inter_tuples = set()
            for (inter, nb_atm_grp) in self._inter_tuples:
                prev_nb_shell = self._manager.get_previous_shell(nb_atm_grp, self.level)
                if prev_nb_shell is None:
                    raise ShellCenterNotFound("No previous shell centered in %s was found." % nb_atm_grp)

                # 1st elem: interaction type.
                # 2nd elem: previous identifier of the neighbor atom group;
                inter_tuples.add((self.feature_mapper[inter.type], prev_nb_shell.identifier))

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

    def __init__(self, num_levels, radius_step, implicit_proximal_inter=False, bucket_size=10, seed=0,
                 np_dtype=np.int64, num_bits=DEFAULT_NBITS):

        self.num_levels = num_levels
        self.radius_step = radius_step
        self.implicit_proximal_inter = implicit_proximal_inter

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

                    # It adds a new shell when there are interactions inside the shell or if implicit proximal
                    # interactions were set on. The latter parameter causes shells containing atoms with no interaction
                    # or that only interact between themselves to produce a shell. It allows the algorithm to find
                    # patterns involving disconnected graphs.
                    if inter_tuples or self.implicit_proximal_inter:
                        # If implicit proximal interactions are set on, the shell neighborhood will always be the atoms
                        # inside a sphere of radius R centered on atom A.
                        if self.implicit_proximal_inter:
                            shell_nb = nb_atm_grps
                        else:
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
                    # If implicit proximal interactions were set on, the limit will be reached when the last shell
                    # comprises all the atom groups provided as parameter (variable neighborhood)
                    if self.implicit_proximal_inter:
                        if len(last_shell.neighborhood) == len(neighborhood):
                            # print("Last shell contains all atom groups and cannot expand anymore.")
                            # If the limit was reached for this centroid, in the next level it can be ignored.
                            skip_atm_grps.add(atm_grp)
                    # Otherwise, the limit will be reached when the last shell already contains all interactions
                    # established by the atom groups inside the last shell. In this case, expanding the radius
                    # will not result in any new shell because a shell is only created when the atoms inside
                    # the last shell establish interactions with the atom groups found after increasing the radius.
                    else:
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


class Fingerprint:

    def __init__(self, indices, fp_length, counts=None, unfolding_map=None, name=None):

        indices = np.asarray(indices, dtype=np.long)

        if np.any(indices >= fp_length):
            raise BitsValueError("Provided indices are in a different bit scale.")

        self._indices = np.unique(indices)
        self._counts = counts

        self.fp_length = fp_length
        self.unfolding_map = unfolding_map

        self._name = name

    @classmethod
    def from_bitstring(cls, bitstring, level=-1, **kwargs):
        pass

    @classmethod
    def from_fingerprint(cls, fp, **kwargs):
        pass

    @classmethod
    def from_rdkit(cls, rdkit_fprint, **kwargs):
        pass

    @property
    def indices(self):
        return self._indices

    @property
    def bit_count(self):
        return self.indices.shape[0]

    @property
    def density(self):
        return self.bit_count / self.fp_length

    @property
    def counts(self):
        if self._counts is None:
            self._counts = dict([(k, 1) for k in self.indices])

        return self._counts

    @property
    def name(self):
        if self._name is None:
            return ""
        else:
            return self._name

    def to_vector(self, compressed=True, dtype=DEFAULT_FP_DTYPE):

        data = [self.counts[i] for i in self.indices]
        if compressed:
            try:
                row = np.zeros(self.bit_count)
                col = self.indices
                vector = csr_matrix((data, (row, col)), shape=(1, self.fp_length), dtype=dtype)
            except ValueError:
                raise BitsValueError("Sparse matrix construction failed. Invalid indices or input data.")
        else:
            vector = np.zeros(self.fp_length, dtype=dtype)
            try:
                vector[self.indices] = data
            except IndexError:
                raise BitsValueError("Some of the provided indices are greater than the fingerprint length.")

        return vector

    def to_bit_vector(self, compressed=True):
        return self.to_vector(compressed=compressed, dtype=DEFAULT_FP_DTYPE)

    def to_bit_string(self):
        bit_vector = self.to_bit_vector(compressed=False).astype(np.int)
        return "".join(map(str, bit_vector))

    def to_rdkit(self, rdkit_fp_cls=None):

        if rdkit_fp_cls is None:
            # Classes to store explicit bit vectors: ExplicitBitVect or SparseBitVect.
            # ExplicitBitVect is most useful for situations where the size of the vector is
            # relatively small (tens of thousands or smaller).
            # For larger vectors, use the _SparseBitVect_ class instead.
            if self.fp_length < 1e5:
                rdkit_fp_cls = ExplicitBitVect
            else:
                rdkit_fp_cls = SparseBitVect

        # RDKit data structure defines fingerprints as a std:set composed of ints (signed int).
        # Since we always have values higher than 0 and since the data structure contains only signed ints,
        # then the max length for a RDKit fingerprint is 2^31 - 1.
        # C signed int (32 bit) ranges: [-2^31, 2^31-1].
        max_rdkit_fp_length = 2**31 - 1
        fp_length = min(self.fp_length, max_rdkit_fp_length)
        indices = self.indices % max_rdkit_fp_length

        rdkit_fp = rdkit_fp_cls(fp_length)
        rdkit_fp.SetBitsFromList(indices.tolist())
        return rdkit_fp

    def fold(self, new_fp_length=DEFAULT_FP_LENGTH):

        if new_fp_length > self.fp_length:
            raise BitsValueError("Fold operation requires the new fingerprint length (%d) "
                                 "to be greater than the current one (%d)." % (new_fp_length, self.fp_length))

        if not np.log2(self.fp_length / new_fp_length).is_integer():
            raise BitsValueError("It is not possible to fold the current fingerprint into the informed new length."
                                 "The current length divided by the new one is not a power of 2 number.")

        folded_indices = self.indices % new_fp_length

        unfolding_map = defaultdict(set)
        for k, v in sorted(zip(folded_indices, self.indices)):
            unfolding_map[k].add(v)

        new_fp = self.__class__(folded_indices, new_fp_length, unfolding_map=unfolding_map)

        return new_fp

    def unfold(self):
        indices = []

        if self.unfolding_map is None:
            logger.warning("This fingerprint was not previously folded.")
        else:
            indices = list(chain.from_iterable(self.unfolding_map.values()))

        return np.array(indices)

    def __repr__(self):
        return ("<%s: indices=%s bits=%d name='%s'>" %
                (self.__class__.__name__, repr(self.indices).replace('\n', '').replace(' ', ''),
                 self.fp_length, self.name))
