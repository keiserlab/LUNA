from util.exceptions import (BitsValueError, ShellCenterNotFound)
from util.default_values import (CHEMICAL_FEATURES_IDS, INTERACTIONS_IDS)

from Bio.KDTree import KDTree

from rdkit.DataStructs.cDataStructs import (ExplicitBitVect, SparseBitVect)

from itertools import (chain, product)
from collections import defaultdict

from mol.wrappers.pymol import (PymolWrapper, mybio_to_pymol_selection)

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


class PymolShellViwer:

    def __init__(self, input_file, show_cartoon=False, bg_color="white", pse_export_version="1.8"):
        self.wrapper = PymolWrapper()
        self.input_file = input_file
        self.show_cartoon = show_cartoon
        self.bg_color = bg_color
        self.pse_export_version = pse_export_version

        # TODO:
        # Accept color parameters

    def create_session(self, shells, output_file):
        self.new_view()
        self.set_view(shells)
        self.save_session(output_file)

    def new_view(self):
        self.wrapper.reset_view()

        self.wrapper.set("pse_export_version", self.pse_export_version)
        self.wrapper.set("transparency_mode", 3)
        self.wrapper.run_cmds([("bg_color", {"color": self.bg_color})])

        self.wrapper.load(self.input_file)
        self.wrapper.color_by_element(["all"])
        self.wrapper.hide_all()

        if self.show_cartoon:
            self.wrapper.show([("cartoon", "all")])

    def set_view(self, shells):
        for (s, shell) in enumerate(shells):
            centroid_name = "sphere_%d" % s

            self.wrapper.add_pseudoatom(centroid_name, {"color": "white", "pos": list(shell.central_atm_grp.centroid)})

            self.wrapper.hide([("nonbonded", centroid_name)])
            self.wrapper.show([("sphere", centroid_name)])
            self.wrapper.show([("nb_spheres", centroid_name)])

            self.wrapper.show([("dots", centroid_name)])
            self.wrapper.set("dot_color", "red")

            # Compound view (residue, ligand, etc)
            self.wrapper.show([("sticks", mybio_to_pymol_selection(shell.central_atm_grp.compound))])
            self.wrapper.color([("gray", mybio_to_pymol_selection(shell.central_atm_grp.compound) + " AND elem C")])

            self.wrapper.set("sphere_scale", shell.radius, {"selection": centroid_name})
            self.wrapper.set("sphere_transparency", 0.7, {"selection": centroid_name})
            self.wrapper.run_cmds([("center", {"selection": centroid_name})])

            for (i, inter) in enumerate(shell.interactions):
                obj1_name = "Group_%s" % hash(inter.comp1)
                obj2_name = "Group_%s" % hash(inter.comp2)

                self.wrapper.show([("sticks", mybio_to_pymol_selection(inter.comp1.compound)),
                                   ("sticks", mybio_to_pymol_selection(inter.comp2.compound))])

                self.wrapper.add_pseudoatom(obj1_name, {"vdw": 1, "pos": list(inter.comp1.centroid)})
                self.wrapper.add_pseudoatom(obj2_name, {"vdw": 1, "pos": list(inter.comp2.centroid)})

                self.wrapper.hide([("nonbonded", obj1_name), ("nonbonded", obj2_name)])
                self.wrapper.show([("sphere", obj1_name), ("sphere", obj2_name)])
                self.wrapper.set("sphere_scale", 0.4, {"selection": obj1_name})
                self.wrapper.set("sphere_scale", 0.4, {"selection": obj2_name})

                self.wrapper.distance("inter_%s_%d" % (s, i), obj1_name, obj2_name)

                # TODO: Remove
                if inter.type != "Proximal":
                    if shell.level == 1:
                        self.wrapper.color([("orange", "inter_%s_%d" % (s, i))])
                    if shell.level == 2:
                        self.wrapper.color([("red", "inter_%s_%d" % (s, i))])
                    if shell.level == 3:
                        self.wrapper.color([("blue", "inter_%s_%d" % (s, i))])
                else:
                    self.wrapper.color([("gray", "inter_%s_%d" % (s, i))])

    def save_session(self, output_file):
        self.wrapper.save_session(output_file)


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

    def __init__(self, num_levels, radius_step, num_bits, shells=None, verbose=False):

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
                        # print("Equal shells containing no interactions")
                        found_shell = added_shell
                        break
                # If both shells have interactions, check if the interactions are equal.
                elif added_shell.interactions and shell.interactions:
                    if added_shell.interactions == shell.interactions:
                        # print("Equal shells containing equal interactions")
                        found_shell = added_shell
                        break

        return found_shell

    def add_shell(self, shell):
        # print()
        # print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
        found_shell = self.find_similar_shell(shell)

        if found_shell:
            # print()
            # print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
            # print(found_shell.central_atm_grp, found_shell.central_atm_grp.compound)
            # print(found_shell.central_atm_grp.chemicalFeatures)
            # print(found_shell.interactions)
            # print()
            # print(shell.central_atm_grp, shell.central_atm_grp.compound)
            # print(shell.central_atm_grp.chemicalFeatures)
            # print(shell.interactions)
            # print()
            # print(found_shell.identifier)
            # print(shell.identifier)
            # print()
            # print(found_shell._data)
            # print(shell._data)
            # print()

            if found_shell.level < shell.level:
                # print("Old level is lower")
                shell.is_valid = False
            elif ((found_shell.level == shell.level) and
                    (found_shell.identifier <= shell.identifier)):
                # print("Old level is equal to the new one, but the older shell has higher precedence.")
                shell.is_valid = False
            else:
                # print("This shell become invalid. %s" % found_shell.identifier)
                found_shell.is_valid = False

            # print()
            # print(shell.level, found_shell.level)
            # print()
            # if shell.level > 0 and found_shell.level > 0:
            #     filename = get_unique_filename("tmp", size=10)

            #     psv = PymolShellViwer("3QQK.pdb")
            #     shells_to_plot = [found_shell]
            #     output_file1 = "%s_Level-%d_Sphere-%s-FOUND" % (filename, found_shell.level, found_shell.identifier)
            #     psv.create_session(shells_to_plot, output_file1)

            #     psv = PymolShellViwer("3QQK.pdb")
            #     shells_to_plot = [shell]
            #     output_file2 = "%s_Level-%d_Sphere-%s-FOUND" % (filename, shell.level, shell.identifier)
            #     psv.create_session(shells_to_plot, output_file2)

            #     print()
            #     print("Check the files %s and %s" % (output_file1, output_file2))

        # print()
        # print("Should I add the new shell sir? %s" % str(shell.is_valid))

        # print("-----------------------------")
        # print()
        # print()

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
            logger.warning("No previous shell centered in '%s' departing from the level %d was found."
                           % (center, level))

        return shell

    def get_identifiers(self, level=None):
        if level:
            identifiers = [s.identifier for s in get_shells_by_level(level)]
        else:
            identifiers = [s.identifier for s in self.shells]

        return sorted(identifiers)

    def _init_controllers(self):
        levels = defaultdict(list)
        centers = defaultdict(lambda: defaultdict(lambda: None))

        for shell in self.shells:
            levels[shell.level].append(shell)
            centers[shell.central_atm_grp][shell.level] = shell

        # TODO: transform everything into a hidden variable _
        self.levels = levels
        self.centers = centers


class Shell:

    def __init__(self, central_atm_grp, level, radius, inter_tuples=None,
                 np_dtype=np.int64, seed=0, manager=None, is_valid=True,
                 feature_mapper=None):

        self.central_atm_grp = central_atm_grp
        self.level = level
        self.radius = radius
        self.is_valid = is_valid

        if not inter_tuples:
            inter_tuples = []
        self._inter_tuples = set(inter_tuples)

        self._interactions = set([x[0] for x in self._inter_tuples])

        self._atom_groups = set([x[1] for x in self._inter_tuples])
        self._atom_groups.add(central_atm_grp)

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
        return self._atom_groups

    @property
    def interactions(self):
        return self._interactions

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
            data = [self.feature_mapper[cf.format_name()] for cf in self.central_atm_grp.chemicalFeatures]

            if len(self.central_atm_grp.atoms) == 1:
                chemicalFeatures = set()
                for a in self.central_atm_grp.atoms:
                    for ga in a.atomGroups:

                        if ga != self.central_atm_grp:
                            chemicalFeatures.update(ga.chemicalFeatures)

                data += [self.feature_mapper[cf.format_name()] for cf in chemicalFeatures]
            data.sort()
        else:
            cent_prev_id = self.previous_shell.identifier

            # Initialization of a new feature vector
            data = [(self.level, cent_prev_id)]

            # print()
            # print("++++++++++++++++++++++++++++++++++=")
            # print(self.central_atm_grp, self.central_atm_grp.compound)
            # print(self.central_atm_grp.chemicalFeatures)

            # print()

            inter_tuples = set()
            for (inter, nb_atm_grp) in self._inter_tuples:
                # print(nb_atm_grp, nb_atm_grp.compound)
                # print(nb_atm_grp.chemicalFeatures)
                # print(inter.get_partner(nb_atm_grp), inter.get_partner(nb_atm_grp).compound)
                # print(inter.type)
                # print()

                prev_nb_shell = self._manager.get_previous_shell(nb_atm_grp, self.level)
                if prev_nb_shell is None:
                    raise ShellCenterNotFound("No previous shell centered in %s was found." % nb_atm_grp)

                # 1st elem: interaction type.
                # 2nd elem: previous identifier of the neighbor atom group;
                inter_tuples.add((self.feature_mapper[inter.type], prev_nb_shell.identifier))

            # print(inter_tuples)


            # for i in self.interactions:
            #     nb_atm_grp = i.comp2 if i.comp1 == self.central_atm_grp else i.comp1
            #     prev_nb_shell = self._manager.get_previous_shell(nb_atm_grp, self.level)
            #     if prev_nb_shell is None:
            #         raise ShellCenterNotFound("No previous shell centered in %s was found." % nb_atm_grp)

            #     prev_nb_id = prev_nb_shell.identifier

            #     # 1st elem: previous identifier of the neighbor atom group;
            #     # 2nd elem: interaction type.
            #     inter_tuples.add((prev_nb_id, self.feature_mapper[i.type]))

            # Sort the tuples to avoid dependence on the order in which tuples are added.
            sorted_list = sorted(inter_tuples)

            # Join the interaction information to the feature vector.
            data += sorted_list

            # print(data)
            # print()
        self._data = data
        np_array = np.array(data, self.np_dtype)

        # TODO: Let user to define hash function
        hashed_shell = mmh3.hash(np_array, self.seed, signed=False)

        return hashed_shell

    def __repr__(self):
        return ("<Shell: level=%d, radius=%d, center=%s, interactions=%d>"
                % (self.level, self.radius, self.central_atm_grp, len(self.interactions)))


class ShellGenerator:

    def __init__(self, num_levels, radius_step, include_proximal=True, bucket_size=10, seed=0,
                 np_dtype=np.int64, num_bits=DEFAULT_NBITS):

        self.num_levels = num_levels
        self.radius_step = radius_step
        self.include_proximal = include_proximal

        self.bucket_size = bucket_size
        self.seed = seed
        self.np_dtype = np_dtype
        self.num_bits = num_bits

    def create_shells(self, neighborhood):
        sm = ShellManager(self.num_levels, self.radius_step, self.num_bits)
        ss = ShellSearch(neighborhood, self.bucket_size)

        for level in range(self.num_levels):
            radius = self.radius_step * level

            print("> Starting Level=%d, Radius=%d" % (level, radius))
            print()

            for atm_grp in neighborhood:
                if radius > 0:
                    nb_atm_grps = set(ss.search(atm_grp.centroid, radius))

                    prev_atm_grps = sm.get_previous_shell(atm_grp, level).neighborhood

                    inter_tuples = set()
                    for prev_atm_grp in prev_atm_grps:
                        # print(prev_atm_grp, prev_atm_grp.compound)
                        # print("Interactions: %d" % len(prev_atm_grp.interactions))
                        # print()
                        for inter in prev_atm_grp.interactions:
                            if inter.get_partner(prev_atm_grp) in nb_atm_grps:
                                # print()
                                # print("-----------")
                                # print(inter.comp1, inter.comp1.compound)
                                # print(inter.comp2, inter.comp2.compound)
                                # print(inter.type)
                                # print("Tuple has: ", inter.get_partner(prev_atm_grp))
                                # print()
                                inter_tuples.add((inter, inter.get_partner(prev_atm_grp)))

                    if inter_tuples:
                        # print("=========================")
                        # print("Neighborhood before: ", prev_atm_grps)
                        # print()

                        shell = Shell(atm_grp, level, radius, inter_tuples=inter_tuples, manager=sm,
                                      seed=self.seed, np_dtype=self.np_dtype)

                        # print("Neighborhood after: ", shell.neighborhood)
                        # print("-----------------------------------------------------")
                        # print()
                        # print()
                        sm.add_shell(shell)

                        # Plot Spheres
                        psv = PymolShellViwer("3QQK.pdb")
                        filename = get_unique_filename("tmp")
                        shells_to_plot = list(sm.get_shells_by_center(atm_grp).values())
                        output_file = "%s_level-%d_Sphere-%d" % (filename, level, shell.identifier)
                        psv.create_session(shells_to_plot, output_file)
                        del(psv)
                    # exit()
                else:
                    shell = Shell(atm_grp, level, radius, manager=sm, seed=self.seed, np_dtype=self.np_dtype)
                    sm.add_shell(shell)

            # print()

        print()
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print()
        print("Number of shells: %d" % sm.num_shells)
        print("Number of unique shells: %d" % sm.num_unique_shells)
        print("Number of invalid shells: %d" % len([s for s in sm.shells if s.is_valid == False]))

        # psv = PymolShellViwer("3QQK.pdb")
        # for s in sm.unique_shells:
        #     filename = get_unique_filename("tmp")
        #     try:
        #         shells_to_plot = list(sm.get_shells_by_center(atm_grp).values())
        #         output_file = "%s_level-%d_Sphere-%d" % (filename, level, shell.identifier)
        #         psv.create_session(shells_to_plot, output_file)
        #     except Exception:
        #         print(filename, level, shell.identifier)
        #         raise
        exit()

    # def create_shells(self, atm_grps, trgt_grps=None):
    #     sm = ShellManager(shell_nbits=self.num_bits)
    #     ss = ShellSearch(atm_grps, self.bucket_size)

    #     for level in range(self.num_levels):
    #         radius = self.radius_step * level

    #         print("> Starting Level=%d, Radius=%d" % (level, radius))
    #         print()

    #         for atm_grp in atm_grps:
    #             if radius > 0:
    #                 print(">>>>>>> New group")
    #                 print(atm_grp, atm_grp.compound)
    #                 print("-------------------------------------")

    #                 nb_atm_grps = ss.search(atm_grp.centroid, radius)
    #                 interactions = self._get_interactions(atm_grp, nb_atm_grps, trgt_grps)

    #                 print("Interactions: %d" % len(interactions))

    #                 if interactions:
    #                     shell = Shell(atm_grp, level, radius, interactions=interactions, manager=sm,
    #                                   seed=self.seed, np_dtype=self.np_dtype)

    #                     print("Len: %d" % len(interactions))
    #                     print()
    #                     teste = sm.find_equal_shell(shell)
    #                     print(teste)

    #                     sm.add_shell(shell)

    #                 print("######################\n")
    #             else:
    #                 shell = Shell(atm_grp, level, radius, manager=sm, seed=self.seed, np_dtype=self.np_dtype)
    #                 sm.add_shell(shell)

    #     return sm

    # def _get_interactions(self, central_atm_grp, nb_atm_grps, trgt_grps=None):
    #     all_interactions = set()

    #     if trgt_grps is not None:
    #         trgt_grps = set(trgt_grps)

    #     for nb_atm_grp in nb_atm_grps:
    #         if nb_atm_grp == central_atm_grp:
    #             continue

    #         if trgt_grps is not None:
    #             if central_atm_grp not in trgt_grps and nb_atm_grp not in trgt_grps:
    #                 continue

    #         print("NB", nb_atm_grp, nb_atm_grp.compound)

    #         interactions = central_atm_grp.get_interactions_with(nb_atm_grp)
    #         print("NCI: %d" % len(interactions))

    #         if self.include_proximal:
    #             iType = InteractionType(central_atm_grp, nb_atm_grp, DEFAULT_PROXIMAL_INTERACTION_LABEL)
    #             interactions.append(iType)

    #         print("NCI + Prox: %d" % len(interactions))
    #         print()

    #         all_interactions.update(interactions)

    #     return all_interactions


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
