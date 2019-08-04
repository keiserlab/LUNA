from itertools import chain, product
from collections import defaultdict
import numpy as np
import mmh3


from util.exceptions import ShellCenterNotFound
from util.default_values import CHEMICAL_FEATURE_IDS, INTERACTION_IDS
from mol.interaction.fp.fingerprint import DEFAULT_SHELL_NBITS, Fingerprint, CountFingerprint
from mol.groups import PseudoAtomGroup, AtomGroupNeighborhood
from mol.features import ChemicalFeature

import logging

logger = logging.getLogger()


LEVEL_TO_PRINT = 1000


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

        if shell.level > LEVEL_TO_PRINT:
            print("~~~~~~~~~~~~~~~~~~")
            print("Searching for similar shells...")
            print()
            print("Current shell: ", sorted(shell.central_atm_grp.atoms))
            print("   Level: ", shell.level)
            print("   Id: ", shell.identifier)
            print()

        for added_shell in self.shells:
            # For each group of similar shells, it will exist only one valid shell.
            if added_shell.is_valid():
                if shell.is_similar(added_shell, "\t"):

                    if shell.level > LEVEL_TO_PRINT:
                        print("\t>>>>>>>>>> FOUND A SIMILAR SHELL <<<<<<<<<<")
                        print()

                    return added_shell

        return None

    def add_shell(self, shell):
        # Any new shell without interactions from level 1 onwards will be automatically considered invalid.
        if shell.level > 0 and len(shell.interactions) == 0:
            shell.valid = False

            if shell.level > LEVEL_TO_PRINT:
                print()
                print("\t@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
                print("\t >>> INVALID SHELL - no interaction information...")
                print()
        else:
            found_shell = self.find_similar_shell(shell)

            if shell.level > LEVEL_TO_PRINT:
                if shell.level > LEVEL_TO_PRINT:
                    print()
                    print("\t@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
                if found_shell:
                    print("\tFound shell: ", sorted(found_shell.central_atm_grp.atoms))
                    print("\t  Level: ", found_shell.level)
                    print("\t  Id: ", found_shell.identifier)
                else:
                    print("\tNo shell was found...")
                print()

            if found_shell:
                if found_shell.level < shell.level:
                    shell.valid = False

                    if shell.level > LEVEL_TO_PRINT:
                        print("\t\t>> Found shell level is lower than the current one.")

                elif found_shell.level == shell.level and found_shell.identifier <= shell.identifier:
                    shell.valid = False

                    if shell.level > LEVEL_TO_PRINT:
                        print("\t\t>> Levels are equal but the found shell id is lower or equal to the current one.")
                else:
                    found_shell.valid = False

                    if shell.level > LEVEL_TO_PRINT:
                        print("\t\t>> Neither conditions accepted.")

                if shell.level > LEVEL_TO_PRINT:
                    print()
                    print("\t\t     -> Current shell: ", shell.valid)
                    print("\t\t     -> Found shell: ", found_shell.valid)
                    print()

        self.shells.append(shell)
        self.levels[shell.level].append(shell)
        self.centers[shell.central_atm_grp][shell.level] = shell

    def get_valid_shells(self):
        return [s for s in self.shells if s.is_valid()]

    def get_shells_by_identifier(self, identifier, unique_shells=False):
        if unique_shells:
            return [s for s in self.shells if s.identifier == identifier and s.is_valid()]
        else:
            return [s for s in self.shells if s.identifier == identifier]

    def get_shells_by_level(self, level, unique_shells=False):
        shells = []

        if level in self.levels:
            if unique_shells:
                shells = [s for s in self.levels[level] if s.is_valid()]
            else:
                shells = self.levels[level]
        elif self.verbose:
            logger.warning("The informed level '%d' does not exist." % level)

        return shells

    def get_shells_by_center(self, center, unique_shells=False):
        shells = {}

        if center in self.centers:
            if unique_shells:
                shells = {i: s for i, s in self.centers[center].items() if s.is_valid()}
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
            shell = shell if shell.is_valid() else None
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
            # Sort by key and get a shell based on the latest value.
            shell = shells[sorted(shells, key=int)[-1]]
        elif self.verbose:
            logger.warning("No shell centered on '%s' was found." % center)

        return shell

    def get_identifiers(self, level=None, unique_shells=False):
        if level is not None:
            if unique_shells:
                identifiers = [s.identifier for s in get_shells_by_level(level) if s.is_valid()]
            else:
                identifiers = [s.identifier for s in get_shells_by_level(level)]
        else:
            if unique_shells:
                identifiers = [s.identifier for s in self.shells if s.is_valid()]
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
                 np_dtype=np.int64, seed=0, manager=None, valid=True, feature_mapper=None):

        self.central_atm_grp = central_atm_grp
        self.level = level
        self.radius = radius
        self.valid = valid

        if not neighborhood:
            # If inter_tuples was not defined initialize the neighborhood as an empty list.
            # Otherwise, define the neighborhood as the atom groups interacting with the central atom group.
            neighborhood = [] if not inter_tuples else [x[1] for x in self._inter_tuples]
        self._neighborhood = set(neighborhood)

        # Always add the central atom in its own neighborhood.
        self._neighborhood.add(central_atm_grp)

        inter_tuples = inter_tuples or []
        self._inter_tuples = set(inter_tuples)

        self._interactions = set([x[0] for x in self._inter_tuples])

        # TODO: FIX IT.
        if not feature_mapper:
            default_dict = {**CHEMICAL_FEATURE_IDS, **INTERACTION_IDS}
            # TODO: Built a class for managing the feature maps.
            # TODO: Attribute a unique id for each new feature.
            # feature_mapper = defaultdict(lambda: -1, default_dict)
            feature_mapper = default_dict
        self.feature_mapper = feature_mapper

        self.np_dtype = np_dtype
        self.seed = seed

        self._manager = manager

        # Encode the interactions. It represents the substructures comprised by this shell.
        self._encoded_data = self._encode_interactions()
        # Defined in the hash function
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
    def encoded_data(self):
        return self._encoded_data

    @property
    def previous_shell(self):
        shell = self._manager.get_previous_shell(self.central_atm_grp, self.level)
        if shell is None:
            logger.exception("No previous shell centered in '%s' was found." % self.central_atm_grp)
            raise ShellCenterNotFound("No previous shell centered in '%s' was found." % self.central_atm_grp)

        return shell

    def is_valid(self):
        return self.valid

    def is_similar(self, shell, indentation):

        if self.level > LEVEL_TO_PRINT and shell.level > LEVEL_TO_PRINT:
            print(indentation, "> Potential shell: ", shell.level, shell.radius, sorted(shell.central_atm_grp.atoms))
            print(indentation, "      Level: ", self.level, "VS", shell.level)
            print(indentation, "      Id: ", self.identifier, "VS", shell.identifier)
            print(indentation, "      Is valid? ", self.valid, "VS", shell.valid)
            print()

        # The shells' identifiers will be equal when the shells encode the same information and were constructed in the same level.
        if self.level == shell.level:
            if self.identifier == shell.identifier:
                if self.level > LEVEL_TO_PRINT:
                    print(indentation, indentation, "    >> They have equal identifiers!!")
                    print()

                return True

            if self.level > LEVEL_TO_PRINT:
                print(indentation, indentation, "    >> They have different identifiers!!")
                print()

            return False

        # If none of the shells have interactions, check if the shells have equal identifiers.
        if not self.interactions and not shell.interactions:
            if self.level > LEVEL_TO_PRINT and shell.level > LEVEL_TO_PRINT:
                print(indentation, indentation, "    >> None of them have interactions!!")
                print(indentation, indentation, "          - Do they have the same identifiers?", shell.identifier == self.identifier)
                print(indentation, indentation, "          - Do they have the same central atom?", shell.central_atm_grp == self.central_atm_grp)
                print()

            # If the identifiers are different, but the shells' central group is the same, it means the shells encode the same
            # information even if their levels are different.
            #
            # OBS: this test is only for fast identification without doing a recursive procedure as the one applied in the ELSE statement.
            if self.central_atm_grp == shell.central_atm_grp:
                return True
            # Although two shells contain different identifiers and their centroids are not the same, they can still be equal if they
            # were obtained in different levels.
            else:
                if self.level > shell.level:
                    if self.level > LEVEL_TO_PRINT and shell.level > LEVEL_TO_PRINT:
                        print(indentation, indentation, "    ................. NEW SHELL HAS HIGHER LEVEL .................")
                        print()
                        print(indentation, indentation, "    >> It will look back to the past...")
                        print(indentation, indentation, "           - ", self.previous_shell)
                        print()

                    is_similar = self.previous_shell.is_similar(shell, "\t\t\t\t")

                    if self.level > LEVEL_TO_PRINT and shell.level > LEVEL_TO_PRINT:
                        print(indentation, indentation, "    >> Is similar??? ", is_similar)
                        print()

                    return is_similar
                else:
                    if self.level > LEVEL_TO_PRINT and shell.level > LEVEL_TO_PRINT:
                        print(indentation, indentation, "    ................. OLD SHELL HAS HIGHER LEVEL .................")
                        print()
                        print(indentation, indentation, "    >> It will look back to the past...")
                        print(indentation, indentation, "    >>     - ", shell.previous_shell)
                        print()

                    is_similar = self.is_similar(shell.previous_shell, "\t\t\t\t")

                    if self.level > LEVEL_TO_PRINT:
                        print(indentation, indentation, "    >> Is similar??? ", is_similar)
                        print()
                    return is_similar
        # If both shells have interactions, check if the interactions are equal.
        elif self.interactions and shell.interactions:
            if self.encoded_data == shell.encoded_data:
                if self.level > LEVEL_TO_PRINT:
                    print(indentation, "> Potential shell: ", shell)
                    print(indentation, "      Level: ", self.level, "VS", shell.level)
                    print(indentation, "      Id: ", self.identifier, "VS", shell.identifier)
                    print(indentation, "      Is valid? ", self.valid, "VS", shell.valid)
                    print()
                    print()
                    print(indentation, indentation, " >> ENCODED DATA - New shell: ", self.encoded_data)
                    print()
                    print(indentation, indentation, " >> ENCODED DATA - Old shell: ", shell.encoded_data)
                    print(indentation, indentation, "       - Are they equal?? ", self.encoded_data == shell.encoded_data)
                    print(indentation, indentation, "       - Are they equal?? ", self.interactions == shell.interactions)
                    print(indentation, indentation, "-------------------------------- // ------------------------------------------")
                    print()
                return True
        return False

    def _encode_interactions(self):
        encoded_data = []
        if self.interactions:
            for inter in self.interactions:
                # print()
                # print(inter.src_grp)
                # print(inter.trgt_grp)
                # print(inter.type)
                # print()

                init_src_shell = self._manager.get_previous_shell(inter.src_grp, 1)
                init_trgt_shell = self._manager.get_previous_shell(inter.trgt_grp, 1)
                shell_ids = tuple(sorted([init_src_shell.identifier, init_trgt_shell.identifier]))
                encoded_data.append((shell_ids, self.feature_mapper[inter.type]))

            #     print(init_src_shell)
            #     print("   - ID: ", init_src_shell.identifier)
            #     print(init_trgt_shell)
            #     print("   - ID: ", init_trgt_shell.identifier)
            #     print()
            #     print("   - DATA: ", (tuple(sorted([init_src_shell.identifier, init_trgt_shell.identifier])), self.feature_mapper[inter.type]))
            #     print("--------")
            #     print()
            # exit()

            # for d in encoded_data:
            #     print(d)
            # exit()
        return sorted(encoded_data)

    def hash_shell(self):
        if self.level == 0:
            data = [self.feature_mapper[cf.format_name()] for cf in self.central_atm_grp.features]

            tmp = [cf.name for cf in self.central_atm_grp.features]

            if len(self.central_atm_grp.atoms) == 1:
                features = set()
                # Given the previous If, this For will loop only through one atom and, therefore, it can be removed without
                # losing information. However, if someday I decide to remove the If, then the code for capturing
                # all atom derived features will be already working.
                for atm in self.central_atm_grp.atoms:
                    for atm_grp in atm.atm_grps:
                        if atm_grp != self.central_atm_grp:
                            features.update(atm_grp.features)

                data += [self.feature_mapper[cf.format_name()] for cf in features]

                tmp += [cf.name for cf in features]
            data.sort()

            # if list(self.central_atm_grp.compounds)[0].id[1] == 145 and self.central_atm_grp.atoms[0].name == "N":
            # if len(self.central_atm_grp.atoms) == 1 and self.central_atm_grp.atoms[0].name == "N":
            #     print("--==---==--==---==--==---==--==--")
            #     print(">>> Hashed information:")
            #     print("        - Data: ", data)
            #     print("        - Size: ", len(self.central_atm_grp.atoms))
            #     print("        - Hash: ", hash(self.central_atm_grp.atoms[0]))
            #     print("        - Id: ", id(self.central_atm_grp.atoms[0]))
            #     print("        - Derived groups: ", [ag for atm in self.central_atm_grp.atoms for ag in atm.atm_grps])
            #     print("---------------------------------")
            #     print()

            if self.level > LEVEL_TO_PRINT:
                print(">>> HASHED INFORMATION: ", sorted(tmp))
        else:
            cent_prev_id = self.previous_shell.identifier

            # Initialization of a new feature vector
            data = [(self.level, cent_prev_id)]

            inter_tuples = []
            for (inter, nb_atm_grp) in self._inter_tuples:
                if isinstance(nb_atm_grp, PseudoAtomGroup):
                    # Pseudo-groups that represents only the atoms involved in a specific interaction don't generate shells.
                    # Only their parents do, i.e., the whole group generates shells. So, use it instead.
                    nb_atm_grp = nb_atm_grp.parent_grp

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

        # TODO: Let the user define a hash function
        hashed_shell = mmh3.hash(np_array, self.seed, signed=False)

        # if self.level > LEVEL_TO_PRINT:
        #     print(">>> HASHED INFORMATION: ", data)
        #     print("        - ID: ", hashed_shell)
        #     print()

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

    def create_shells(self, atm_grps_mngr):
        sm = ShellManager(self.num_levels, self.radius_step, self.num_bits)

        all_interactions = atm_grps_mngr.get_all_interactions()

        neighborhood = set(atm_grps_mngr.atm_grps)

        # Create pseudo-groups for hydrophobic interactions.
        pseudo_grps = set(self._generate_pseudo_grps(atm_grps_mngr))

        # Create a neighborhood scope for the searches
        search_scope = neighborhood | pseudo_grps
        nbs = AtomGroupNeighborhood(search_scope, self.bucket_size)

        # Map an atom to all of its pseudo-groups.
        pseudo_grps_mapping = defaultdict(set)
        for pseudo_grp in pseudo_grps:
            atoms = tuple(sorted(pseudo_grp.atoms))
            for atm in atoms:
                # Add an atom to the mapping
                pseudo_grps_mapping[(atm,)].add(pseudo_grp)
            # Add a list of atoms
            pseudo_grps_mapping[atoms].add(pseudo_grp)

        skip_atm_grps = set()

        sorted_neighborhood = sorted(neighborhood)

        for level in range(self.num_levels):
            radius = self.radius_step * level

            for atm_grp in sorted_neighborhood:
                # Ignore centroids that already reached the limit of possible substructures.
                if atm_grp in skip_atm_grps:
                    continue

                # It stores all possible expansions each group can do. Each expansion is a derived group.
                # Initially, the list contains only the groups derived from the central atom group.
                # But, later it may also contain derived groups from interacting partner groups.
                all_derived_atm_grps = self._get_derived_grps(atm_grp, pseudo_grps_mapping)

                # In level 0, the number of unique derived groups is 0 as the shell initially only contain information of the centroid.
                unique_derived_atm_grps = []

                shell = None

                if radius > 0:

                    if level > LEVEL_TO_PRINT:
                        print()
                        print()
                        print("##############################################################")
                        print()
                        print(">>>", sorted(atm_grp.atoms))
                        print(">>>", sorted(atm_grp.features))
                        print(">>>", atm_grp.centroid)
                        print()
                        print("Level: ", level)
                        print("Radius: ", radius)
                        print()

                    prev_shell = sm.get_previous_shell(atm_grp, level)
                    if not prev_shell:
                        logger.exception("No previous shell centered in %s was found." % atm_grp)
                        raise ShellCenterNotFound("There are no shells initialized to the atom group '%s'." % atm_grp)

                    prev_atm_grps = prev_shell.neighborhood
                    prev_interactions = prev_shell.interactions

                    nb_atm_grps = set(nbs.search(atm_grp.centroid, radius))

                    inter_tuples = set()
                    interactions_to_add = set()

                    # if len(atm_grp.features) == 1 and atm_grp.features[0].name == "Amide":
                    # if len(atm_grp.atoms) == 1 and atm_grp.atoms[0].name == "N":
                    #     print("++++++++++++++++++++++++++++")
                    #     print()
                    #     print(prev_atm_grps)
                    #     print()
                    #     print(prev_interactions)
                    #     print()
                    #     print(prev_shell)
                    #     print()
                    #     print(nb_atm_grps)
                    #     print()
                    #     print()
                    #     print("ATOMS: ")
                    #     for atm in atm_grp.atoms:
                    #         print(atm)
                    #         print("        >>", atm.atm_grps)
                    #     print()
                    #     print("Evaluating inside the FOR...")

                    # For each atom group from the previous shell
                    for prev_atm_grp in sorted(prev_atm_grps):

                        # if len(atm_grp.features) == 1 and atm_grp.features[0].name == "Amide":
                        # if len(atm_grp.atoms) == 1 and atm_grp.atoms[0].name == "N":
                        #     print(">> ", prev_atm_grp)
                        #     print("\t     - Derived grps:")
                        #     print("\t        ", self._get_derived_grps(prev_atm_grp, pseudo_grps_mapping))
                        #     print("..................")
                        #     print()

                        # Include the partner's derived groups to the set of all derived groups.
                        all_derived_atm_grps.update(self._get_derived_grps(prev_atm_grp, pseudo_grps_mapping))

                        for inter in prev_atm_grp.interactions:

                            # Only PseudoAtomGroup objects should deal with hydrophobic interactions.
                            if isinstance(prev_atm_grp, PseudoAtomGroup) is False and inter.type == "Hydrophobic":
                                continue

                            partner_grp = self._recover_partner_grp(prev_atm_grp, inter, pseudo_grps_mapping)

                            if partner_grp is not None:
                                if inter not in interactions_to_add and partner_grp in nb_atm_grps:
                                    new_tuple = (inter, partner_grp)

                                    # Ignore interactions that already exists in the previous shell to avoid duplications
                                    # in the list of interactions. For example, without this control, an interaction I
                                    # between atom A1 and A2 would appear twice in the list: (I, A1) and (I, A2).
                                    # Thus, it keeps only the first interaction that appears while increasing the shell.
                                    if inter in prev_interactions and new_tuple not in prev_shell.inter_tuples:
                                        continue

                                    inter_tuples.add(new_tuple)
                                    interactions_to_add.add(inter)

                    # if len(atm_grp.features) == 1 and atm_grp.features[0].name == "Amide":
                    # if len(atm_grp.atoms) == 1 and atm_grp.atoms[0].name == "N":
                    #     print("++++++++++++++++++++++++++++ // ++++++++++++++++++++++++++++++")
                    #     print(">> All derived groups:")
                    #     print("\t", all_derived_atm_grps)
                    #     print()

                    # Get valid derived groups, which are those ones inside the current sphere.
                    valid_derived_atm_grps = set([ag for ag in all_derived_atm_grps if ag in nb_atm_grps])

                    unique_derived_atm_grps = valid_derived_atm_grps - prev_atm_grps

                    unique_interactions = interactions_to_add - prev_interactions

                    # if set([80, 81]) == set([c.id[1] for c in atm_grp.compounds]):
                    # if len(atm_grp.features) == 1 and atm_grp.features[0].name == "Amide":
                    # if len(atm_grp.atoms) == 1 and atm_grp.atoms[0].name == "N":
                    #     print("--------------------------")
                    #     print("Curr level: ", level)
                    #     print("Prev level: ", prev_shell.level)
                    #     print()
                    #     print("# unique inters", len(unique_interactions))
                    #     print("# unique grps", len(unique_derived_atm_grps))
                    #     print(len(unique_interactions or unique_derived_atm_grps))
                    #     print(len(unique_interactions or unique_derived_atm_grps) != 0)
                    #     print()
                    #     print()
                    #     print(" -> All derived groups:")
                    #     print("        - ", all_derived_atm_grps)
                    #     print()
                    #     print(" -> Valid derived groups:")
                    #     print("        - ", valid_derived_atm_grps)
                    #     print()
                    #     print(" -> Previous groups:")
                    #     print("        - ", prev_atm_grps)
                    #     print()

                    if level > LEVEL_TO_PRINT:
                        print("# interactions: ", len(inter_tuples))
                        print("# unique derived groups: ", len(unique_derived_atm_grps))
                        print()

                    # It adds a new shell when there are new interactions and derived atom group inside the shell.
                    if len(unique_interactions or unique_derived_atm_grps) != 0:
                        shell_nb = set([x[1] for x in inter_tuples]) | valid_derived_atm_grps
                        shell_nb.add(atm_grp)

                        shell = Shell(atm_grp, level, radius, neighborhood=shell_nb, inter_tuples=inter_tuples,
                                      manager=sm, seed=self.seed, np_dtype=self.np_dtype)
                    # else:

                    #     if len(atm_grp.features) == 1 and atm_grp.features[0].name == "Amide":
                    #         print(">>> GROUP DIDNT GENERATE A SHELL...")
                    #         print()

                    #         # shell_nb = set([x[1] for x in inter_tuples]) | valid_derived_atm_grps
                    #         # shell_nb.add(atm_grp)

                    #         # tmp_shell = Shell(atm_grp, level, radius, neighborhood=shell_nb, inter_tuples=inter_tuples,
                    #         #                   manager=sm, seed=self.seed, np_dtype=self.np_dtype)
                    #         # print("        - If it has generated a shell it would be: ")
                    #         # print("                ", tmp_shell)
                    #         # print("                ", tmp_shell.identifier)
                    #         print("---------------------------------------------------------------------------------")
                    #         print()

                else:
                    if level > LEVEL_TO_PRINT:
                        print()
                        print()
                        print("##############################################################")
                        print()
                        print(">>>", sorted(atm_grp.atoms))
                        print(">>>", sorted(atm_grp.features))
                        print(">>>", atm_grp.centroid)
                        print()
                        print("Level: ", level)
                        print("Radius: ", radius)
                        print()

                    shell = Shell(atm_grp, level, radius, manager=sm, seed=self.seed, np_dtype=self.np_dtype)

                if level > LEVEL_TO_PRINT:
                    print("........................................")
                    print("New shell...")
                    print()

                    if shell:
                        print(sorted(shell.central_atm_grp.atoms))
                        print("Level", shell.level)
                        print("Id", shell.identifier)
                        print()
                        # print("..... ### .....")
                        # print("Interactions:")
                        # for i in shell.interactions:
                        #     print(i)
                        #     print("    ---> ", hash(i))
                        #     print()
                        # print("...............")
                    else:
                        print(shell)
                    print()

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
                    all_nb_interactions = set(chain.from_iterable([g.interactions for g in last_shell.neighborhood]))
                    # It considers only interactions whose atom groups exist in the neighborhood.
                    valid_interactions = set([i for i in all_nb_interactions if i.src_grp in neighborhood and i.trgt_grp in neighborhood])

                    # It identifies if the convergence for the interactions was reached.
                    interactions_converged = valid_interactions == last_shell.interactions

                    # It identifies if the convergence for the expansion of atom groups was reached, which happens when all derived groups
                    # were already included in the centroid neighborhood. However, it may happen that this requirement was fulfilled
                    # by the time new unique derived groups were included in the centroid neighborhood. Thus, the second test provides
                    # a chance to all of these recently discovered groups to expand.
                    grp_expansions_converged = all_derived_atm_grps == last_shell.neighborhood and len(unique_derived_atm_grps) == 0

                    # The local convergence is reached when no atom group inside the current sphere can expand or all of its interactions
                    # were already included to the sphere.
                    local_convergence = interactions_converged and grp_expansions_converged
                    # The global convergence occurs when all interactions in the binding site were already included in the sphere.
                    global_convergence = all_interactions == last_shell.interactions

                    if level > LEVEL_TO_PRINT:
                        print("Latest shell: ", last_shell.level, last_shell.radius, sorted(last_shell.central_atm_grp.atoms))
                        print("    - ID: ", last_shell.identifier)
                        print()
                        print("   - Limit for interactions was reached?", interactions_converged)
                        print("   - Limit for expansions was reached? Test 1: ", all_derived_atm_grps == last_shell.neighborhood)
                        print("   - Limit for expansions was reached? Test 2 ", len(unique_derived_atm_grps) == 0)
                        print()
                        print("   - LOCAL limit was reached? ", interactions_converged and grp_expansions_converged)
                        print("   - GLOBAL limit was reached?", global_convergence)
                        print()
                        print("# curr Interactions: ", len(last_shell.interactions))
                        print("# valid Interactions: ", len(valid_interactions))
                        print("# all possible Interactions in the NB: ", len(all_nb_interactions))
                        print("# all interactions: ", len(all_interactions))
                        print()
                        print("# curr NB size: ", len(last_shell.neighborhood))
                        print("# max NB size: ", len(all_derived_atm_grps))
                        print()
                        for ag in all_derived_atm_grps:
                            print("-> ", ag, ag.features)
                        print()

                        if level > 0:
                            print("^^^^^^^^^^^^^^^^^^^")
                            for g in sorted(last_shell.neighborhood - prev_atm_grps, key=lambda x: (sorted(x.atoms), sorted(x.feature_names))):
                                print(">>>", sorted(g.atoms), sorted(g.features))
                                print()
                            print()
                            print(">> # NBs", len(last_shell.neighborhood))
                            print()

                    # If the limit was reached for this centroid, in the next level it can be ignored.
                    if local_convergence or global_convergence:

                        if level > LEVEL_TO_PRINT:
                            print("!!!!!! ####### >>>>>>> LIMIT REACHED TO: ", sorted(atm_grp.atoms))
                            print("!!!!!! ####### >>>>>>>     MAX LEVEL: ", level)
                            print()

                        skip_atm_grps.add(atm_grp)

            # If all atom groups reached the limit of possible substructures, just leave the loop.
            if len(skip_atm_grps) == len(neighborhood):
                logger.warning("The list of shells cannot be expanded anymore. The maximum number "
                               "of substructures were reached.")
                break

            # if level == 0:
            #     # exit()
            #     break

        logger.info("Shells creation finished.")
        logger.info("The last level executed was: %d." % level)
        logger.info("The number of levels defined was: %d." % self.num_levels)
        logger.info("Total number of shells created: %d" % sm.num_shells)
        logger.info("Total number of unique shells created: %d" % sm.num_unique_shells)

        # import glob
        # output_file = "tmp/errors/identifiers%s" % (len(glob.glob("tmp/errors/identifiers*")) + 1)
        # print(output_file)
        # with open(output_file, "w") as OUT:
        #     ids = defaultdict(list)
        #     for s in sm.shells:
        #         ids[s.identifier].append(s)

        #     for i in sorted(ids):
        #         for s in ids[i]:
        #             OUT.write("%d,%s" % (i, str((s))))
        #             OUT.write("\n")

        # print("---------------------- FINAL ----------------------")
        # print()
        # print(sm.num_shells)
        # print()
        # exit()

        # clusters = defaultdict(list)

        # for s1 in sorted(sm.shells, key=lambda x: x.identifier):
        #     ref = s1

        #     for s2 in sorted(clusters, key=lambda x: x.identifier):

        #         if s1.is_similar(s2, "\t"):
        #             if s2.level <= ref.level and s2.identifier <= ref.identifier:
        #                 ref = s2
        #                 break
        #     clusters[ref].append(s1)

        # # print("# clusters: ", len(clusters))
        # # print()
        # # exit()

        # for k in sorted(clusters, key=lambda x: x.identifier):
        #     print("Clusters: ", k.identifier)
        #     print(sorted([s.identifier for s in clusters[k]]))
        #     print()
        #     print()

        # exit()

        return sm

    def _generate_pseudo_grps(self, atm_grps_mngr):
        pseudo_grps = {}
        mapped_interactions = set()

        for hydrop_grp in atm_grps_mngr.filter_by_types(["Hydrophobe"]):
            for inter in hydrop_grp.interactions:

                # Ignore non-hydrophobic interactions or interactions already mapped.
                if inter.type != "Hydrophobic" or inter in mapped_interactions:
                    continue

                src_tuple = (inter.src_grp, tuple(sorted(inter.src_interacting_atms)))
                trgt_tuple = (inter.trgt_grp, tuple(sorted(inter.trgt_interacting_atms)))

                for atm_grp, atms in [src_tuple, trgt_tuple]:
                    # It will get a pseudo-group already created or create a new one.
                    pseudo_grp = pseudo_grps.get(atms, PseudoAtomGroup(atm_grp, atms, [ChemicalFeature("Hydrophobe")]))

                    pseudo_grp.add_interactions([inter])

                    pseudo_grps[atms] = pseudo_grp

                mapped_interactions.add(inter)

        return pseudo_grps.values()

    def _get_derived_grps(self, atm_grp, pseudo_grps_mapping):

        # Get derived groups for the informed atom group.
        derived_atm_grps = set([ag for a in atm_grp.atoms for ag in a.atm_grps])

        # Get derived pseudo-groups for the informed atom group.
        derived_atm_grps.update(set([pseudo_grp for a in atm_grp.atoms for pseudo_grp in pseudo_grps_mapping.get((a,), [])]))

        return derived_atm_grps

    def _recover_partner_grp(self, atm_grp, interaction, pseudo_grps_mapping):

        if isinstance(atm_grp, PseudoAtomGroup):
            atoms = sorted(atm_grp.atoms)

            src_atms = sorted(interaction.src_interacting_atms)
            trgt_atms = sorted(interaction.trgt_interacting_atms)

            if atoms == src_atms:
                partner_atms = trgt_atms
            elif atoms == trgt_atms:
                partner_atms = src_atms
            else:
                return None

            partner_grp = None
            # Recover the partner groups based on the interacting atoms. Then, it will check which returned
            # pseudo-group corresponds to the interacting atoms, i.e., the one that contains exactly the same atoms.
            # This verification is mainly necessary for pseudo-groups composed only by one atom as such groups tend
            # to belong to more than one pseudo-group.
            for pseudo_grp in pseudo_grps_mapping.get(tuple(partner_atms), []):
                if sorted(pseudo_grp.atoms) == partner_atms:
                    partner_grp = pseudo_grp
                    break

            return partner_grp
        else:
            return interaction.get_partner(atm_grp)
