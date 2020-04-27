from itertools import chain
from collections import defaultdict
from enum import Enum, auto
import numpy as np
import mmh3

from luna.version import __version__
from luna.util.exceptions import ShellCenterNotFound
from luna.util.default_values import CHEMICAL_FEATURE_IDS, INTERACTION_IDS
from luna.mol.interaction.fp.fingerprint import DEFAULT_FP_LENGTH, Fingerprint, CountFingerprint
from luna.mol.groups import PseudoAtomGroup, AtomGroupNeighborhood
from luna.mol.features import ChemicalFeature

import logging

logger = logging.getLogger()


PAD_DEFAULT = -9


class IFPType(Enum):
    # auto() creates automatic identifiers for each type as this value is not important.
    # Therefore, always remember to compare IFP types using their names and not their values.
    EIFP = auto()
    EIFP_WITH_PHARM_FOR_GROUPS = auto()
    FIFP = auto()


class CompoundClassIds(Enum):
    HETATM = 1
    RESIDUE = 2
    NUCLEOTIDE = 3
    WATER = 4
    UNKNOWN = 5


class ShellManager:

    def __init__(self, num_levels, radius_step, num_bits, ifp_type, shells=None, full_control=True, verbose=False):
        self.num_levels = num_levels
        self.radius_steps = radius_step
        self.num_bits = num_bits
        self.ifp_type = ifp_type

        self.shells = shells or []
        self.verbose = verbose

        self.version = __version__

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

        for added_shell in self.shells:
            # For each group of similar shells, it will exist only one valid shell.
            if added_shell.is_valid():
                if shell.is_similar(added_shell):
                    return added_shell
        return None

    def add_shell(self, shell):
        # Any new shell without interactions from level 1 onwards will be automatically considered invalid.
        if shell.level > 0 and len(shell.interactions) == 0:
            shell.valid = False
        else:
            found_shell = self.find_similar_shell(shell)

            if found_shell:
                if found_shell.level < shell.level:
                    shell.valid = False

                elif found_shell.level == shell.level and found_shell.identifier <= shell.identifier:
                    shell.valid = False

                else:
                    found_shell.valid = False

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
                identifiers = [s.identifier for s in self.get_shells_by_level(level) if s.is_valid()]
            else:
                identifiers = [s.identifier for s in self.get_shells_by_level(level)]
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
                 diff_comp_classes=True, np_dtype=np.int64, seed=0, manager=None,
                 valid=True, feature_mapper=None):

        self.central_atm_grp = central_atm_grp
        self.level = level
        self.radius = radius
        self.valid = valid
        self.diff_comp_classes = diff_comp_classes

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

    def is_similar(self, shell):

        # The shells' identifiers will be equal when the shells encode the same information and were constructed in the same level.
        if self.level == shell.level:
            if self.identifier == shell.identifier:
                return True
            return False

        # If none of the shells have interactions, check if the shells have equal identifiers.
        if not self.interactions and not shell.interactions:
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
                    return self.previous_shell.is_similar(shell)
                return self.is_similar(shell.previous_shell)

        # If both shells have interactions, check if the interactions are equal.
        elif self.interactions and shell.interactions:
            if self.encoded_data == shell.encoded_data:
                return True
        return False

    def hash_shell(self):
        if self.level == 0:
            # Get the initial data for the first shell according to the IFP type the user has chosen.
            data = self._initial_shell_data()
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

        return hashed_shell

    def _initial_shell_data(self):

        data = []

        # EIFP uses atomic invariants.
        if self.manager.ifp_type == IFPType.EIFP:
            # Shells use atomic invariants as data. In case of atom groups, the data consists of a list of invariants.
            data = sorted([atm.invariants for atm in self.central_atm_grp.atoms])

        # EIFP_WITH_PHARM_FOR_GROUPS uses atomic invariants for atoms and pharmacophore information for atom groups.
        elif self.manager.ifp_type == IFPType.EIFP_WITH_PHARM_FOR_GROUPS:

            if len(self.central_atm_grp.atoms) == 1:
                # Shells whose centroid are atoms use invariants as data.
                data = sorted([atm.invariants for atm in self.central_atm_grp.atoms])
            else:
                # Shells whose centroid are atoms' group use pharmacophore as data.
                data = [[self.feature_mapper[cf.format_name()] for cf in self.central_atm_grp.features]]

        # FIFP uses pharmacophore properties for atoms and atoms' group.
        elif self.manager.ifp_type == IFPType.FIFP:

            atm_grp_data = [self.feature_mapper[cf.format_name()] for cf in self.central_atm_grp.features]

            if len(self.central_atm_grp.atoms) == 1:
                features = set()
                # It will loop only through one atom and, therefore, it can be removed without
                # losing information. However, if someday I decide to remove the If, then the code for capturing
                # all atom derived features will be already working.
                for atm in self.central_atm_grp.atoms:
                    for atm_grp in atm.atm_grps:
                        if atm_grp != self.central_atm_grp:
                            features.update(atm_grp.features)

                atm_grp_data += [self.feature_mapper[cf.format_name()] for cf in features]

            data = [sorted(atm_grp_data)]

        #
        #
        # Include differentiation between compound classes, i.e., groups belonging to Residues, Nucleotides, Ligands, or Waters
        # are treated as being different even when the group would be considered the same.
        if self.diff_comp_classes:
            # Classes is a list as multiple classes (so multiple residues) can exist in a group.
            # That can happen, for instance, in amide groups from the backbone.
            classes = sorted([CompoundClassIds[c.get_class().upper()].value for c in self.central_atm_grp.compounds])

            np_data = np.array(data, self.np_dtype)
            np_classes = np.array(classes, self.np_dtype)

            # Add padding to np_classes
            if np_data.shape[1] > np_classes.shape[0]:
                np_classes = np.pad(np_classes, (0, (np_data.shape[1] - np_classes.shape[0])), 'constant', constant_values=PAD_DEFAULT)
            # Add padding to np_data.
            elif np_data.shape[1] < np_classes.shape[0]:
                np_data = np.pad(np_data, ((0, 0), (0, np_classes.shape[0] - np_data.shape[1])), 'constant', constant_values=PAD_DEFAULT)

            data = np.vstack((np_data, np_classes))

        return data

    def _encode_interactions(self):
        encoded_data = []
        for inter in self.interactions:
            init_src_shell = self._manager.get_previous_shell(inter.src_grp, 1)
            init_trgt_shell = self._manager.get_previous_shell(inter.trgt_grp, 1)
            shell_ids = tuple(sorted([init_src_shell.identifier, init_trgt_shell.identifier]))
            encoded_data.append((shell_ids, self.feature_mapper[inter.type]))
        return sorted(encoded_data)

    def __repr__(self):
        return ("<Shell: level=%d, radius=%d, center=%s, interactions=%d>"
                % (self.level, self.radius, self.central_atm_grp, len(self.interactions)))


class ShellGenerator:

    def __init__(self, num_levels, radius_step, diff_comp_classes=True, num_bits=DEFAULT_FP_LENGTH, ifp_type=IFPType.FIFP,
                 bucket_size=10, seed=0, np_dtype=np.int64):

        self.num_levels = num_levels
        self.radius_step = radius_step
        self.diff_comp_classes = diff_comp_classes
        self.num_bits = num_bits
        self.ifp_type = ifp_type

        self.bucket_size = bucket_size
        self.seed = seed
        self.np_dtype = np_dtype

    def create_shells(self, atm_grps_mngr):
        sm = ShellManager(self.num_levels, self.radius_step, self.num_bits, self.ifp_type)

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

        # It controls which atom groups already converged.
        skip_atm_grps = set()

        # Sort the atom groups for avoiding order dependence.
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

                # In level 0, the number of unique derived groups is 0 as the shell initially only contains information of the centroid.
                unique_derived_atm_grps = []

                shell = None

                if radius > 0:
                    prev_shell = sm.get_previous_shell(atm_grp, level)
                    if not prev_shell:
                        logger.exception("No previous shell centered in %s was found." % atm_grp)
                        raise ShellCenterNotFound("There are no shells initialized to the atom group '%s'." % atm_grp)

                    prev_atm_grps = prev_shell.neighborhood
                    prev_interactions = prev_shell.interactions

                    nb_atm_grps = set(nbs.search(atm_grp.centroid, radius))

                    inter_tuples = set()
                    interactions_to_add = set()

                    # For each atom group from the previous shell
                    for prev_atm_grp in sorted(prev_atm_grps):
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

                    # Get valid derived groups, which are those ones inside the current sphere.
                    valid_derived_atm_grps = set([ag for ag in all_derived_atm_grps if ag in nb_atm_grps])

                    unique_derived_atm_grps = valid_derived_atm_grps - prev_atm_grps

                    # It adds a new shell only if there are new interactions and derived atom groups inside the shell.
                    shell_nb = set([x[1] for x in inter_tuples]) | valid_derived_atm_grps
                    shell_nb.add(atm_grp)

                    shell = Shell(atm_grp, level, radius, neighborhood=shell_nb, inter_tuples=inter_tuples,
                                  manager=sm, diff_comp_classes=self.diff_comp_classes,
                                  seed=self.seed, np_dtype=self.np_dtype)
                else:
                    shell = Shell(atm_grp, level, radius, manager=sm, diff_comp_classes=self.diff_comp_classes,
                                  seed=self.seed, np_dtype=self.np_dtype)

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

                    # If the limit was reached for this centroid, in the next level it can be ignored.
                    if local_convergence or global_convergence:
                        skip_atm_grps.add(atm_grp)

            # If all atom groups reached the limit of possible substructures, just leave the loop.
            if len(skip_atm_grps) == len(neighborhood):
                logger.debug("The list of shells cannot be expanded anymore. The maximum number "
                               "of substructures were reached.")
                break

        logger.debug("Shells creation finished.")
        logger.debug("The last level executed was: %d." % level)
        logger.debug("The number of levels defined was: %d." % self.num_levels)
        logger.debug("Total number of shells created: %d" % sm.num_shells)
        logger.debug("Total number of unique shells created: %d" % sm.num_unique_shells)

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
