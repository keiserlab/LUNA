from itertools import chain
from collections import defaultdict
from enum import Enum
import numpy as np
import mmh3

from luna.version import __version__
from luna.util.exceptions import ShellCenterNotFound
from luna.util.default_values import CHEMICAL_FEATURE_IDS, INTERACTION_IDS
from luna.interaction.fp.fingerprint import (DEFAULT_FP_LENGTH, Fingerprint,
                                             CountFingerprint)
from luna.interaction.fp.type import IFPType
from luna.mol.groups import PseudoAtomGroup, AtomGroupNeighborhood
from luna.mol.features import ChemicalFeature


import logging

logger = logging.getLogger()


PAD_DEFAULT = -9


class CompoundClassIds(Enum):
    """An enumeration of compound classes."""
    HETATM = 1
    RESIDUE = 2
    NUCLEOTIDE = 3
    WATER = 4
    UNKNOWN = 5


class ShellManager:

    """Store and manage :class:`Shell` objects.

    Parameters
    ----------
    num_levels : int
        The maximum number of iterations for fingerprint generation.
    radius_step : float
        The multiplier used to increase shell size at each iteration.
        At iteration 0, shell radius is 0 * ``radius_step``, at iteration 1,
        radius is 1 * ``radius_step``, etc.
    fp_length : int
        The fingerprint length (total number of bits).
    ifp_type : :class:`~luna.interaction.fp.type.IFPType`
        The fingerprint type (EIFP, FIFP, or HIFP).
    shells : iterable of :class:`Shell`, optional
        An initial sequence of :class:`Shell` objects (fingerprint features).
    verbose : bool
        If True, warnings issued during the usage of this `ShellManager` will
        be displayed. The default value is False.

    Attributes
    ----------
    num_levels : int
        The maximum number of iterations for fingerprint generation.
    radius_step : float
        The multiplier used to increase shell size at each iteration.
    fp_length : int
        The fingerprint length (total number of bits).
    ifp_type : :class:`~luna.interaction.fp.type.IFPType`
        The fingerprint type (EIFP, FIFP, or HIFP).
    shells : iterable of :class:`Shell`
        The sequence of shells (fingerprint features).
    verbose : bool
        The verbosity state.
    version : str
        The LUNA's version with which shells were generated.
    levels : dict of {int: list of `Shell`}
        Register shells by level, where keys are levels and values are
        lists of `Shell` objects.

        .. note::
            Levels are 0-indexed. So, the first level is 0, second is 1, etc.
            That means if ``num_levels`` is 5, the last level will be 4.
    centers : dict of dict of {int: `Shell`}
        Register shells by center, where keys are \
        :class:`~luna.mol.groups.AtomGroup` objects and values are dict that
        store all shells generated for that center at each iteration (level).
    """

    def __init__(self,
                 num_levels,
                 radius_step,
                 fp_length,
                 ifp_type,
                 shells=None,
                 verbose=False):
        self.num_levels = num_levels
        self.radius_step = radius_step
        self.fp_length = fp_length
        self.ifp_type = ifp_type

        self.shells = shells or []
        self.verbose = verbose

        self.version = __version__

        self._init_controllers()

    @property
    def num_shells(self):
        """int, read-only: Total number of shells in ``shells``."""
        return len(self.shells)

    @property
    def num_unique_shells(self):
        """int, read-only: Total number of unique shells in ``shells``."""
        return len(self.unique_shells)

    @property
    def unique_shells(self):
        """iterable of `Shell`, read-only: Unique shells. \
        Return the same as :meth:`get_valid_shells`."""
        return self.get_valid_shells()

    def find_similar_shell(self, shell):
        """Find a shell in ``shells`` similar to ``shell``.

        Two shells are similar if they represent the same substructural
        information.

        Parameters
        ----------
        shell : `Shell`

        Returns
        -------
         : `Shell` or None
            Return a similar shell or None if it does not find any.
        """
        for added_shell in self.shells:
            # For each group of similar shells, it will exist only
            # one valid shell.
            if added_shell.is_valid():
                if shell.is_similar(added_shell):
                    return added_shell
        return None

    def add_shell(self, shell):
        """Add a new shell to ``shells``.

        Parameters
        ----------
        shell : `Shell`
        """

        # Any new shell without interactions from level 1 onwards will be
        # automatically considered invalid.
        if shell.level > 0 and len(shell.interactions) == 0:
            shell.valid = False
        else:
            found_shell = self.find_similar_shell(shell)

            if found_shell:
                if found_shell.level < shell.level:
                    shell.valid = False

                elif (found_shell.level == shell.level
                        and found_shell.identifier <= shell.identifier):
                    shell.valid = False

                else:
                    found_shell.valid = False

        self.shells.append(shell)
        self.levels[shell.level].append(shell)
        self.centers[shell.central_atm_grp][shell.level] = shell

    def get_valid_shells(self):
        """Return only valid shells.

        A shell is considered invalid if, by the time it is added in
        ``shells``, there is another shell representing the same
        substructural information. That means this shell is not unique
        and does not contributes to any new information.

        On the other hand, if the shell contributes by adding new information
        to ``shells``, then it will be considered valid and unique. So, the
        first shell of a series of shells containing the same information is
        considered valid and the others invalid.

        Returns
        -------
         : list of `Shell`
        """
        return [s for s in self.shells if s.is_valid()]

    def get_shells_by_identifier(self, identifier, unique_shells=False):
        """Get shells by identifier.

        Parameters
        ----------
        identifier : int
            The shell identifier.
        unique_shells : bool
            If True, return only unique shells.
            Otherwise, return all shells having the identifier
            ``identifier`` (the default).

        Returns
        -------
         : list of `Shell`
        """
        if unique_shells:
            return [s for s in self.shells
                    if s.identifier == identifier and s.is_valid()]
        else:
            return [s for s in self.shells if s.identifier == identifier]

    def get_shells_by_level(self, level, unique_shells=False):
        """Get shells by level (iteration number).

        Parameters
        ----------
        level : int
        unique_shells : bool
            If True, return only unique shells.
            Otherwise, return all shells at level ``level`` (the default).

        Returns
        -------
         : list of `Shell`
        """
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
        """Get shells by center (:class:`~luna.mol.groups.AtomGroup` object).

        Parameters
        ----------
        center : :class:`~luna.mol.groups.AtomGroup`
            The center of a shell, which consists of an
            :class:`~luna.mol.groups.AtomGroup` object.
        unique_shells : bool
            If True, return only unique shells.
            Otherwise, return all shells generated for center
            ``center`` (the default).


        Returns
        -------
         : dict of {int: `Shell`}
            All shells generated for center ``center`` at each iteration (key).
        """
        shells = {}

        if center in self.centers:
            if unique_shells:
                shells = {i: s for i, s in self.centers[center].items()
                          if s.is_valid()}
            else:
                shells = self.centers[center]
        elif self.verbose:
            logger.warning("The informed center '%s' does not exist." % center)

        return shells

    def get_shell_by_center_and_level(self,
                                      center,
                                      level,
                                      unique_shells=False):
        """Get the shell generated for center ``center``
        (:class:`~luna.mol.groups.AtomGroup` object) at level (iteration)
        ``level``.

        Parameters
        ----------
        center : :class:`~luna.mol.groups.AtomGroup`
            The center of a shell, which consists of an
            :class:`~luna.mol.groups.AtomGroup` object.
        level : int
            The target level (iteration).
        unique_shells : bool
            If True, return the `Shell` object if it is unique and None
            otherwise. The default value is False.

        Returns
        -------
         : `Shell` or None
            The shell generated for center ``center`` at level ``level``.
            If the `Shell` object is not unique, return None.
        """
        shell = self.centers.get(center, {}).get(level)

        if shell is None and self.verbose:
            logger.warning("The informed center '%s' does not exist in "
                           "the level '%d'." % (center, level))

        if unique_shells:
            shell = shell if shell.is_valid() else None
            if shell is None and self.verbose:
                logger.warning("The shell found with center '%s' and level "
                               "'%d' is not unique." % (center, level))

        return shell

    def get_previous_shell(self, center, curr_level, unique_shells=False):
        """Get the last shell having center ``center`` that was generated
        before level ``curr_level``. For instance, if the current level
        (iteration) is 5 and the last valid shell generated for center
        :math:`C` was at level 4, then :meth:`get_previous_shell` would
        return that shell at level 4.

        Parameters
        ----------
        center : :class:`~luna.mol.groups.AtomGroup`
            The center of a shell, which consists of an
            :class:`~luna.mol.groups.AtomGroup` object.
        curr_level : int
            The current level (iteration).
        unique_shells : bool
            If True, ignore non-valid shells and go down to inferior
            levels until a valid shell is found. If level 0 was reached
            and no valid shell was found, then return None.
            The default value is False.

        Returns
        -------
         : `Shell` or None
            The first previous valid shell or None if no valid shell was found.
        """
        shell = None
        while curr_level != 0 and shell is None:
            level = curr_level - 1
            shell = self.get_shell_by_center_and_level(center,
                                                       level,
                                                       unique_shells)
            curr_level = level

        if shell is None and self.verbose:
            logger.warning("No previous shell centered on '%s' departing from "
                           "the level '%d' was found." % (center, level))

        return shell

    def get_last_shell(self, center, unique_shells=False):
        """Get the last shell generated for center ``center``.

        Parameters
        ----------
        center : :class:`~luna.mol.groups.AtomGroup`
            The center of a shell, which consists of an
            :class:`~luna.mol.groups.AtomGroup` object.
        unique_shells : bool
            If True, ignore non-valid shells.
            That means shells generated at superior levels may be ignored
            if they are not valid. The default value is False.

        Returns
        -------
         : `Shell` or None
            The last shell generated for center ``center`` or None if
            no valid shell was found.
        """
        shell = None

        shells = self.get_shells_by_center(center, unique_shells)
        if shells:
            # Sort by key and get a shell based on the latest value.
            shell = shells[sorted(shells, key=int)[-1]]
        elif self.verbose:
            logger.warning("No shell centered on '%s' was found." % center)

        return shell

    def get_identifiers(self, level=None, unique_shells=False):
        """Get all shells' identifier.

        Parameters
        ----------
        level : int, optional
            If provided, only return identifiers of shells at level ``level``.
        unique_shells : bool
            If True, ignore identifiers of non-valid shells.
            The default value is False.

        Returns
        -------
         : list of int
        """
        if level is not None:
            if unique_shells:
                identifiers = [s.identifier
                               for s in self.get_shells_by_level(level)
                               if s.is_valid()]
            else:
                identifiers = [s.identifier
                               for s in self.get_shells_by_level(level)]
        else:
            if unique_shells:
                identifiers = [s.identifier for s in self.shells
                               if s.is_valid()]
            else:
                identifiers = [s.identifier for s in self.shells]

        return sorted(identifiers)

    def to_fingerprint(self,
                       fold_to_length=None,
                       count_fp=False,
                       unique_shells=False):
        """Encode shells into an interaction fingerprint.

        Parameters
        ----------
        fold_to_length : int, optional
            If provided, fold the fingerprint to length ``fold_to_length``.
        count_fp : bool
            If True, create a count fingerprint
            (:class:`~luna.interaction.fp.fingerprint.CountFingerprint`).
            Otherwise, return a bit fingerprint
            (:class:`~luna.interaction.fp.fingerprint.Fingerprint`).
        unique_shells : bool
            If True, only unique shells are used to create the fingerprint.
            The default value is False.

        Returns
        -------
         : :class:`~luna.interaction.fp.fingerprint.CountFingerprint` or \
                :class:`~luna.interaction.fp.fingerprint.Fingerprint`
        """
        indices = self.get_identifiers(unique_shells=unique_shells)
        props = {"num_levels": self.num_levels,
                 "fp_length": self.fp_length,
                 "radius_step": self.radius_step}

        if count_fp:
            fp = CountFingerprint(indices,
                                  fp_length=self.fp_length,
                                  props=props)
        else:
            fp = Fingerprint(indices,
                             fp_length=self.fp_length,
                             props=props)

        if fold_to_length:
            return fp.fold(fold_to_length)
        return fp

    def trace_back_feature(self,
                           feature_id, ifp,
                           unique_shells=False):
        """Trace a feature from a fingerprint back to the shells
        that originated that feature.

        .. note::
            Due to fingerprint folding, multiple substructures may end up
            encoded in the same bit, the so-called collision problem.
            So, if the provided feature contains collisions, shells
            representing different substructures may be returned by
            :meth:`trace_back_feature`.

        Parameters
        ----------
        feature_id: int
            The target feature id.
        ifp : `Fingerprint`
            The fingerprint containing the feature ``feature_id``.
        unique_shells : bool
            If True, ignore identifiers of non-valid shells.
            The default value is False.

        Yields
        -------
         : `Shell`

        Examples
        --------

        In the below example, we will assume a LUNA project object named
        ``proj_obj`` already exists. Then, we will generate an EIFP
        fingerprint for the first :class:`~luna.mol.groups.AtomGroupsManager`
        object at ``proj_obj``.

        >>> from luna.interaction.fp.shell import ShellGenerator
        >>> from luna.interaction.fp.type import IFPType
        >>> atm_grps_mngr = list(proj_obj.atm_grps_mngrs)[0]
        >>> num_levels, radius_step = 2, 3
        >>> sg = ShellGenerator(num_levels, radius_step, ifp_type=IFPType.EIFP)
        >>> sm = sg.create_shells(atm_grps_mngr)
        >>> fp = sm.to_fingerprint(fold_to_length=1024, count_fp=True)
        >>> print(fp.indices)
        [   2   19   22   23   34   37   39   45   54   67   71   75   83   84
           93  109  138  140  157  162  181  182  186  187  191  194  206  209
          211  237  246  251  263  271  281  296  304  315  323  358  370  374
          388  392  399  400  419  439  476  481  487  509  519  527  532  578
          587  592  604  605  629  635  645  661  668  698  711  713  732  736
          740  753  764  795  813  815  820  824  825  831  836  855  873  882
          911  926  967  975  976  984  990  996 1020]

        Now, we can trace features back to original identifiers and investigate
        its substructural information.

        >>> ori_indices = list(sm.trace_back_feature(34, fp, \
unique_shells=True))
        >>> print(ori_indices)
        [(494318626, [<Shell: level=0, radius=0.000000, center=<AtomGroup: \
[<ExtendedAtom: 3QQK/0/A/GLN/85/CD>, <ExtendedAtom: 3QQK/0/A/GLN/85/NE2>, \
<ExtendedAtom: 3QQK/0/A/GLN/85/OE1>]>, interactions=0>])]

        """
        for ori_feature in ifp.unfolding_map[feature_id]:
            yield (ori_feature,
                   self.get_shells_by_identifier(ori_feature, unique_shells))

    def _init_controllers(self):
        levels = defaultdict(list)
        centers = defaultdict(lambda: defaultdict(lambda: None))

        for shell in self.shells:
            levels[shell.level].append(shell)
            centers[shell.central_atm_grp][shell.level] = shell

        self.levels = levels
        self.centers = centers


class Shell:

    """ A container to store substructural information, which is the base for
    LUNA fingerprints.

    Shells are centered on an atom or atom group
    (:class:`~luna.mol.groups.AtomGroup` objects) and represent all atoms and
    interactions explicitly within it.

    Parameters
    ----------
    central_atm_grp : :class:`~luna.mol.groups.AtomGroup`
        The shell center.
    level : int
        The level (iteration) at which the shell was generated.
    radius : float
        The shell radius.
    neighborhood : iterable of :class:`~luna.mol.groups.AtomGroup`
        All atoms and atom groups within a shell of radius ``radius`` centered
        on ``central_atm_grp``.
    inter_tuples : iterable of \
            (:class:`~luna.interaction.type.InteractionType`, \
             :class:`~luna.mol.groups.AtomGroup`)
        All interactions within a shell of radius ``radius`` centered on
        ``central_atm_grp``. Each tuple contains an
        :class:`~luna.interaction.type.InteractionType` object and one of the
        :class:`~luna.mol.groups.AtomGroup` objects participating to the
        interaction.

        .. note::
            As an interaction involves two participants, it would be expected
            that each interaction produces two tuples. However, by default,
            ShellGenerator sorts atom groups and considers only the first tuple
            that appears, which guarantees that only one of the possible tuples
            is added to avoid information duplication.
    diff_comp_classes : bool
        If True (the default), include differentiation between
        compound classes. That means structural information originated from
        :class:`~luna.mol.groups.AtomGroup` objects belonging to residues,
        nucleotides,  ligands, or water molecules will be considered different
        even if their structural information are the same. This is useful for
        example to differentiate protein-ligand interactions from
        residue-residue ones.
    dtype : data-type
        Use arrays of type ``dtype`` to store information.
        The default value is np.int64.
    seed : int
        A seed to generate shell identifiers through the MurmurHash3 hash
        function. The default value is 0.
    manager : `ShellManager`
        The `ShellManager` object that stores and controls this `Shell` object.
    valid : bool
        If the shell is valid or not. By default, all shells are considered
        valid.
    feature_mapper : dict, optional
        A dict that maps atoms and interactions to unique values.
        If not provided, ``feature_mapper`` will inherit from the default
        mappings ``CHEMICAL_FEATURE_IDS`` and ``INTERACTION_IDS``.

    Attributes
    ----------
    central_atm_grp : :class:`~luna.mol.groups.AtomGroup`
    level : int
    radius : float
    diff_comp_classes : bool
    dtype : data-type
    seed : int
    valid : bool
    feature_mapper : dict
    """

    def __init__(self, central_atm_grp,
                 level, radius, neighborhood=None, inter_tuples=None,
                 diff_comp_classes=True, dtype=np.int64, seed=0, manager=None,
                 valid=True, feature_mapper=None):

        self.central_atm_grp = central_atm_grp
        self.level = level
        self.radius = radius
        self.valid = valid
        self.diff_comp_classes = diff_comp_classes

        if not neighborhood:
            # If inter_tuples was not defined, initialize the neighborhood
            # as an empty list. Otherwise, define the neighborhood as the
            # atom groups interacting with the central atom group.
            neighborhood = ([] if not inter_tuples
                            else [x[1] for x in self._inter_tuples])
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

        self.dtype = dtype
        self.seed = seed

        self._manager = manager

        # Encode the interactions. It represents the substructures comprised by
        # this shell.
        self._encoded_data = self._encode_interactions()
        # Defined in the hash function
        self._identifier = self.hash_shell()

    @property
    def neighborhood(self):
        """iterable of :class:`~luna.mol.groups.AtomGroup`, read-only: \
                All atoms and atom groups within this shell."""
        return self._neighborhood

    @property
    def interactions(self):
        """iterable of :class:`~luna.interaction.type.InteractionType`, \
                read-only: All interactions within this shell."""
        return self._interactions

    @property
    def inter_tuples(self):
        """iterable of tuple, read-only: Each tuple contains an \
                :class:`~luna.interaction.type.InteractionType` object \
                and one of the :class:`~luna.mol.groups.AtomGroup` objects \
                participating to the interaction."""
        return self._inter_tuples

    @property
    def manager(self):
        """`ShellManager`, read-only: The `ShellManager` object that stores \
        and controls this `Shell` object."""
        return self._manager

    @property
    def identifier(self):
        """int, read-only: This shell identifier, which is generated by \
        hashing its encoded data with a hash function. By default, LUNA \
        uses MurmurHash3 as the hash function."""
        return self._identifier

    @property
    def encoded_data(self):
        """iterable of tuple, read-only: The data encoded in this shell."""
        return self._encoded_data

    @property
    def previous_shell(self):
        """`Shell`, read-only: The previous shell, i.e., a shell centered \
        on the same central :class:`~luna.mol.groups.AtomGroup` object from \
        a previous level. For example, if this shell is in level 5, return \
        a shell from level 4 having the same center."""
        shell = self._manager.get_previous_shell(self.central_atm_grp,
                                                 self.level)
        if shell is None:
            error_msg = ("No previous shell centered in '%s' was found."
                         % self.central_atm_grp)
            logger.exception(error_msg)
            raise ShellCenterNotFound(error_msg)
        return shell

    def is_valid(self):
        """If the shell is valid or not.

        Returns
        -------
         : bool
        """
        return self.valid

    def is_similar(self, shell):
        """ If this shell is similar to ``shell``.

        Two shells are similar if they represent the same substructural
        information.

        Parameters
        ----------
        shell : `Shell`

        Returns
        -------
         : bool
        """

        # The shells' identifiers will be equal when the shells encode the same
        # information and were constructed in the same level.
        if self.level == shell.level:
            if self.identifier == shell.identifier:
                return True
            return False

        # If none of the shells have interactions, check if the shells have
        # equal identifiers.
        if not self.interactions and not shell.interactions:
            # If the identifiers are different, but the shells' central group
            # is the same, it means the shells encode the same information
            # even if their levels are different.
            #
            # OBS: this test is only for fast identification without doing a
            # recursive procedure as the one applied in the ELSE statement.
            if self.central_atm_grp == shell.central_atm_grp:
                return True
            # Although two shells contain different identifiers and their
            # centroids are not the same, they can still be equal if they
            # were obtained in different levels.
            else:
                if self.level > shell.level:
                    return self.previous_shell.is_similar(shell)
                return self.is_similar(shell.previous_shell)

        # If both shells have interactions, check if the interactions are
        # equal.
        elif self.interactions and shell.interactions:
            if self.encoded_data == shell.encoded_data:
                return True
        return False

    def hash_shell(self):
        """Hash this shells' substructural information into a 32-bit integer
        using MurmurHash3.

        Returns
        -------
         : int
            A 32-bit integer representing this shell's substructural
            information.
        """

        if self.level == 0:
            # Get the initial data for the first shell according to the IFP
            # type the user has chosen.
            data = self._initial_shell_data()
        else:
            cent_prev_id = self.previous_shell.identifier

            # Initialization of a new feature vector
            data = [(self.level, cent_prev_id)]

            inter_tuples = []
            for (inter, nb_atm_grp) in self._inter_tuples:
                if isinstance(nb_atm_grp, PseudoAtomGroup):
                    # Pseudo-groups that represents only the atoms involved in
                    # a specific interaction don't generate shells.
                    # Only their parents do, i.e., the whole group generates
                    # shells. So, use it instead.
                    nb_atm_grp = nb_atm_grp.parent_grp

                prev_nb_shell = \
                    self._manager.get_previous_shell(nb_atm_grp, self.level)
                if prev_nb_shell is None:
                    error_msg = ("No previous shell centered in %s was found."
                                 % nb_atm_grp)
                    logger.exception(error_msg)
                    raise ShellCenterNotFound(error_msg)

                # 1st elem: interaction type.
                # 2nd elem: previous identifier of the neighbor atom group;
                inter_tuples.append((self.feature_mapper[inter.type],
                                     prev_nb_shell.identifier))

            # Sort the tuples to avoid dependence on the order in which tuples
            # are added.
            sorted_list = sorted(inter_tuples)

            # Join the interaction information to the feature vector.
            data += sorted_list

        np_array = np.array(data, self.dtype)

        # TODO: Let the user define a hash function
        hashed_shell = mmh3.hash(np_array, self.seed, signed=False)

        return hashed_shell

    def _initial_shell_data(self):
        data = []

        # EIFP uses atomic invariants.
        if self.manager.ifp_type == IFPType.EIFP:
            # Shells use atomic invariants as data. In case of atom groups,
            # the data consists of a list of invariants.
            data = sorted([atm.invariants
                           for atm in self.central_atm_grp.atoms])

        # HIFP uses atomic invariants for atoms and pharmacophore information
        # for atom groups.
        elif self.manager.ifp_type == IFPType.HIFP:

            if len(self.central_atm_grp.atoms) == 1:
                # Shells whose centroid are atoms use invariants as data.
                data = sorted([atm.invariants
                               for atm in self.central_atm_grp.atoms])
            else:
                # Shells whose centroid are atoms' group use pharmacophore
                # as data.
                data = [[self.feature_mapper[cf.format_name()]
                         for cf in self.central_atm_grp.features]]

        # FIFP uses pharmacophore properties for atoms and atoms' group.
        elif self.manager.ifp_type == IFPType.FIFP:

            atm_grp_data = [self.feature_mapper[cf.format_name()]
                            for cf in self.central_atm_grp.features]

            if len(self.central_atm_grp.atoms) == 1:
                features = set()
                # It will loop only through one atom and, therefore, it can be
                # removed without losing information. However, if someday I
                # decide to remove the If, then the code for capturing all atom
                # derived features will be already working.
                for atm in self.central_atm_grp.atoms:
                    for atm_grp in atm.atm_grps:
                        if atm_grp != self.central_atm_grp:
                            features.update(atm_grp.features)

                atm_grp_data += [self.feature_mapper[cf.format_name()]
                                 for cf in features]

            data = [sorted(atm_grp_data)]

        # Include differentiation between compound classes, i.e., groups
        # belonging to Residues, Nucleotides, Ligands, or Waters are
        # treated as being different even when the group would be considered
        # the same.
        if self.diff_comp_classes:
            # `classes` is a list because an atom group can be formed by
            # multiple compounds, which can be from different (e.g., residue
            # and ligand bound covalently) or same (e.g., amide groups formed
            # by backbone residues) classes.
            classes = sorted([CompoundClassIds[c.get_class().upper()].value
                              for c in self.central_atm_grp.compounds])

            np_data = np.array(data, self.dtype)
            np_classes = np.array(classes, self.dtype)

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
            init_src_shell = \
                self._manager.get_previous_shell(inter.src_grp, 1)
            init_trgt_shell = \
                self._manager.get_previous_shell(inter.trgt_grp, 1)
            shell_ids = tuple(sorted([init_src_shell.identifier,
                                      init_trgt_shell.identifier]))
            encoded_data.append((shell_ids, self.feature_mapper[inter.type]))
        return sorted(encoded_data)

    def __repr__(self):
        return ("<Shell: level=%d, radius=%f, center=%s, interactions=%d>"
                % (self.level, self.radius,
                   self.central_atm_grp, len(self.interactions)))


class ShellGenerator:

    """Generate shells, the base information of LUNA fingerprints.

    Parameters
    ----------
    num_levels : int
        The maximum number of iterations for fingerprint generation.
    radius_step : float
        The multiplier used to increase shell size at each iteration.
        At iteration 0, shell radius is 0 * ``radius_step``, at iteration 1,
        radius is 1 * ``radius_step``, etc.
    fp_length : int
        The fingerprint length (total number of bits).
        The default value is :math:`2^{32}`.
    ifp_type : :class:`~luna.interaction.fp.type.IFPType`
        The fingerprint type (EIFP, FIFP, or HIFP). The default value is EIFP.
    diff_comp_classes : bool
        If True (the default), include differentiation between compound
        classes. That means structural information originated from
        :class:`~luna.mol.groups.AtomGroup` objects belonging to residues,
        nucleotides,  ligands, or water molecules will be considered different
        even if their structural information are the same. This is useful for
        example to differentiate protein-ligand interactions from
        residue-residue ones.
    dtype : data-type
        Use arrays of type ``dtype`` to store information.
        The default value is np.int64.
    seed : int
        A seed to generate shell identifiers through the MurmurHash3
        hash function. The default value is 0.
    bucket_size : int
        Bucket size of KD tree.
        You can play around with this to optimize speed if you feel like it.
        The default value is 10.

    Attributes
    ----------
    num_levels : int
    radius_step : float
    fp_length : int
    ifp_type : :class:`~luna.interaction.fp.type.IFPType`
    diff_comp_classes : bool
    dtype : data-type
    seed : int
    bucket_size : int

    Examples
    --------

    In the below example, we will assume a LUNA project object named
    ``proj_obj`` already exists.

    First, let's define a `ShellGenerator` object that will create shells over
    2 iterations (levels). At each iteration, the shell radius will be
    increased by 3 and substructural information will be encoded following EIFP
    definition.

    >>> from luna.interaction.fp.shell import ShellGenerator
    >>> from luna.interaction.fp.type import IFPType
    >>> num_levels, radius_step = 2, 3
    >>> sg = ShellGenerator(num_levels, radius_step, ifp_type=IFPType.EIFP)

    After defining the generator, we can create shells by calling
    :meth:`create_shells`, which expects an
    :class:`~luna.mol.groups.AtomGroupsManager` object. In this example, we
    will the first :class:`~luna.mol.groups.AtomGroupsManager` object from an
    existing LUNA project (``proj_obj``).

    >>> atm_grps_mngr = list(proj_obj.atm_grps_mngrs)[0]
    >>> sm = sg.create_shells(atm_grps_mngr)
    >>> print(sm.num_shells)
    528

    Now, with shells stored in the `ShellManager` object you can, for instance:

        * Generate fingerprints:

            >>> fp = sm.to_fingerprint(fold_to_length=1024)
            >>> print(fp.indices)
            [   2   19   22   23   34   37   39   45   54   67   71   75   83   84
               93  109  138  140  157  162  181  182  186  187  191  194  206  209
              211  237  246  251  263  271  281  296  304  315  323  358  370  374
              388  392  399  400  419  439  476  481  487  509  519  527  532  578
              587  592  604  605  629  635  645  661  668  698  711  713  732  736
              740  753  764  795  813  815  820  824  825  831  836  855  873  882
              911  926  967  975  976  984  990  996 1020]

        * Visualize substructural information in Pymol\:

            >>> from luna.interaction.fp.view import ShellViewer
            >>> shell_tuples = [(atm_grps_mngr.entry,
            ...                  sm.unique_shells,
            ...                  proj_obj.pdb_path)]
            >>> sv = ShellViewer()
            >>> sv.new_session(shell_tuples, "example.pse")
    """

    def __init__(self, num_levels, radius_step,
                 fp_length=DEFAULT_FP_LENGTH, ifp_type=IFPType.EIFP,
                 diff_comp_classes=True, dtype=np.int64, seed=0,
                 bucket_size=10):

        self.num_levels = num_levels
        self.radius_step = radius_step
        self.fp_length = fp_length
        self.ifp_type = ifp_type
        self.diff_comp_classes = diff_comp_classes

        self.seed = seed
        self.dtype = dtype
        self.bucket_size = bucket_size

    def create_shells(self, atm_grps_mngr):
        """Perceive substructural information from
        :class:`~luna.mol.groups.AtomGroup` objects and their interactions,
        and represent such information as shells.

        Parameters
        ----------
        atm_grps_mngr : :class:`~luna.mol.groups.AtomGroupsManager`
            Container of :class:`~luna.mol.groups.AtomGroup` objects and their
            interactions.

        Returns
        -------
         : `ShellManager`

        Raises
        ------
        ShellCenterNotFound
            If it fails to recover a shell having a given center.
        """

        sm = ShellManager(self.num_levels,
                          self.radius_step,
                          self.fp_length,
                          self.ifp_type)

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

        # Sort the atom groups to avoid order dependence.
        sorted_neighborhood = sorted(neighborhood)

        level = -1
        for level in range(self.num_levels):
            radius = self.radius_step * level

            for atm_grp in sorted_neighborhood:

                # Ignore centroids that already reached the limit of possible
                # substructures.
                if atm_grp in skip_atm_grps:
                    continue

                # It stores all possible expansions each group can do.
                # Each expansion is a derived group. Initially, the list
                # contains only the groups derived from the central atom
                # group. But, later it may also contain derived groups
                # from interacting partner groups.
                all_derived_atm_grps = \
                    self._get_derived_grps(atm_grp, pseudo_grps_mapping)

                # In level 0, the number of unique derived groups is 0 as the
                # shell initially only contains information of the centroid.
                unique_derived_atm_grps = []

                shell = None

                if radius > 0:
                    prev_shell = sm.get_previous_shell(atm_grp, level)
                    if not prev_shell:
                        error_msg = ("No previous shell centered in %s "
                                     "was found." % atm_grp)
                        logger.exception(error_msg)
                        raise ShellCenterNotFound(error_msg)

                    prev_atm_grps = prev_shell.neighborhood
                    prev_interactions = prev_shell.interactions

                    nb_atm_grps = set(nbs.search(atm_grp.centroid, radius))

                    inter_tuples = set()
                    interactions_to_add = set()

                    # For each atom group from the previous shell.
                    for prev_atm_grp in sorted(prev_atm_grps):
                        # Include the partner's derived groups to the set of
                        # all derived groups.
                        all_derived_atm_grps.update(self._get_derived_grps(prev_atm_grp, pseudo_grps_mapping))

                        for inter in prev_atm_grp.interactions:
                            # Only PseudoAtomGroup objects should deal with
                            # hydrophobic interactions.
                            if isinstance(prev_atm_grp, PseudoAtomGroup) is False and inter.type == "Hydrophobic":
                                continue

                            partner_grp = self._recover_partner_grp(prev_atm_grp, inter, pseudo_grps_mapping)

                            if partner_grp is not None:
                                if inter not in interactions_to_add and partner_grp in nb_atm_grps:
                                    new_tuple = (inter, partner_grp)

                                    # Ignore interactions that already exists
                                    # in the previous shell to avoid
                                    # duplications in the list of interactions.
                                    # For example, without this control, an
                                    # interaction I between atom A1 and A2
                                    # would appear twice in the list: (I, A1)
                                    # and (I, A2). Thus, it keeps only the
                                    # first interaction that appears while
                                    # increasing the shell.
                                    if inter in prev_interactions and new_tuple not in prev_shell.inter_tuples:
                                        continue

                                    inter_tuples.add(new_tuple)
                                    interactions_to_add.add(inter)

                    # Get valid derived groups, which are those ones inside the
                    # current sphere.
                    valid_derived_atm_grps = \
                        set([ag for ag in all_derived_atm_grps
                             if ag in nb_atm_grps])

                    unique_derived_atm_grps = \
                        valid_derived_atm_grps - prev_atm_grps

                    # It adds a new shell only if there are new interactions
                    # and derived atom groups inside the shell.
                    shell_nb = set([x[1] for x in inter_tuples]) | valid_derived_atm_grps
                    shell_nb.add(atm_grp)

                    shell = Shell(atm_grp, level, radius,
                                  neighborhood=shell_nb,
                                  inter_tuples=inter_tuples,
                                  manager=sm,
                                  diff_comp_classes=self.diff_comp_classes,
                                  seed=self.seed, dtype=self.dtype)
                else:
                    shell = Shell(atm_grp, level, radius, manager=sm,
                                  diff_comp_classes=self.diff_comp_classes,
                                  seed=self.seed, dtype=self.dtype)

                if shell:
                    sm.add_shell(shell)
                    last_shell = shell
                else:
                    last_shell = sm.get_last_shell(atm_grp)

                # Evaluate if the limit of possible substructures for the
                # current centroid (atom group) was reached.
                if last_shell:
                    # The limit will be reached when the last shell already
                    # contains all interactions established by the atom
                    # groups inside the shell. In this case, expanding the
                    # radius will not result in any new shell because a shell
                    # is only created when the atoms inside the last shell
                    # establish interactions with the atom groups found after
                    # increasing the radius.
                    all_nb_interactions = set(chain.from_iterable([g.interactions for g in last_shell.neighborhood]))
                    # It considers only interactions whose atom groups exist in
                    # the neighborhood.
                    valid_interactions = set([i for i in all_nb_interactions if i.src_grp in neighborhood and i.trgt_grp in neighborhood])

                    # It identifies if the convergence for the interactions
                    # was reached.
                    interactions_converged = valid_interactions == last_shell.interactions

                    # It identifies if the convergence for the expansion of
                    # atom groups was reached, which happens when all derived
                    # groups were already included in the centroid
                    # neighborhood. However, it may happen that this
                    # requirement was fulfilled by the time new unique derived
                    # groups were included in the centroid neighborhood.
                    # Thus, the second test provides a chance to all of these
                    # recently discovered groups to expand.
                    grp_expansions_converged = all_derived_atm_grps == last_shell.neighborhood and len(unique_derived_atm_grps) == 0

                    # The local convergence is reached when no atom group
                    # inside the current sphere can expand or all of its
                    # interactions were already included to the sphere.
                    local_convergence = \
                        interactions_converged and grp_expansions_converged
                    # The global convergence occurs when all interactions in
                    # the binding site were already included in the sphere.
                    global_convergence = \
                        all_interactions == last_shell.interactions

                    # If the limit was reached for this centroid, in the next
                    # level it can be ignored.
                    if local_convergence or global_convergence:
                        skip_atm_grps.add(atm_grp)

            # If all atom groups reached the limit of possible substructures,
            # just leave the loop.
            if len(skip_atm_grps) == len(neighborhood):
                logger.debug("The list of shells cannot be expanded anymore. "
                             "The maximum number of substructures were "
                             "reached.")
                break

        logger.debug("Shells creation finished.")
        logger.debug("The last level executed was: %d." % level)
        logger.debug("The number of levels defined was: %d." % self.num_levels)
        logger.debug("Total number of shells created: %d" % sm.num_shells)
        logger.debug("Total number of unique shells created: %d"
                     % sm.num_unique_shells)

        return sm

    def _generate_pseudo_grps(self, atm_grps_mngr):
        pseudo_grps = {}
        mapped_interactions = set()

        for hydrop_grp in atm_grps_mngr.filter_by_types(["Hydrophobe"]):
            for inter in hydrop_grp.interactions:

                # Ignore non-hydrophobic interactions or interactions already
                # mapped.
                if inter.type != "Hydrophobic" or inter in mapped_interactions:
                    continue

                src_tuple = (inter.src_grp, tuple(sorted(inter.src_interacting_atms)))
                trgt_tuple = (inter.trgt_grp, tuple(sorted(inter.trgt_interacting_atms)))

                for atm_grp, atms in [src_tuple, trgt_tuple]:
                    # It will get a pseudo-group already created or create
                    # a new one.
                    feats = [ChemicalFeature("Hydrophobe")]
                    pseudo_grp = \
                        pseudo_grps.get(atms,
                                        PseudoAtomGroup(atm_grp, atms, feats))

                    pseudo_grp.add_interactions([inter])

                    pseudo_grps[atms] = pseudo_grp

                mapped_interactions.add(inter)

        return pseudo_grps.values()

    def _get_derived_grps(self, atm_grp, pseudo_grps_mapping):
        # Get derived groups for the informed atom group.
        derived_atm_grps = set([ag
                                for a in atm_grp.atoms
                                for ag in a.atm_grps])

        # Get derived pseudo-groups for the informed atom group.
        derived_atm_grps.update(set([pseudo_grp
                                     for a in atm_grp.atoms
                                     for pseudo_grp in pseudo_grps_mapping.get((a,), [])]))

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
            # Recover the partner groups based on the interacting atoms.
            # Then, it will check which returned pseudo-group corresponds to
            # the interacting atoms, i.e., the one that contains exactly the
            # same atoms. This verification is mainly necessary for
            # pseudo-groups composed only by one atom as such groups tend
            # to belong to more than one pseudo-group.
            for pseudo_grp in pseudo_grps_mapping.get(tuple(partner_atms), []):
                if sorted(pseudo_grp.atoms) == partner_atms:
                    partner_grp = pseudo_grp
                    break

            return partner_grp
        else:
            return interaction.get_partner(atm_grp)
