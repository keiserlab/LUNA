from openbabel import openbabel as ob
from operator import le, ge
from itertools import combinations, product
from collections import defaultdict
import json

from luna.interaction.config import DefaultInteractionConfig, InteractionConfig
from luna.interaction.filter import InteractionFilter
from luna.interaction.type import InteractionType
from luna.mol.features import ChemicalFeature
from luna.wrappers.base import BondType
from luna.analysis.summary import count_interaction_types
import luna.util.math as im
from luna.util.default_values import BOUNDARY_CONFIG
from luna.util.exceptions import IllegalArgumentError
from luna.mol.groups import AtomGroupNeighborhood
from luna.util.file import pickle_data, unpickle_data
from luna.version import __version__

import logging
logger = logging.getLogger()


CATIONS = ("PositivelyIonizable", "PosIonizable", "Positive")
ANIONS = ("NegativelyIonizable", "NegIonizable", "Negative")

COV_BONDS_MAPPING = {
    BondType.SINGLE: "Single bond",
    BondType.DOUBLE: "Double bond",
    BondType.TRIPLE: "Triple bond",
    BondType.AROMATIC: "Aromatic bond"
}

WATER_NAMES = ['HOH', 'DOD', 'WAT', 'H2O', 'OH2']
DEFAULT_LAZY_LIST = WATER_NAMES + ["NH3", "NH4", "CMO", "SCN"]


class InteractionsManager:

    """Store and manage :class:`~luna.interaction.type.InteractionType`
    objects.

    Parameters
    ----------
    interactions : iterable of \
                :class:`~luna.interaction.type.InteractionType`, optional
        An initial sequence of :class:`~luna.interaction.type.InteractionType`
        objects.
    entry : :class:`~luna.mol.entry.Entry`, optional
        The chain or compound used as reference to calculate interactions.

    """

    def __init__(self, interactions=None, entry=None):
        if interactions is None:
            interactions = []

        self.entry = entry
        self._interactions = list(interactions)

        self.version = __version__

    @property
    def interactions(self):
        """ list of :class:`~luna.interaction.type.InteractionType`, \
            read-only: The list of interactions.\
        Additional interactions should be added using the method
        :py:meth:`add_interactions`."""
        return self._interactions

    @property
    def size(self):
        """int, read-only: The number of interactions."""
        return len(self._interactions)

    def get_all_atm_grps(self):
        """Get all atom groups establishing interactions.

        Returns
        -------
         : set of :class:`~luna.mol.groups.AtomGroup`
        """
        atm_grps = set()
        for inter in self.interactions:
            atm_grps.add(inter.src_grp)
            atm_grps.add(inter.trgt_grp)
        return atm_grps

    def count_interations(self, must_have_target=False):
        """Count the number of each type of interaction in ``interactions``.

        Parameters
        ----------
        must_have_target : bool
                If True, count only interactions involving the target ligand.
                The default value is False, which implies all interactions
                will be considered.

        Returns
        -------
         : dict
        """
        return count_interaction_types(self.interactions,
                                       must_have_target=must_have_target)

    def add_interactions(self, interactions):
        """Add one or more :class:`~luna.interaction.type.InteractionType`
        objects to ``interactions``."""

        self._interactions = list(set(self.interactions + list(interactions)))

    def remove_interactions(self, interactions):
        """Remove one or more :class:`~luna.interaction.type.InteractionType`
        objects from ``interactions``.

        Any recursive references to the removed objects will also be cleared.
        """
        self._interactions = list(set(self.interactions) - set(interactions))

        for inter in interactions:
            inter.clear_refs()

    def filter_by_types(self, types):
        """Filter :class:`~luna.interaction.type.InteractionType` objects by
        their types.

        Parameters
        ----------
        types : iterable of str
            A sequence of interaction types.

        Yields
        ------
        :class:`~luna.interaction.type.InteractionType`
        """
        for inter in self.interactions:
            if inter.type in types:
                yield inter

    def filter_out_by_binding_mode(self, binding_modes_filter):
        """Filter out interactions based on binding modes.

        **Note:** this method modifies ``interactions``.

        Parameters
        ----------
        binding_modes_filter : \
                :class:`~luna.interaction.filter.BindingModeFilter`
            A :class:`~luna.interaction.filter.BindingModeFilter` object that
            defines binding mode conditions to decide which interactions
            are valid.

        Returns
        -------
         : set of :class:`~luna.interaction.type.InteractionType`
            The interactions that were filtered out.
        """
        inters_to_remove = set()
        for inter in self.interactions:
            if not binding_modes_filter.is_valid(inter):
                inters_to_remove.add(inter)

        self.remove_interactions(inters_to_remove)

        return inters_to_remove

    def to_csv(self, output_file):
        """Write interactions to a comma-separated values (csv) file.

        Parameters
        ----------
        output_file : str
            The output CSV file.
        """
        interactions_set = set()
        for inter in self.interactions:
            grp1 = ";".join(sorted(["/".join(a.full_atom_name.split("/"))
                                    for a in inter.src_grp.atoms]))
            grp2 = ";".join(sorted(["/".join(a.full_atom_name.split("/"))
                                    for a in inter.trgt_grp.atoms]))

            grp1, grp2 = sorted([grp1, grp2])
            interactions_set.add((grp1, grp2, inter.type))

        with open(output_file, "w") as OUT:
            OUT.write("atom_group1,atom_group2,interaction\n")
            # Sort lines before writing to always keep the same order.
            OUT.write("\n".join([",".join(k)
                                 for k in sorted(interactions_set)]))

    def to_json(self, output_file=None, indent=None):
        """Write interactions to a_initial_shell_data JSON file.

        Parameters
        ----------
        output_file : str
            The output JSON file.
        indent : int or str, optional
            Indent level for pretty-printed JSON files.
            An indent level of 0, negative, or '' only insert newlines.
            Positive integers indent that many spaces per level.
            If a string is provided (e.g., '\\\\t'), it will be used to indent
            each level. The default value is None, which selects the most
            compact representation.
        """
        with open(output_file, 'w') as OUT:
            inter_objs = [inter.as_json() for inter in self.interactions]
            json.dump(inter_objs, OUT, indent=indent)

    def save(self, output_file, compressed=True):
        """Write the pickled representation of the `InteractionsManager` object
        to the file ``output_file``.

        Parameters
        ----------
        output_file : str
            The output file.
        compressed : bool, optional
            If True (the default), compress the pickled representation as a
            gzip file (.gz).

        Raises
        -------
        FileNotCreated
            If the file could not be created.
        """
        pickle_data(self, output_file, compressed)

    @staticmethod
    def load(input_file):
        """Load the pickled representation of an `InteractionsManager` object
        saved at the file ``input_file``.

        Returns
        ----------
         : `InteractionsManager`
            The reconstituted `InteractionsManager` object.

        Raises
        -------
        PKLNotReadError
            If the file could not be loaded.
        """
        return unpickle_data(input_file)

    def __len__(self):
        # Number of interactions
        return self.size

    def __iter__(self):
        """Iterate over children."""
        for inter in self.interactions:
            yield inter


class InteractionCalculator:

    """Calculate interactions.

    .. note::
        This class provides default LUNA methods to calculate interactions.
        However, one can provide their own methods without modifying this
        class. In the **Examples** section, we will show how to define
        custom functions.

    .. note::
        In case you want to disable specific parameters (e.g., angles) used
        during the calculation of interactions, you do not need to define a
        custom function for it. You could just delete the parameter from the
        configuration and LUNA will automatically recognize that a given
        parameter is not necessary anymore.

        Check **Examples 3** to see how to do it and how to implement this
        automatic behavior on your custom functions.


    Parameters
    ----------
    inter_config : :class:`~luna.interaction.config.InteractionConfig`
        An :class:`~luna.interaction.config.InteractionConfig` object with
        all parameters and cutoffs necessary to compute interactions defined
        in ``inter_funcs``. If not provided, the default LUNA configuration
        will be used instead \
        (:class:`~luna.interaction.config.DefaultInteractionConfig`).
    inter_filter : :class:`~luna.interaction.filter.InteractionFilter`, \
            optional
        An :class:`~luna.interaction.filter.InteractionFilter` object to filter
        out interactions on-the-fly. The default value is None, which implies
        no interaction will be filtered out.
    inter_funcs : dict of {tuple : iterable of callable}
        A dict to define custom functions to calculate interactions,
        where keys are tuples of feature names \
        (e.g. ``("Hydrophobic", "Hydrophobic")``) and values are lists of
        references to custom functions (see Examples for more details).
        If not provided, the default LUNA methods will be used instead.
    add_non_cov : bool
         If True (the default), compute non-covalent interactions.
         If you are providing custom functions to compute non-covalent
         interactions and want to make them controllable by this flag, make
         sure to verify the state of ``add_non_cov`` at the beginning of the
         function and return an empty list in case it is False.
    add_cov : bool
        If True (the default), compute covalent interactions.
        If you are providing custom functions to compute covalent interactions
        and want to make them controllable by this flag, make sure to verify
        the state of ``add_cov`` at the beginning of the function and return
        an empty list in case it is False.
    add_proximal : bool
        If True, compute proximal interactions, which are only distance-based
        contacts between atoms or atom groups that, therefore, only imply
        proximity. The default value is False. If you are providing custom
        functions to compute proximal interactions and want to make them
        controllable by this flag, make sure to verify the state of
        ``add_proximal`` at the beginning of the function and return an empty
        list in case it is False.
    add_atom_atom : bool
        If True (the default), compute atom-atom interactions,
        which, as the name suggests, are interactions that only involve atoms
        no matter their features. If you are providing custom functions to
        compute atom-atom interactions and want to make them controllable by
        this flag, make sure to verify the state of ``add_atom_atom`` at the
        beginning of the function and return an empty list in case it is False.

        .. note::
            In LUNA, we consider the following interactions as atom-atom:
            `Van der Waals`, `Van der Waals clash`, and `Atom overlap`.
            We opted to separate `Van der Waals` from other non-covalent
            interactions because LUNA may generate an unnecessary number of
            additional interactions that are usually already represented by
            other non-covalent interactions as weak hydrogen bonds,
            hydrophobic, or dipole-dipole interactions. Thus, to give users
            a fine-grain control over which interactions to calculate, we
            provided this additional flag to turn off the calculation of
            Van der Waals interactions.
    add_dependent_inter : bool
        If True, compute interactions that depend on other interactions.
        Currently, only water-bridged hydrogen bonds and salt bridges have a
        dependency on other interactions. The first, depends on two or more
        hydrogen bonds, while the second depends on an ionic and a hydrogen
        bond. The default value is False, which implies no dependent
        interaction will be computed.
    add_h2o_pairs_with_no_target : bool
        If True, keep interactions of water with atoms and atom groups that
        do not belong to the target of LUNA's analysis, which are chains or
        molecules defined as an :class:`~luna.mol.entry.Entry` instance.
        For example, if the target is a ligand and
        ``add_h2o_pairs_with_no_target`` is False, then water-water and
        water-residue hydrogen bonds will be removed because the ligand is not
        participating in the interactions. The default value is False.
    strict_donor_rules : bool
        If True (the default), hydrogen bonds will only be considered for donor
        atoms with explicit hydrogens bound to them. In that case, angles and
        distances will be evaluated. However, if the molecule containing the
        donor atom is in ``lazy_comps_list``, then angles and hydrogens will be
        ignored and LUNA will proceed with the determination of hydrogen bonds
        based only on donor-acceptor distances. Another exception occurs for
        solvent molecules in which the donor atom is only bound to hydrogens
        atoms (e.g., water, ammonia, and hydrogen sulfide). In that case,
        hydrogens can be positioned in many different ways by Open Babel, which
        may cause LUNA to detect different hydrogen bonds at each run.
        So, to circumvent this problem, by default, LUNA always ignores the
        explicit hydrogen position for donor atoms that only contain hydrogens
        bound to it.
    strict_weak_donor_rules : bool
        If True (the default), weak hydrogen bonds will only be considered for
        donor atoms with explicit hydrogens bound to them. In that case, angles
        and distances will be evaluated.
        The same exceptions described for ``strict_donor_rules`` apply here.
    lazy_comps_list : iterable
         A sequence of molecule names to ignore explicit hydrogen position
         during the calculation of hydrogen bonds and weak hydrogen bonds.
         The default list is
         ['HOH', 'DOD', 'WAT', 'H2O', 'OH2', 'NH3', 'NH4'], which only contains
         water, ammonia, and ammonium ion, including water name variations
         used by different programs.


    Examples
    --------

    **Example 1) How to define custom interactions:**

    In this example, we will define a custom function to calculate
    hydrogen bonds.

    First, let's start importing the classes and the function we will use.

    >>> from luna.interaction.type import InteractionType
    >>> from luna.interaction.calc import InteractionCalculator
    >>> from luna.util.math import euclidean_distance

    Now, we define the custom function, which simply calculates hydrogen bonds
    based on donor-acceptor distances. If it is less than 3.5, then a new
    :class:`~luna.interaction.type.InteractionType` object is created with type
    `Hydrogen bond`.

    .. code-block::

        def custom_hbond_function(self, params):
            if not self.add_non_cov:
                return []

            group1, group2, feat1, feat2 = params
            interactions = []

            cc_dist = euclidean_distance(group1.centroid, group2.centroid)
            if cc_dist <= 3.5:
                params = {"dist_hbond_inter": cc_dist}
                inter = InteractionType(group1, group2, "Hydrogen bond", \
params=params)
                interactions.append(inter)
            return interactions

    .. note::
        Observe that the function checks if ``add_non_cov`` has been turned off
        and if so returns an empty list. That's a recommended strategy because
        it allows one to turn off all non-covalent interactions with a single
        flag.

        Also, observe that `InteractionCalculator` always expects functions to
        return a list at the end. That means multiple interactions may be
        detected for a single pair of :class:`~luna.mol.groups.AtomGroup`
        objects. For example, a donor atom containing 2 hydrogens could,
        in theory, form two different hydrogen bonds with an acceptor atom.


    Now, we have two options to set the custom function to an
    `InteractionCalculator` object:

        1) Define a new dict with the custom functions:

        >>> custom_funcs = {("Donor", "Acceptor"): [custom_hbond_function]}
        >>> ic = InteractionCalculator(inter_funcs=custom_funcs)

        2) Overwrite the default dict in InteractionCalculator:

        >>> ic = InteractionCalculator()
        >>> ic.funcs[("Donor", "Acceptor")] = [custom_hbond_function]

    **Example 2)** How to modify parameters to calculate interactions:

    If you just want to modify specific values from the default configuration,
    you can create a new
    :class:`~luna.interaction.config.DefaultInteractionConfig`,
    alter parameters, and pass it to `InteractionCalculator`.

    >>> from luna.interaction.calc import InteractionCalculator
    >>> from luna.interaction.config import DefaultInteractionConfig
    >>> custom_config = DefaultInteractionConfig()
    >>> custom_config.config["min_dha_ang_hb_inter"] = 120
    >>> ic = InteractionCalculator(inter_config=custom_config)

    Alternatively, you can initiate a new `InteractionCalculator` without
    providing an :class:`~luna.interaction.config.InteractionConfig` object,
    which will cause `InteractionCalculator` to initiate the default
    configuration. Then, you can modify it directly as we did before.

    >>> from luna.interaction.calc import InteractionCalculator
    >>> ic = InteractionCalculator()
    >>> print(ic.inter_config["min_dha_ang_hb_inter"])
    90
    >>> ic.inter_config["min_dha_ang_hb_inter"] = 120
    >>> print(ic.inter_config["min_dha_ang_hb_inter"])
    120

    Finally, if you want to define a custom configuration that will be
    used in your custom functions, you first need to define the
    parameters as a dict and then initiate a new
    :class:`~luna.interaction.config.InteractionConfig`. See below:

    >>> from luna.interaction.config import InteractionConfig
    >>> config = {"param1": 2.5, "param2": 90}
    >>> custom_config = InteractionConfig(config)
    >>> ic = InteractionCalculator(inter_config=custom_config)
    >>> print(ic.inter_config["param1"])
    2.5
    >>> print(ic.inter_config["param2"])
    90

    **Example 3)** How to disable specific parameters and how to enable
    automatic recognition of disabled parameters in custom functions:

    To automatically disable, for instance, verification of angles during the
    calculation of interactions using LUNA's default functions, we just need
    to remove the parameters related to angles from ``inter_config``.
    Let's see an example where we disable angles from hydrogen bonds:

    >>> ic = InteractionCalculator()
    >>> del ic.inter_config["min_dha_ang_hb_inter"]
    >>> del ic.inter_config["min_har_ang_hb_inter"]
    >>> del ic.inter_config["min_dar_ang_hb_inter"]

    This simple behavior is possible thanks to the function
    `is_within_boundary`, which always returns True if the parameter does not
    exist in ``inter_config``. Thus, you can take advantage of this system when
    implementing custom functions so others will also have the possibility to
    turn off specific parameters without modifying the code directly.
    Let's see that in practice.

    First, let's start importing the classes and the function we will use.

    >>> from luna.interaction.config import InteractionConfig
    >>> from luna.interaction.type import InteractionType
    >>> from luna.interaction.calc import InteractionCalculator
    >>> from luna.util.math import euclidean_distance
    >>> from operator import le

    Now, we define a custom function that calculates hydrogen bonds based on
    donor-acceptor distances only. Observe at line #9 that instead of checking
    the cutoff directly, we call `is_within_boundary` with the value,
    parameter, and a comparison function (``le``: less than or equal), which
    will make it possible to modify or even disable the parameter
    automatically.

    .. code-block::
        :linenos:
        :emphasize-lines: 9

        def custom_hbond_function(self, params):
            if not self.add_non_cov:
                return []

            group1, group2, feat1, feat2 = params
            interactions = []

            cc_dist = euclidean_distance(group1.centroid, group2.centroid)
            if self.is_within_boundary(cc_dist, "max_hb_dist", le):
                params = {"dist_hbond_inter": cc_dist}
                inter = InteractionType(group1, group2, "Hydrogen bond", \
params=params)
                interactions.append(inter)
            return interactions

    Finally, we are ready to provide the function and custom parameters
    to `InteractionCalculator`.

    >>> custom_config = InteractionConfig({"max_hb_dist": 3})
    >>> custom_funcs = {("Donor", "Acceptor"): [custom_hbond_function]}
    >>> ic = InteractionCalculator(inter_funcs=custom_funcs, \
inter_config=custom_config)

    By doing so, we can now alter the new parameter or turn it off.

    >>> ic.inter_config["max_hb_dist"] = 3.5
    >>> del ic.inter_config["max_hb_dist"]

    """

    def __init__(self, inter_config=DefaultInteractionConfig(),
                 inter_filter=None, inter_funcs=None, add_non_cov=True,
                 add_cov=True, add_proximal=False, add_atom_atom=True,
                 add_dependent_inter=False,
                 add_h2o_pairs_with_no_target=False,
                 strict_donor_rules=True,
                 strict_weak_donor_rules=True,
                 lazy_comps_list=DEFAULT_LAZY_LIST):

        if (inter_config is not None
                and isinstance(inter_config, InteractionConfig) is False):
            msg = ("The informed interaction configuration "
                   "must be an instance of '%s'." % InteractionConfig)
            raise IllegalArgumentError(msg)

        if (inter_filter is not None
                and isinstance(inter_filter, InteractionFilter) is False):
            msg = ("The informed interaction filter must be "
                   "an instance of '%s'." % InteractionFilter)
            raise IllegalArgumentError(msg)

        self.inter_config = inter_config
        self.inter_filter = inter_filter
        self._inter_funcs = inter_funcs or self._default_functions()

        self.add_non_cov = add_non_cov
        self.add_cov = add_cov
        self.add_proximal = add_proximal
        self.add_atom_atom = add_atom_atom

        self.add_dependent_inter = add_dependent_inter
        self.add_h2o_pairs_with_no_target = add_h2o_pairs_with_no_target

        self.strict_donor_rules = strict_donor_rules
        self.strict_weak_donor_rules = strict_weak_donor_rules
        self.lazy_comps_list = lazy_comps_list or []

    @property
    def funcs(self):
        """dict: The dict that defines functions to calculate interactions."""
        return self._inter_funcs

    @funcs.setter
    def funcs(self, funcs):
        self._inter_funcs = funcs

    def calc_interactions(self, trgt_atm_grps, nb_atm_grps=None):
        """Calculate interactions established by atoms and atoms groups in
        ``trgt_atm_grps`` using methods available in ``funcs``.

        The functions in ``funcs`` are chosen based on the features of each
        atom or atom group. For example, consider that a pair of
        :class:`~luna.mol.groups.AtomGroup` objects have both the features
        'Hydrophobic'. Then, `calc_interactions` will call any interaction
        function defined for the tuple ``("Hydrophobic", "Hydrophobic")`` in
        ``funcs``. Consider now a pair of :class:`~luna.mol.groups.AtomGroup`
        objects whose features are 'Donor' and 'Hydrophobic'. Once again,
        `calc_interactions` will evaluate if there is any function defined for
        the tuple ``("Donor", "Hydrophobic")`` (the order does not matter).
        If there is none, nothing is done and the pair is skipped.

        Parameters
        ----------
        trgt_atm_grps : iterable of :class:`~luna.mol.groups.AtomGroup`
            Compute interactions involving these
            :class:`~luna.mol.groups.AtomGroup` objects.
        nb_atm_grps : iterable of :class:`~luna.mol.groups.AtomGroup`
            If defined, only compute interactions between
            :class:`~luna.mol.groups.AtomGroup` objects from ``trgt_atm_grps``
            with :class:`~luna.mol.groups.AtomGroup` objects from
            ``nb_atm_grps``. If not provided, set ``nb_atm_grps`` to be the
            same as ``trgt_atm_grps``, which implies that interactions will be
            calculated using only pairs of :class:`~luna.mol.groups.AtomGroup`
            objects from ``trgt_atm_grps``.
        """

        # TODO: Water-bridged interaction with weak hydrogen bond
        # TODO: threshold for including slightly out of limit interactions.
        #       For example, a hydrogen bond not included for 0.01A.
        #           Distances: 0.2 and Angles: 5

        # If nb_atm_grps was not informed, it uses the trgt_atm_grps as the
        # neighbors. In this case, the interactions will be target x target.
        nb_comp_grps = nb_atm_grps or trgt_atm_grps

        # Define the scope of the neighborhood search.
        ss = AtomGroupNeighborhood(nb_comp_grps, 10)

        computed_pairs = set()
        all_interactions = []

        bsite_param = "bsite_cutoff"
        bsite_cutoff = self.inter_config.get(bsite_param,
                                             BOUNDARY_CONFIG[bsite_param])

        for trgt_atm_grp in trgt_atm_grps:
            for nb_atm_grp in ss.search(trgt_atm_grp.centroid, bsite_cutoff):

                # It will always ignore interactions involving the same atom
                # groups. Loops in the graph is not permitted and does not make
                # any sense.
                if trgt_atm_grp == nb_atm_grp:
                    continue

                # If the pair has already been calculated.
                if ((trgt_atm_grp, nb_atm_grp) in computed_pairs
                        or (nb_atm_grp, trgt_atm_grp) in computed_pairs):
                    continue

                # If no filter was informed, it will accept everything.
                if (self.inter_filter is not None
                        and not self.inter_filter.is_valid_pair(trgt_atm_grp,
                                                                nb_atm_grp)):
                    continue

                computed_pairs.add((trgt_atm_grp, nb_atm_grp))

                feat_pairs = list(product(trgt_atm_grp.features,
                                          nb_atm_grp.features))
                feat_pairs = filter(lambda x: self.is_feature_pair_valid(*x),
                                    feat_pairs)

                # If the groups belongs to the same molecule
                # (intramolecule interaction).
                is_intramol_inter = \
                    self._is_intramol_inter(trgt_atm_grp, nb_atm_grp)
                shortest_path_length = None

                for pair in feat_pairs:
                    # It will ignore interactions for atoms in the same
                    # molecule that are separated from each other by only
                    # N bonds. Covalent bonds keep atoms very tightly,
                    # producing distances lower than their sum of Van der Waals
                    # radius. As a consequence the algorithm will find a lot of
                    # false interactions.
                    #
                    # But, it will never skip pairs of Atom features because
                    # they are used to calculate covalent interactions.
                    if (pair[0].name != "Atom" and pair[1].name != "Atom"
                            and is_intramol_inter):
                        # Compute the shortest path only once. The reason not
                        # to precompute it outside the For is to avoid
                        # computing the algorithm for groups containing only
                        # Atom features.
                        if shortest_path_length is None:
                            # By providing a cutoff, it will force the
                            # algorithm to return paths only for groups
                            # connected by less than or equal N paths.
                            # So, if two groups return INF, it means
                            # they are a valid combination as they
                            # match the minimum bond separation.
                            cutoff = self.inter_config.get("min_bond_separation", 0)
                            shortest_path_length = trgt_atm_grp.get_shortest_path_length(nb_atm_grp, cutoff)

                        # If get_shortest_path_length() returns any value that
                        # is not infinite (INF), it means these two groups
                        # contain a path with at less than or equal to the
                        # cutoff 'min_bond_separation'.
                        # Therefore, ignore them.
                        if shortest_path_length != float('inf'):
                            continue

                    calc_inter_params = (trgt_atm_grp, nb_atm_grp) + pair
                    interactions = \
                        self._resolve_interactions(*calc_inter_params)
                    all_interactions.extend(interactions)

        if self.add_dependent_inter:
            dependent_interactions = \
                self.find_dependent_interactions(all_interactions)
            all_interactions.extend(dependent_interactions)

        # Get only unique interactions.
        all_interactions = set(all_interactions)

        # Remove potential inconsistences. For example: a hydrogen bond and
        # an unfavorable interation between the same atoms.
        self.remove_inconsistencies(all_interactions)

        if not self.add_h2o_pairs_with_no_target:
            self.remove_h2o_pairs_with_no_target(all_interactions)

        logger.debug("Number of potential interactions found: %d"
                     % len(all_interactions))

        return InteractionsManager(all_interactions)

    def _resolve_interactions(self, group1, group2, feat1, feat2):
        funcs = self.get_functions(feat1.name, feat2.name)
        if len(funcs) == 0:
            raise IllegalArgumentError("It does not exist a corresponding "
                                       "function to the features: '%s' and "
                                       "'%s'." % (feat1, feat2))

        interactions = []
        for func in funcs:
            result = func(self, (group1, group2, feat1, feat2)) or []
            interactions.extend(result)

        return interactions

    def find_dependent_interactions(self, interactions):
        """ Compute interactions that depend on other interactions.
        Currently, only water-bridged hydrogen bonds and salt bridges have a
        dependency on other interactions. The first, depends on two or more
        hydrogen bonds, while the second depends on an ionic and a hydrogen
        bond. The default value is False, which implies no dependent
        interaction will be computed.

        Parameters
        ----------
        interactions : :class:`~luna.interaction.type.InteractionType`
            Use these interactions to compute dependent interactions.
        """
        hbond_set = set()
        ionic_set = set()
        h2o_pairs = defaultdict(dict)
        dependent_interactions = set()

        # Save all hydrogen bonds involving waters and ionic interactions.
        for inter in interactions:
            if inter.type == "Hydrogen bond":
                comp1 = next(iter(inter.src_grp.compounds))
                comp2 = next(iter(inter.trgt_grp.compounds))

                if comp1.is_water():
                    h2o_key = inter.src_grp
                    secondKey = inter.trgt_grp

                    if secondKey not in h2o_pairs[h2o_key]:
                        h2o_pairs[h2o_key][secondKey] = inter

                if comp2.is_water():
                    h2o_key = inter.trgt_grp
                    secondKey = inter.src_grp

                    if secondKey not in h2o_pairs[h2o_key]:
                        h2o_pairs[h2o_key][secondKey] = inter

                hbond_set.add(inter)
            elif inter.type == "Ionic":
                ionic_set.add(inter)

        for h2o_key in h2o_pairs:
            pairs = combinations(h2o_pairs[h2o_key].keys(), 2)

            for src_grp, trgt_grp in pairs:

                # It ignores intramolecular interactions.
                if self._is_intramol_inter(src_grp, trgt_grp):
                    pass

                comp1 = next(iter(src_grp.compounds))
                comp2 = next(iter(trgt_grp.compounds))

                if isinstance(self.inter_filter, InteractionFilter):
                    if not self.inter_filter.is_valid_pair(src_grp, trgt_grp):
                        continue

                params = {"depends_on": [h2o_pairs[h2o_key][src_grp],
                                         h2o_pairs[h2o_key][trgt_grp]]}

                inter = InteractionType(src_grp, trgt_grp,
                                        "Water-bridged hydrogen bond",
                                        directional=True, params=params)
                dependent_interactions.add(inter)

        # It will try to match Hydrogen bonds and Ionic interactions
        # involving the same chemical groups to attribute salt bridges.
        sb_groups = set()
        for hbond, ionic in product(hbond_set, ionic_set):

            condA = (ionic.src_grp.has_atom(hbond.src_grp.atoms[0])
                     and ionic.trgt_grp.has_atom(hbond.trgt_grp.atoms[0]))

            condB = (ionic.src_grp.has_atom(hbond.trgt_grp.atoms[0])
                     and ionic.trgt_grp.has_atom(hbond.src_grp.atoms[0]))

            # If an acceptor atom belongs to a negative group, and the donor
            # to a positive group (and vice-versa), it means that the
            # interaction occurs between the same meioties. However, just one
            # condition should occur. For example, it is not possible that an
            # acceptor atom belongs to a negative and positive group at the
            # same time.
            if condA ^ condB:
                key1 = (ionic.src_grp, ionic.trgt_grp)
                key2 = (ionic.trgt_grp, ionic.src_grp)

                if key1 in sb_groups or key2 in sb_groups:
                    continue

                if isinstance(self.inter_filter, InteractionFilter):
                    if not self.inter_filter.is_valid_pair(ionic.src_grp,
                                                           ionic.trgt_grp):
                        continue

                sb_groups.add(key1)
                params = {"depends_on": [hbond, ionic]}

                inter = InteractionType(ionic.src_grp,
                                        ionic.trgt_grp,
                                        "Salt bridge",
                                        params=params)
                dependent_interactions.add(inter)

        return dependent_interactions

    def remove_inconsistencies(self, interactions):
        """Remove conflicts between interactions in ``interactions``.

        .. note::
            By default, LUNA defines conflicts as any unfavorable dipole
            interaction involving an atom establishing a hydrogen bond.
            Due to the strength of a hydrogen bond, atoms may approximate more
            to each other, which sometimes may cause unfavorable interactions
            involving dipoles to be detected. To avoid such conflicts, LUNA
            removes the unfavorable interactions.

            Also, it may occur that unfavorable dipole interactions are
            detected for amide-aromatic stackings, in which the aromatic
            ring contains heteroatoms. To avoid such conflicts, LUNA also
            removes those unfavorable interactions.

        Parameters
        ----------
        interactions : :class:`~luna.interaction.type.InteractionType`
        """
        amide_inconsistences = defaultdict(list)
        hbond_inconsistences = defaultdict(list)
        for inter in interactions:
            if inter.type == "Unfavorable nucleophile-nucleophile":
                # A nucleophile may have only 1 atom (water oxygen).
                atm1 = inter.src_grp.atoms[0]
                # If a nucleophile has 2 atoms, it will select the partially
                # negative atom based on the electronegativity.
                if len(inter.src_grp.atoms) == 2:
                    atm1 = (inter.src_grp.atoms[0]
                            if (inter.src_grp.atoms[0].electronegativity
                                > inter.src_grp.atoms[1].electronegativity)
                            else inter.src_grp.atoms[1])

                # A nucleophile may have only 1 atom (water oxygen).
                atm2 = inter.trgt_grp.atoms[0]
                # If a nucleophile has 2 atoms, it will select the partially
                # negative atom based on the electronegativity.
                if len(inter.trgt_grp.atoms) == 2:
                    atm2 = (inter.trgt_grp.atoms[0]
                            if (inter.trgt_grp.atoms[0].electronegativity
                                > inter.trgt_grp.atoms[1].electronegativity)
                            else inter.trgt_grp.atoms[1])
                key = (atm1, atm2)
                if (atm2, atm1) in hbond_inconsistences:
                    key = (atm2, atm1)
                hbond_inconsistences[key].append(inter)
            elif inter.type == "Unfavorable anion-nucleophile":
                nucl_grp = (inter.src_grp
                            if any([f.name == "Nucleophile"
                                    for f in inter.src_grp.features])
                            else inter.trgt_grp)
                # A nucleophile may have only 1 atom (water oxygen).
                nucl_atm = nucl_grp.atoms[0]
                # If a nucleophile has 2 atoms, it will select the partially
                # negative atom based on the electronegativity.
                if len(nucl_grp.atoms) == 2:
                    nucl_atm = (nucl_grp.atoms[0]
                                if (nucl_grp.atoms[0].electronegativity
                                    > nucl_grp.atoms[1].electronegativity)
                                else nucl_grp.atoms[1])
                anion_grp = inter.get_partner(nucl_grp)
                for anion_atm in anion_grp.atoms:
                    key = (nucl_atm, anion_atm)
                    if (anion_atm, nucl_atm) in hbond_inconsistences:
                        key = (anion_atm, nucl_atm)
                    hbond_inconsistences[key].append(inter)
            elif inter.type == "Hydrogen bond":
                # Acceptor and donor atoms always have only one atom.
                atm1, atm2 = (inter.src_grp.atoms[0], inter.trgt_grp.atoms[0])

                key = (atm1, atm2)
                if (atm2, atm1) in hbond_inconsistences:
                    key = (atm2, atm1)
                hbond_inconsistences[key].append(inter)
            elif inter.type == "Amide-aromatic stacking":
                amide_grp = (inter.src_grp
                             if any([f.name == "Amide"
                                     for f in inter.src_grp.features])
                             else inter.trgt_grp)
                arom_grp = inter.get_partner(amide_grp)
                for amide_atm in amide_grp.atoms:
                    amide_inconsistences[(amide_atm, arom_grp)].append(inter)
            elif inter.type == "Unfavorable cation-electrophile":
                elect_grp = (inter.src_grp
                             if any([f.name == "Nucleophile"
                                     for f in inter.src_grp.features])
                             else inter.trgt_grp)
                # A nucleophile may have only 1 atom (water oxygen).
                elect_atm = elect_grp.atoms[0]
                # If a nucleophile has 2 atoms, it will select the partially
                # negative atom based on the electronegativity.
                if len(elect_grp.atoms) == 2:
                    elect_atm = (elect_grp.atoms[0]
                                 if (elect_grp.atoms[0].electronegativity
                                     > elect_grp.atoms[1].electronegativity)
                                 else elect_grp.atoms[1])

                cation_grp = inter.get_partner(elect_grp)
                amide_inconsistences[(elect_atm, cation_grp)].append(inter)

        inconsistencies = set()
        for (atm1, atm2), inters in hbond_inconsistences.items():
            if (len(inters) > 1
                    and any([i.type == "Hydrogen bond" for i in inters])):
                inconsistencies.update([i for i in inters
                                        if i.type != "Hydrogen bond"])

        for (amide_atm, arom_grp), inters in amide_inconsistences.items():
            if len(inters) > 1:
                inter_name = "Amide-aromatic stacking"
                inconsistencies.update([i for i in inters
                                        if i.type != inter_name])

        interactions -= inconsistencies

        # Clear the references of each interaction from the AtomGroup objects.
        for inter in inconsistencies:
            inter.clear_refs()

    def remove_h2o_pairs_with_no_target(self, interactions):
        """Remove interactions of water with atoms and atom groups that do not
        belong to the target of LUNA's analysis, which are chains or molecules
        defined as an :class:`~luna.mol.entry.Entry` instance.

        Parameters
        ----------
        interactions : :class:`~luna.interaction.type.InteractionType`
        """
        valid_h2o_set = set()
        h2o_interactions = set()
        for inter in interactions:

            if inter.src_grp.is_water() or inter.trgt_grp.is_water():
                h2o_interactions.add(inter)

                # If a water is interacting with the ligand.
                if inter.src_grp.is_water() and inter.trgt_grp.has_target():
                    valid_h2o_set.add(inter.src_grp)

                if inter.trgt_grp.is_water() and inter.src_grp.has_target():
                    valid_h2o_set.add(inter.trgt_grp)

        inters_to_remove = set()
        for inter in h2o_interactions:
            if inter.src_grp.is_water() and inter.trgt_grp.is_water():
                if (inter.src_grp not in valid_h2o_set
                        or inter.trgt_grp not in valid_h2o_set):
                    inters_to_remove.add(inter)
            else:
                if (inter.src_grp.is_water()
                        and inter.src_grp not in valid_h2o_set):
                    inters_to_remove.add(inter)

                if (inter.trgt_grp.is_water()
                        and inter.trgt_grp not in valid_h2o_set):
                    inters_to_remove.add(inter)

        interactions -= inters_to_remove

        # Clear the references of each interaction from the AtomGroup objects.
        for inter in inters_to_remove:
            inter.clear_refs()

    def _default_functions(self):
        #         Hydrophobic interaction
        return {("Hydrophobic", "Hydrophobic"): [self.calc_hydrop],
                ("Hydrophobe", "Hydrophobe"): [self.calc_hydrop],

                # Hydrogen bond
                ("Donor", "Acceptor"): [self.calc_hbond],

                # Weak hydrogen bond
                ("WeakDonor", "Acceptor"): [self.calc_weak_hbond],
                ("WeakDonor", "WeakAcceptor"): [self.calc_weak_hbond],
                ("Donor", "Aromatic"): [self.calc_hbond_pi],
                ("WeakDonor", "Aromatic"): [self.calc_hbond_pi],

                # Halogen bond
                ("HalogenDonor", "Acceptor"): [self.calc_xbond],
                ("HalogenDonor", "Aromatic"): [self.calc_xbond_pi],

                # Chalcogen bond
                ("ChalcogenDonor", "Acceptor"): [self.calc_chalc_bond],
                ("ChalcogenDonor", "Aromatic"): [self.calc_chalc_bond_pi],

                # Stackings
                ("Aromatic", "Aromatic"): [self.calc_pi_pi],
                ("Amide", "Aromatic"): [self.calc_amide_pi],
                ("Positive", "Aromatic"): [self.calc_cation_pi],
                ("PosIonizable", "Aromatic"): [self.calc_cation_pi],
                ("PositivelyIonizable", "Aromatic"): [self.calc_cation_pi],

                # Ionic interaction
                ("NegativelyIonizable",
                 "PositivelyIonizable"): [self.calc_ionic],
                ("NegIonizable", "PosIonizable"): [self.calc_ionic],
                ("Negative", "Positive"): [self.calc_ionic],

                # Repulsive interaction
                ("NegativelyIonizable",
                 "NegativelyIonizable"): [self.calc_repulsive],
                ("PositivelyIonizable",
                 "PositivelyIonizable"): [self.calc_repulsive],
                ("NegIonizable", "NegIonizable"): [self.calc_repulsive],
                ("PosIonizable", "PosIonizable"): [self.calc_repulsive],
                ("Negative", "Negative"): [self.calc_repulsive],
                ("Positive", "Positive"): [self.calc_repulsive],

                # Favorable multipolar interactions.
                ("Nucleophile", "Electrophile"): [self.calc_multipolar],

                # Unfavorable multipolar interactions.
                ("Nucleophile", "Nucleophile"): [self.calc_multipolar],
                ("Electrophile", "Electrophile"): [self.calc_multipolar],

                # Favorable ion-dipole interactions
                ("Nucleophile",
                 "PositivelyIonizable"): [self.calc_ion_multipole],
                ("Nucleophile",
                 "PosIonizable"): [self.calc_ion_multipole],
                ("Nucleophile",
                 "Positive"): [self.calc_ion_multipole],
                ("Electrophile",
                 "NegativelyIonizable"): [self.calc_ion_multipole],
                ("Electrophile",
                 "NegIonizable"): [self.calc_ion_multipole],
                ("Electrophile",
                 "Negative"): [self.calc_ion_multipole],

                # Unfavorable ion-dipole interactions
                ("Nucleophile",
                 "NegativelyIonizable"): [self.calc_ion_multipole],
                ("Nucleophile",
                 "NegIonizable"): [self.calc_ion_multipole],
                ("Nucleophile",
                 "Negative"): [self.calc_ion_multipole],
                ("Electrophile",
                 "PositivelyIonizable"): [self.calc_ion_multipole],
                ("Electrophile",
                 "PosIonizable"): [self.calc_ion_multipole],
                ("Electrophile",
                 "Positive"): [self.calc_ion_multipole],

                ("Atom", "Metal"): [self.calc_metal_coord],

                # Proximal, covalent, vdw, clash
                ("Atom", "Atom"): [self.calc_atom_atom, self.calc_proximal]}

        # TODO: Include:

        # Anion - pi system

        # disulfide bond

        # Weak donor - weak acceptor

        # Agostic and Hydrogen-Bonding X–H· · · M
        # agnostic, anagostic

        # Metalic complex

        # REF: https://onlinelibrary.wiley.com/doi/epdf/10.1002/anie.200390319
        # aromatic between hbond arrays

    @staticmethod
    def calc_cation_pi(self, params):
        """Default method to calculate cation-pi interactions.

        Parameters
        ----------
        params : tuple of (:class:`~luna.mol.groups.AtomGroup`,\
                           :class:`~luna.mol.groups.AtomGroup`,\
                           :class:`~luna.mol.features.ChemicalFeature`,\
                           :class:`~luna.mol.features.ChemicalFeature`)
            The tuple follows the order (:math:`A`, :math:`B`, :math:`A_f`,
            :math:`B_f`), where :math:`A` and :math:`B` are two
            :class:`~luna.mol.groups.AtomGroup` objects, and :math:`A_f` and
            :math:`B_f` are their features
            (:class:`~luna.mol.features.ChemicalFeature` objects),
            respectively.

        Returns
        -------
         : list
        """
        if not self.add_non_cov:
            return []

        group1, group2, feat1, feat2 = params
        interactions = []

        cc_dist = im.euclidean_distance(group1.centroid, group2.centroid)

        if (self.is_within_boundary(cc_dist, "bsite_cutoff", le)
                and self.is_within_boundary(cc_dist,
                                            "max_dist_cation_pi_inter",
                                            le)):
            params = {"dist_cation_pi_inter": cc_dist}
            inter = InteractionType(group1, group2, "Cation-pi", params=params)

            interactions.append(inter)
        return interactions

    @staticmethod
    def calc_pi_pi(self, params):
        """Default method to calculate aromatic stackings.

        Parameters
        ----------
        params : tuple of (:class:`~luna.mol.groups.AtomGroup`,\
                           :class:`~luna.mol.groups.AtomGroup`,\
                           :class:`~luna.mol.features.ChemicalFeature`,\
                           :class:`~luna.mol.features.ChemicalFeature`)
            The tuple follows the order (:math:`A`, :math:`B`, :math:`A_f`,
            :math:`B_f`), where :math:`A` and :math:`B` are two
            :class:`~luna.mol.groups.AtomGroup` objects, and :math:`A_f` and
            :math:`B_f` are their features
            (:class:`~luna.mol.features.ChemicalFeature` objects),
            respectively.

        Returns
        -------
         : list
        """
        if not self.add_non_cov:
            return []

        ring1, ring2, feat1, feat2 = params
        interactions = []

        cc_dist = im.euclidean_distance(ring1.centroid, ring2.centroid)

        if (self.is_within_boundary(cc_dist, "bsite_cutoff", le)
                and self.is_within_boundary(cc_dist, "max_cc_dist_pi_pi_inter",
                                            le)):

            dihedral_angle = im.to_quad1(im.angle(ring1.normal, ring2.normal))

            min_disp_angle = float("Inf")
            for r1, r2 in [(ring1, ring2), (ring2, ring1)]:
                cc_vect = r2.centroid - r1.centroid
                disp_angle = im.to_quad1(im.angle(r1.normal, cc_vect))

                if disp_angle < min_disp_angle:
                    ring1, ring2 = r1, r2
                    min_disp_angle = disp_angle

            criteria = ["min_dihed_ang_slope_pi_pi_inter",
                        "max_dihed_ang_slope_pi_pi_inter",
                        "min_disp_ang_offset_pi_pi_inter",
                        "max_disp_ang_offset_pi_pi_inter"]

            # If the angle criterion were not defined, a specific Pi-stacking
            # definition is not possible as it depends on angle criterion.
            # Therefore, a more general classification is used instead, i.e.,
            # all interactions will be Pi-stacking.
            if any([c not in self.inter_config for c in criteria]):
                inter_type = "Pi-stacking"
            elif self.is_within_boundary(min_disp_angle,
                                         "min_disp_ang_offset_pi_pi_inter",
                                         le):
                if self.is_within_boundary(dihedral_angle,
                                           "min_dihed_ang_slope_pi_pi_inter",
                                           le):
                    inter_type = "Face-to-face pi-stacking"
                elif self.is_within_boundary(dihedral_angle,
                                             "max_dihed_ang_slope_pi_pi_inter",
                                             ge):
                    inter_type = "Face-to-edge pi-stacking"
                else:
                    inter_type = "Face-to-slope pi-stacking"
            elif self.is_within_boundary(min_disp_angle,
                                         "max_disp_ang_offset_pi_pi_inter",
                                         ge):
                if self.is_within_boundary(dihedral_angle,
                                           "min_dihed_ang_slope_pi_pi_inter",
                                           le):
                    inter_type = "Edge-to-edge pi-stacking"
                elif self.is_within_boundary(dihedral_angle,
                                             "max_dihed_ang_slope_pi_pi_inter",
                                             ge):
                    inter_type = "Edge-to-face pi-stacking"
                else:
                    inter_type = "Edge-to-slope pi-stacking"
            else:
                if self.is_within_boundary(dihedral_angle,
                                           "min_dihed_ang_slope_pi_pi_inter",
                                           le):
                    inter_type = "Displaced face-to-face pi-stacking"
                elif self.is_within_boundary(dihedral_angle,
                                             "max_dihed_ang_slope_pi_pi_inter",
                                             ge):
                    inter_type = "Displaced face-to-edge pi-stacking"
                else:
                    inter_type = "Displaced face-to-slope pi-stacking"

            params = {"cc_dist_pi_pi_inter": cc_dist,
                      "dihed_ang_pi_pi_inter": dihedral_angle,
                      "disp_ang_pi_pi_inter": min_disp_angle}

            inter = InteractionType(ring1,
                                    ring2,
                                    inter_type,
                                    directional=True,
                                    params=params)
            interactions.append(inter)

        return interactions

    @staticmethod
    def calc_amide_pi(self, params):
        """Default method to calculate amide-pi stackings.

        Parameters
        ----------
        params : tuple of (:class:`~luna.mol.groups.AtomGroup`,\
                           :class:`~luna.mol.groups.AtomGroup`,\
                           :class:`~luna.mol.features.ChemicalFeature`,\
                           :class:`~luna.mol.features.ChemicalFeature`)
            The tuple follows the order (:math:`A`, :math:`B`, :math:`A_f`,
            :math:`B_f`), where :math:`A` and :math:`B` are two
            :class:`~luna.mol.groups.AtomGroup` objects, and :math:`A_f` and
            :math:`B_f` are their features
            (:class:`~luna.mol.features.ChemicalFeature` objects),
            respectively.

        Returns
        -------
         : list
        """
        if not self.add_non_cov:
            return []

        group1, group2, feat1, feat2 = params
        interactions = []

        if feat1.name == "Aromatic" and feat2.name == "Amide":
            ring_grp = group1
            amide_grp = group2
        elif feat2.name == "Aromatic" and feat1.name == "Amide":
            ring_grp = group2
            amide_grp = group1
        else:
            logger.warning("Amide-aromatic interactions require an aromatic "
                           "and an amide group. However, the informed groups "
                           "have the features %s and %s."
                           % (group1.feature_names, group2.feature_names))
            return []

        # Distance between the amide and ring centroids.
        cc_dist = im.euclidean_distance(ring_grp.centroid, amide_grp.centroid)

        if (self.is_within_boundary(cc_dist, "bsite_cutoff", le)
                and self.is_within_boundary(cc_dist,
                                            "max_cc_dist_amide_pi_inter", le)):

            dihedral_angle = im.to_quad1(im.angle(ring_grp.normal,
                                                  amide_grp.normal))
            cc_vect = amide_grp.centroid - ring_grp.centroid
            disp_angle = im.to_quad1(im.angle(ring_grp.normal, cc_vect))

            if (self.is_within_boundary(dihedral_angle,
                                        "max_dihed_ang_amide_pi_inter", le)
                    and self.is_within_boundary(disp_angle,
                                                "max_disp_ang_pi_pi_inter",
                                                le)):

                params = {"cc_dist_amide_pi_inter": cc_dist,
                          "dihed_ang_amide_pi_inter": dihedral_angle,
                          "disp_ang_amide_pi_inter": disp_angle}

                inter = InteractionType(group1,
                                        group2,
                                        "Amide-aromatic stacking",
                                        params=params)
                interactions.append(inter)
        return interactions

    @staticmethod
    def calc_hydrop(self, params):
        """Default method to calculate hydrophobic interactons.

        Parameters
        ----------
        params : tuple of (:class:`~luna.mol.groups.AtomGroup`,\
                           :class:`~luna.mol.groups.AtomGroup`,\
                           :class:`~luna.mol.features.ChemicalFeature`,\
                           :class:`~luna.mol.features.ChemicalFeature`)
            The tuple follows the order (:math:`A`, :math:`B`, :math:`A_f`,
            :math:`B_f`), where :math:`A` and :math:`B` are two
            :class:`~luna.mol.groups.AtomGroup` objects, and :math:`A_f` and
            :math:`B_f` are their features
            (:class:`~luna.mol.features.ChemicalFeature` objects),
            respectively.

        Returns
        -------
         : list
        """
        if not self.add_non_cov:
            return []

        group1, group2, feat1, feat2 = params
        interactions = []

        if ((feat1.name != "Hydrophobic" and feat1.name != "Hydrophobe")
                or (feat2.name != "Hydrophobic"
                    and feat2.name != "Hydrophobe")):
            logger.warning("Hydrophobic interactions require hydrophobic "
                           "atoms or hydrophobes (group of hydrophobic "
                           "atoms). However, the informed groups have the "
                           "features %s and %s."
                           % (group1.feature_names, group2.feature_names))
            return []

        # Check if the interaction involves the same compound.
        # For these cases, we ignore hydrophobic interactions.
        if self._is_intramol_inter(group1, group2):
            return []

        # Verify if the groups contain the required number of atoms to
        # form a valid surface.
        if (not self.is_within_boundary(len(group1.atoms), "min_surf_size", ge)
                or not self.is_within_boundary(len(group2.atoms),
                                               "min_surf_size", ge)):
            return []

        interacting_atms_in_surf1 = set()
        interacting_atms_in_surf2 = set()
        min_cc_dist = float('Inf')
        for atm1, atm2 in product(group1.atoms, group2.atoms):
            cc_dist = atm1 - atm2

            if self.is_within_boundary(cc_dist, "max_dist_hydrop_inter", le):
                interacting_atms_in_surf1.add(atm1)
                interacting_atms_in_surf2.add(atm2)

                if cc_dist < min_cc_dist:
                    min_cc_dist = cc_dist

        # Verify if the number of interacting atoms attends the required number
        # of interating atoms per surface.
        if (not self.is_within_boundary(len(interacting_atms_in_surf1), "min_inter_atom_in_surf", ge)
                or not self.is_within_boundary(len(interacting_atms_in_surf2), "min_inter_atom_in_surf", ge)):
            return []

        if (self.is_within_boundary(min_cc_dist, "bsite_cutoff", le)
                and self.is_within_boundary(min_cc_dist,
                                            "max_dist_hydrop_inter", le)):

            params = {"dist_hydrop_inter": min_cc_dist}
            inter = InteractionType(group1,
                                    group2,
                                    "Hydrophobic",
                                    params=params)
            interactions.append(inter)

        return interactions

    @staticmethod
    def calc_ion_multipole(self, params):
        """Default method to calculate favorable and unfavorable ion-dipole
        interactions.

        Parameters
        ----------
        params : tuple of (:class:`~luna.mol.groups.AtomGroup`,\
                           :class:`~luna.mol.groups.AtomGroup`,\
                           :class:`~luna.mol.features.ChemicalFeature`,\
                           :class:`~luna.mol.features.ChemicalFeature`)
            The tuple follows the order (:math:`A`, :math:`B`, :math:`A_f`,
            :math:`B_f`), where :math:`A` and :math:`B` are two
            :class:`~luna.mol.groups.AtomGroup` objects, and :math:`A_f` and
            :math:`B_f` are their features
            (:class:`~luna.mol.features.ChemicalFeature` objects),
            respectively.

        Returns
        -------
         : list
        """
        if not self.add_non_cov:
            return []

        group1, group2, feat1, feat2 = params
        interactions = []

        # Favorable cation-nucleophile interaction
        if feat1.name == "Nucleophile" and feat2.name in CATIONS:
            dipole_grp, dipole_type = group1, "Nucleophile"
            ion_grp, ion_type = group2, "Cation"
        elif feat1.name in CATIONS and feat2.name == "Nucleophile":
            dipole_grp, dipole_type = group2, "Nucleophile"
            ion_grp, ion_type = group1, "Cation"
        # Favorable anion-electrophile interaction
        elif feat1.name == "Electrophile" and feat2.name in ANIONS:
            dipole_grp, dipole_type = group1, "Electrophile"
            ion_grp, ion_type = group2, "Anion"
        elif feat1.name in ANIONS and feat2.name == "Electrophile":
            dipole_grp, dipole_type = group2, "Electrophile"
            ion_grp, ion_type = group1, "Anion"
        # Unfavorable anion-nucleophile interaction
        elif feat1.name == "Nucleophile" and feat2.name in ANIONS:
            dipole_grp, dipole_type = group1, "Nucleophile"
            ion_grp, ion_type = group2, "Anion"
        elif feat1.name in ANIONS and feat2.name == "Nucleophile":
            dipole_grp, dipole_type = group2, "Nucleophile"
            ion_grp, ion_type = group1, "Anion"
        # Unfavorable cation-electrophile interaction
        elif feat1.name == "Electrophile" and feat2.name in CATIONS:
            dipole_grp, dipole_type = group1, "Electrophile"
            ion_grp, ion_type = group2, "Cation"
        elif feat1.name in CATIONS and feat2.name == "Electrophile":
            dipole_grp, dipole_type = group2, "Electrophile"
            ion_grp, ion_type = group1, "Cation"
        else:
            logger.warning("Ion-dipole interactions require a dipole and an "
                           "ion group. However, the informed groups have the "
                           "features %s and %s."
                           % (group1.feature_names, group2.feature_names))
            return []

        # A nucleophile may have only 1 atom (water oxygen).
        part_charged_atm = dipole_grp.atoms[0]
        # If a nucleophile has 2 atoms, it will select the partially negative
        # atom based on the electronegativity.
        if len(dipole_grp.atoms) == 2 and dipole_type == "Nucleophile":
            part_charged_atm = \
                (dipole_grp.atoms[0]
                 if (dipole_grp.atoms[0].electronegativity
                     > dipole_grp.atoms[1].electronegativity)
                 else dipole_grp.atoms[1])

        # If an electrophile has two atoms. It will select the partially
        # negative atom based on the electronegativity.
        elif len(dipole_grp.atoms) == 2 and dipole_type == "Electrophile":
            part_charged_atm = \
                (dipole_grp.atoms[0]
                 if (dipole_grp.atoms[0].electronegativity
                     < dipole_grp.atoms[1].electronegativity)
                 else dipole_grp.atoms[1])

        # Distance between the ion and the dipole.
        id_dist = im.euclidean_distance(part_charged_atm.coord,
                                        ion_grp.centroid)

        if (self.is_within_boundary(id_dist, "bsite_cutoff", le)
                and self.is_within_boundary(id_dist,
                                            "max_id_dist_ion_multipole_inter",
                                            le)):

            idy_angle = -1
            if len(dipole_grp.atoms) == 2:
                # Model: I ... D-Y, where I is the ion, D the dipole atom of
                # interest (the electrophile or nucleophile),
                # and Y is its counterpart.
                y_atm = dipole_grp.atoms[1] if dipole_grp.atoms[0] == part_charged_atm else dipole_grp.atoms[0]
                di_vect = ion_grp.centroid - part_charged_atm.coord
                dy_vect = y_atm.coord - part_charged_atm.coord
                idy_angle = im.angle(di_vect, dy_vect)

            # Dipoles containing only one atom are allowed to pass without
            # checking the angle IDY.
            if len(dipole_grp.atoms) == 1 or self.is_within_boundary(idy_angle, "min_idy_ang_ion_multipole_inter", ge):

                dipole_nb_coords = [nbi.coord
                                    for nbi in part_charged_atm.neighbors_info
                                    if nbi.atomic_num != 1]
                params = {}
                if len(dipole_nb_coords) > 1:
                    dipole_normal \
                        = im.calc_normal(dipole_nb_coords
                                         + [part_charged_atm.coord])
                    disp_angle = im.to_quad1(im.angle(dipole_normal, di_vect))

                    prop_name = "max_disp_ang_ion_multipole_inter"
                    if self.is_within_boundary(disp_angle, prop_name, le):
                        params = {"id_dist_ion_multipole_inter": id_dist,
                                  "idy_ang_ion_multipole_inter": idy_angle,
                                  "disp_ang_ion_multipole_inter": disp_angle}
                else:
                    params = {"id_dist_ion_multipole_inter": id_dist,
                              "idy_ang_ion_multipole_inter": idy_angle,
                              "disp_ang_ion_multipole_inter": -1}

                if params:
                    if dipole_type == "Nucleophile" and ion_type == "Cation":
                        inter = InteractionType(dipole_grp, ion_grp,
                                                "Cation-nucleophile",
                                                params=params)
                        interactions.append(inter)
                    elif dipole_type == "Nucleophile" and ion_type == "Anion":
                        inter_name = "Unfavorable anion-nucleophile"
                        inter = InteractionType(dipole_grp, ion_grp,
                                                inter_name, params=params)
                        interactions.append(inter)
                    elif dipole_type == "Electrophile" and ion_type == "Anion":
                        inter = InteractionType(ion_grp, dipole_grp,
                                                "Anion-electrophile",
                                                params=params)
                        interactions.append(inter)
                    elif (dipole_type == "Electrophile"
                            and ion_type == "Cation"):
                        inter_name = "Unfavorable cation-electrophile"
                        inter = InteractionType(ion_grp, dipole_grp,
                                                inter_name, params=params)
                        interactions.append(inter)
        return interactions

    @staticmethod
    def calc_multipolar(self, params):
        """Default method to calculate favorable and unfavorable dipole-dipole
        interactions.

        Parameters
        ----------
        params : tuple of (:class:`~luna.mol.groups.AtomGroup`,\
                           :class:`~luna.mol.groups.AtomGroup`,\
                           :class:`~luna.mol.features.ChemicalFeature`,\
                           :class:`~luna.mol.features.ChemicalFeature`)
            The tuple follows the order (:math:`A`, :math:`B`, :math:`A_f`,
            :math:`B_f`), where :math:`A` and :math:`B` are two
            :class:`~luna.mol.groups.AtomGroup` objects, and :math:`A_f` and
            :math:`B_f` are their features
            (:class:`~luna.mol.features.ChemicalFeature` objects),
            respectively.

        Returns
        -------
         : list
        """
        if not self.add_non_cov:
            return []

        group1, group2, feat1, feat2 = params
        interactions = []

        if (len(group1.atoms) != 1 and len(group1.atoms) != 2
                and len(group2.atoms) != 1 and len(group2.atoms) != 2):
            logger.warning("A dipole group should have 1 (for cases when the "
                           "atom has only hydrogens bonded to it) or 2 atoms. "
                           "However, the informed groups '%s' and '%s' have "
                           "%d and %d atoms, respectively."
                           % (group1, group2,
                              len(group1.atoms), len(group2.atoms)))
            return []

        # The reference dipole will always be the second one, i.e., one of its
        # atom will be the center in the angle NEY.
        #
        # Favorable interactions: in these cases, the Dipole 1 will always be
        # the nucleophile and the Dipole 2 the electrophile in order to
        # represent the nucleophile atack, i.e., the angles calculated using
        # the dipole 2 as reference represents how the nucleophile aproximate
        # the electrophile.
        if feat1.name == "Nucleophile" and feat2.name == "Electrophile":
            dipole_grp1, dipole_type1 = group1, feat1.name
            dipole_grp2, dipole_type2 = group2, feat2.name
        elif feat2.name == "Nucleophile" and feat1.name == "Electrophile":
            dipole_grp1, dipole_type1 = group2, feat2.name
            dipole_grp2, dipole_type2 = group1, feat1.name
        # Unfavorable interactions: in these cases, the reference dipole will
        # depend on the number of atoms in the dipoles. With dipoles containing
        # 1 atom, it takes a generous approach by ignoring angles and accepting
        # everything.
        #
        # With dipoles containing two atoms, it requires that at least one of
        # the angles fits the rules to be accepted.
        elif (feat1.name == feat2.name
                and (feat1.name == "Nucleophile"
                     or feat1.name == "Electrophile")):
            # If only one group contains 1 atom, use it as the dipole 2 because
            # it is used as the reference to calculate the NEY angle. Since we
            # take a generous approach, with one atom no angle will be
            # calculated and the interaction will be accepted.
            if len(group1.atoms) == 1 and len(group2.atoms) == 2:
                dipole_grp1, dipole_type1 = group2, feat2.name
                dipole_grp2, dipole_type2 = group1, feat1.name
            # All the other number combinations ([2,1], [1,1], [2, 2]) come
            # here.
            else:
                dipole_grp1, dipole_type1 = group1, feat1.name
                dipole_grp2, dipole_type2 = group2, feat2.name
        else:
            logger.warning("Multipolar interactions require a nucleophile and "
                           "an electrophile group. However, the informed "
                           "groups have the features %s and %s."
                           % (group1.feature_names, group2.feature_names))
            return []

        # Ignore dipoles containing at least one common atom, which can happen
        # to covalently bound dipoles. An example of it is the C-S-C
        # substructure that contains two dipoles.
        if any(atm in group2.atoms for atm in group1.atoms):
            return []

        # Atom 1 => Dipole 1
        #
        # A nucleophile may have only 1 atom (water oxygen).
        dipole_atm1 = dipole_grp1.atoms[0]
        # If it has 2 atoms, it will select the nucleophilic atom based on the
        # electronegativity.
        if len(dipole_grp1.atoms) == 2 and dipole_type1 == "Nucleophile":
            dipole_atm1 = (dipole_grp1.atoms[0]
                           if (dipole_grp1.atoms[0].electronegativity
                               > dipole_grp1.atoms[1].electronegativity)
                           else dipole_grp1.atoms[1])
        # Or, it will select the nucleophilic atom based on the
        # electronegativity.
        elif len(dipole_grp1.atoms) == 2 and dipole_type1 == "Electrophile":
            dipole_atm1 = (dipole_grp1.atoms[0]
                           if (dipole_grp1.atoms[0].electronegativity
                               < dipole_grp1.atoms[1].electronegativity)
                           else dipole_grp1.atoms[1])

        # Atom 2 => Dipole 2
        #
        # An electrophile may have only 1 atom. E.g.: NH4, although by default
        # we consider it as an ion.
        dipole_atm2 = dipole_grp2.atoms[0]
        # If it has 2 atoms, it will select the nucleophilic atom based on the
        # electronegativity.
        if len(dipole_grp2.atoms) == 2 and dipole_type2 == "Nucleophile":
            dipole_atm2 = (dipole_grp2.atoms[0]
                           if (dipole_grp2.atoms[0].electronegativity
                               > dipole_grp2.atoms[1].electronegativity)
                           else dipole_grp2.atoms[1])
        # Or, it will select the nucleophilic atom based on the
        # electronegativity.
        elif len(dipole_grp2.atoms) == 2 and dipole_type2 == "Electrophile":
            dipole_atm2 = (dipole_grp2.atoms[0]
                           if (dipole_grp2.atoms[0].electronegativity
                               < dipole_grp2.atoms[1].electronegativity)
                           else dipole_grp2.atoms[1])

        # Model for favorable interactions: A-N ... E-Y
        # Model for unfavorable interactions: A-N ... N-A, Y-E ... E-Y.
        #
        # Although there are two different models for unfavorable interactions,
        # the method for them are equal to the favorable interaction.
        # So, from now on, we will deal with them as if it was the first model.
        #
        # Distance between the nucleophile and electrophile.
        ne_dist = im.euclidean_distance(dipole_atm1.coord, dipole_atm2.coord)

        if (self.is_within_boundary(ne_dist, "bsite_cutoff", le)
                and self.is_within_boundary(ne_dist,
                                            "max_ne_dist_multipolar_inter",
                                            le)):

            # No angle can be calculated if the electrophile (dipole 2) has
            # only one atom.
            if len(dipole_grp2.atoms) == 1:
                params = {"ne_dist_multipolar_inter": ne_dist,
                          "ney_ang_multipolar_inter": -1,
                          "disp_ang_multipolar_inter": -1,
                          "an_ey_ang_multipolar_inter": -1}

                inter_type = ("Multipolar" if not dipole_type1 == dipole_type2
                              else "Unfavorable %s-%s"
                              % (dipole_type1.lower(), dipole_type2.lower()))
                inter = InteractionType(dipole_grp1, dipole_grp2, inter_type,
                                        directional=True, params=params)
                interactions.append(inter)

            else:
                dipole1 = (dipole_grp1, dipole_atm1, dipole_type1)
                dipole2 = (dipole_grp2, dipole_atm2, dipole_type2)

                combinations = [(dipole1, dipole2)]
                # For unfavorable interactions, it is necessary to evaluate
                # each combination of dipoles. So, it can produce two
                # interactions.
                if feat1.name == feat2.name:
                    combinations = [(dipole1, dipole2), (dipole2, dipole1)]

                for d1, d2 in combinations:
                    dipole_grp1, dipole_atm1, dipole_type1 = d1
                    dipole_grp2, dipole_atm2, dipole_type2 = d2

                    # Model: A-N ... E-Y
                    y_atm = (dipole_grp2.atoms[1]
                             if dipole_grp2.atoms[0] == dipole_atm2
                             else dipole_grp2.atoms[0])
                    en_vect = dipole_atm1.coord - dipole_atm2.coord
                    ey_vect = y_atm.coord - dipole_atm2.coord
                    ney_angle = im.angle(en_vect, ey_vect)

                    prop_name1 = "min_ney_ang_multipolar_inter"
                    prop_name2 = "max_ney_ang_multipolar_inter"
                    if (self.is_within_boundary(ney_angle, prop_name1, ge)
                            and self.is_within_boundary(ney_angle,
                                                        prop_name2, le)):

                        elect_nb_coords = \
                            [nbi.coord
                             for nbi in dipole_atm2.neighbors_info
                             if nbi.atomic_num != 1]
                        elect_normal = \
                            im.calc_normal(elect_nb_coords
                                           + [dipole_atm2.coord])
                        disp_angle = im.to_quad1(im.angle(elect_normal,
                                                          en_vect))

                        prop_name = "max_disp_ang_multipolar_inter"
                        if self.is_within_boundary(disp_angle, prop_name, le):

                            # If the nucleophile has two atoms, then we will be
                            # able to calculate the angle between the vectors AN
                            # and EY. This angle is necessary to define the
                            # orientation of the dipole.
                            if len(dipole_grp1.atoms) == 2:
                                # Model: A-N ... E-Y
                                a_atm = (dipole_grp1.atoms[1]
                                         if dipole_grp1.atoms[0] == dipole_atm1
                                         else dipole_grp1.atoms[0])
                                an_vect = dipole_atm1.coord - a_atm.coord
                                # Angle between vectors AN and EY
                                an_ey_vect_angle = im.angle(an_vect, ey_vect)

                                params = \
                                    {"ne_dist_multipolar_inter": ne_dist,
                                     "ney_ang_multipolar_inter": ney_angle,
                                     "disp_ang_multipolar_inter": disp_angle,
                                     "an_ey_ang_multipolar_inter": an_ey_vect_angle}

                                if not dipole_type1 == dipole_type2:
                                    if self.is_within_boundary(an_ey_vect_angle,
                                                               "max_an_ey_ang_para_multipolar_inter",
                                                               le):
                                        inter_name = "Parallel multipolar"
                                        inter = \
                                            InteractionType(dipole_grp1,
                                                            dipole_grp2,
                                                            inter_name,
                                                            directional=True,
                                                            params=params)
                                        interactions.append(inter)

                                    elif self.is_within_boundary(an_ey_vect_angle,
                                                                 "min_an_ey_ang_antipara_multipolar_inter",
                                                                 ge):
                                        inter_name = "Antiparallel multipolar"
                                        inter = \
                                            InteractionType(dipole_grp1,
                                                            dipole_grp2,
                                                            inter_name,
                                                            directional=True,
                                                            params=params)
                                        interactions.append(inter)

                                    elif (self.is_within_boundary(an_ey_vect_angle,
                                                                  "min_an_ey_ang_ortho_multipolar_inter",
                                                                  ge)
                                            and self.is_within_boundary(an_ey_vect_angle,
                                                                        "max_an_ey_ang_ortho_multipolar_inter",
                                                                        le)):
                                        inter_name = "Orthogonal multipolar"
                                        inter = \
                                            InteractionType(dipole_grp1,
                                                            dipole_grp2,
                                                            inter_name,
                                                            directional=True,
                                                            params=params)
                                        interactions.append(inter)

                                    else:
                                        inter_name = "Tilted multipolar"
                                        inter = \
                                            InteractionType(dipole_grp1,
                                                            dipole_grp2,
                                                            inter_name,
                                                            directional=True,
                                                            params=params)
                                        interactions.append(inter)

                                else:
                                    inter_type = ("Unfavorable %s-%s"
                                                  % (dipole_type1.lower(),
                                                     dipole_type2.lower()))
                                    inter = InteractionType(dipole_grp1,
                                                            dipole_grp2,
                                                            inter_type,
                                                            directional=True,
                                                            params=params)
                                    interactions.append(inter)

                            # Otherwise, ignore the angle AN and EY and add a
                            # general interaction (Multipolar) without a
                            # specific definition of the orientation. It will
                            # happen only with Water molecules.
                            else:
                                params = \
                                    {"ne_dist_multipolar_inter": ne_dist,
                                     "ney_ang_multipolar_inter": ney_angle,
                                     "disp_ang_multipolar_inter": disp_angle,
                                     "an_ey_ang_multipolar_inter": -1}
                                inter_type = \
                                    ("Multipolar"
                                     if (not dipole_type1 == dipole_type2)
                                     else "Unfavorable %s-%s"
                                     % (dipole_type1.lower(),
                                        dipole_type2.lower()))
                                inter = InteractionType(dipole_grp1,
                                                        dipole_grp2,
                                                        inter_type,
                                                        directional=True,
                                                        params=params)
                                interactions.append(inter)

        return interactions

    @staticmethod
    def calc_xbond_pi(self, params):
        """Default method to calculate halogen bonds between halogens and
        aromatic rings.

        Parameters
        ----------
        params : tuple of (:class:`~luna.mol.groups.AtomGroup`,\
                           :class:`~luna.mol.groups.AtomGroup`,\
                           :class:`~luna.mol.features.ChemicalFeature`,\
                           :class:`~luna.mol.features.ChemicalFeature`)
            The tuple follows the order (:math:`A`, :math:`B`, :math:`A_f`,
            :math:`B_f`), where :math:`A` and :math:`B` are two
            :class:`~luna.mol.groups.AtomGroup` objects, and :math:`A_f` and
            :math:`B_f` are their features
            (:class:`~luna.mol.features.ChemicalFeature` objects),
            respectively.

        Returns
        -------
         : list
        """
        if not self.add_non_cov:
            return []

        group1, group2, feat1, feat2 = params
        interactions = []

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
        xa_dist = im.euclidean_distance(donor_grp.centroid, ring_grp.centroid)

        if (self.is_within_boundary(xa_dist, "bsite_cutoff", le)
                and self.is_within_boundary(xa_dist,
                                            "max_xc_dist_xbond_inter", le)):

            ax_vect = donor_grp.centroid - ring_grp.centroid
            disp_angle = im.to_quad1(im.angle(ring_grp.normal, ax_vect))

            if (self.is_within_boundary(disp_angle,
                                        "max_disp_ang_xbond_inter", le)):
                # Interaction model: C-X ---- A
                # XA vector is always the same
                xa_vect = ring_grp.centroid - donor_grp.centroid

                # Defining angle CXA, in which A is the ring center
                # It may happen that X is covalently bound to more than one
                # group. In such cases the halogen may also form more than
                # one halogen bond.
                #
                # Ref: Cavallo, G. et al. The Halogen Bond. (2016).
                #
                carbon_coords = [nbi.coord
                                 for nbi in donor_atm.neighbors_info
                                 if nbi.atomic_num == 6]
                for c_coord in carbon_coords:
                    xc_vect = c_coord - donor_grp.centroid
                    cxa_angle = im.angle(xc_vect, xa_vect)

                    if (self.is_within_boundary(cxa_angle,
                                                "min_cxa_ang_xbond_inter",
                                                ge)):
                        params = {"xc_dist_xbond_inter": xa_dist,
                                  "disp_ang_xbond_inter": disp_angle,
                                  "cxa_ang_xbond_inter": cxa_angle}

                        inter = InteractionType(donor_grp,
                                                ring_grp,
                                                "Halogen-pi",
                                                directional=True,
                                                params=params)
                        interactions.append(inter)
        return interactions

    @staticmethod
    def calc_xbond(self, params):
        """Default method to calculate halogen bonds.

        Parameters
        ----------
        params : tuple of (:class:`~luna.mol.groups.AtomGroup`,\
                           :class:`~luna.mol.groups.AtomGroup`,\
                           :class:`~luna.mol.features.ChemicalFeature`,\
                           :class:`~luna.mol.features.ChemicalFeature`)
            The tuple follows the order (:math:`A`, :math:`B`, :math:`A_f`,
            :math:`B_f`), where :math:`A` and :math:`B` are two
            :class:`~luna.mol.groups.AtomGroup` objects, and :math:`A_f` and
            :math:`B_f` are their features
            (:class:`~luna.mol.features.ChemicalFeature` objects),
            respectively.

        Returns
        -------
         : list
        """
        if not self.add_non_cov:
            return []

        group1, group2, feat1, feat2 = params
        interactions = []

        if len(group1.atoms) != 1 or len(group2.atoms) != 1:
            logger.warning("One or more invalid atom groups were informed: "
                           "%s and %s. In halogen bonds, halogen donor "
                           "and acceptor groups should always contain only "
                           "one atom." % (group1, group2))
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
        xa_dist = im.euclidean_distance(donor_grp.centroid,
                                        acceptor_grp.centroid)

        if (self.is_within_boundary(xa_dist, "bsite_cutoff", le)
                and self.is_within_boundary(xa_dist,
                                            "max_xa_dist_xbond_inter", le)):

            # Interaction model: C-X ---- A-R
            # XA vector is always the same
            xa_vect = acceptor_grp.centroid - donor_grp.centroid

            # Defining the angle CXA
            # It may happen that X is covalently bound to more than one group.
            # In such cases the halogen may also form more than one
            # halogen bond. Ref: Cavallo, G. et al. The Halogen Bond. (2016).
            carbon_coords = [nbi.coord
                             for nbi in donor_atm.neighbors_info
                             if nbi.atomic_num == 6]

            # Interaction model: C-X ---- A-R.
            # R coordinates, in which R is a heavy atom.
            r_coords = [nbi.coord
                        for nbi in acceptor_atm.neighbors_info
                        if nbi.atomic_num != 1]

            for c_coord in carbon_coords:
                xc_vect = c_coord - donor_grp.centroid
                cxa_angle = im.angle(xc_vect, xa_vect)

                if self.is_within_boundary(cxa_angle,
                                           "min_cxa_ang_xbond_inter", ge):

                    # If no heavy atom is bonded to the acceptor, it means that
                    # only hydrogens may be bound to it. Then, we do not
                    # calculate the angles because hydrogens are too dynamic,
                    # i.e., the acceptor could be ionized or not at a specific
                    # moment in time and its hydrogens may be positioned in
                    # different ways.
                    if len(r_coords) == 0:
                        params = {"xa_dist_xbond_inter": xa_dist,
                                  "cxa_ang_xbond_inter": cxa_angle,
                                  "xar_ang_xbond_inter": -1}

                        inter = InteractionType(donor_grp,
                                                acceptor_grp,
                                                "Halogen bond",
                                                directional=True,
                                                params=params)
                        interactions.append(inter)
                    else:
                        # AX vector is always the same.
                        # Obs: the donor_grp is the Halogen (X).
                        ax_vect = donor_grp.centroid - acceptor_grp.centroid

                        lowest_xar_angle = None
                        for r_coord in r_coords:
                            ar_vect = r_coord - acceptor_grp.centroid
                            xar_angle = im.angle(ax_vect, ar_vect)

                            # Update the XAR angle with the lowest value.
                            if (lowest_xar_angle is None
                                    or xar_angle < lowest_xar_angle):
                                lowest_xar_angle = xar_angle

                        # The angle will be None when any R (heavy atom) atom
                        # was found. In this case, the criteria must always
                        # fail.
                        if lowest_xar_angle is None:
                            lowest_xar_angle = -1

                        prop_name = "min_xar_ang_xbond_inter"
                        if self.is_within_boundary(lowest_xar_angle,
                                                   prop_name, ge):
                            params = {"xa_dist_xbond_inter": xa_dist,
                                      "cxa_ang_xbond_inter": cxa_angle,
                                      "xar_ang_xbond_inter": lowest_xar_angle}

                            inter = InteractionType(donor_grp,
                                                    acceptor_grp,
                                                    "Halogen bond",
                                                    directional=True,
                                                    params=params)
                            interactions.append(inter)

        return interactions

    @staticmethod
    def calc_chalc_bond(self, params):
        """Default method to calculate chalcogen bonds.

        Parameters
        ----------
        params : tuple of (:class:`~luna.mol.groups.AtomGroup`,\
                           :class:`~luna.mol.groups.AtomGroup`,\
                           :class:`~luna.mol.features.ChemicalFeature`,\
                           :class:`~luna.mol.features.ChemicalFeature`)
            The tuple follows the order (:math:`A`, :math:`B`, :math:`A_f`,
            :math:`B_f`), where :math:`A` and :math:`B` are two
            :class:`~luna.mol.groups.AtomGroup` objects, and :math:`A_f`
            and :math:`B_f` are their features
            (:class:`~luna.mol.features.ChemicalFeature` objects),
            respectively.

        Returns
        -------
         : list
        """
        if not self.add_non_cov:
            return []

        group1, group2, feat1, feat2 = params
        interactions = []

        if len(group1.atoms) != 1 or len(group2.atoms) != 1:
            logger.warning("One or more invalid atom groups were informed: "
                           " %s and %s. In chalcogen bonds, chalcogen donor "
                           "and acceptor groups should always contain "
                           "only one atom." % (group1, group2))
            return []

        if (feat1.name == "Acceptor" and feat2.name == "ChalcogenDonor"):
            donor_grp = group2
            acceptor_grp = group1
        else:
            donor_grp = group1
            acceptor_grp = group2

        # There are always just one donor/acceptor atom.
        donor_atm = donor_grp.atoms[0]
        acceptor_atm = acceptor_grp.atoms[0]

        # Interaction model: R-Y --- A-N
        # Distance YA.
        ya_dist = im.euclidean_distance(donor_grp.centroid,
                                        acceptor_grp.centroid)

        prop_name = "max_ya_dist_ybond_inter"
        if (self.is_within_boundary(ya_dist, "bsite_cutoff", le)
                and self.is_within_boundary(ya_dist, prop_name, le)):

            # Interaction model: R-Y --- A-N
            # YA vector is always the same
            ya_vect = acceptor_grp.centroid - donor_grp.centroid

            # Defining the angle RYA
            r_atms = [nbi for nbi in donor_atm.neighbors_info
                      if nbi.atomic_num != 1]
            r_elems = sorted([nbi.atomic_num
                              for nbi in donor_atm.neighbors_info
                              if nbi.atomic_num != 1])

            # Isothiazoles have only one sigma-hole located on the oposite
            # site of the N. Therefore, we should only evaluate this
            # sigma-hole by ignoring the angle formed with the carbon.
            # Ref: https://doi.org/10.1021/jm501853m.
            ignore_carbon = False
            if len(r_elems) == 2 and r_elems[0] == 6 and r_elems[1] == 7:
                ignore_carbon = True

            # Interaction model: R-Y --- A-N.
            # N coordinates, in which N is a heavy atom.
            n_coords = [nbi.coord
                        for nbi in acceptor_atm.neighbors_info
                        if nbi.atomic_num != 1]

            for r_atm in r_atms:
                # Isothiazoles have only one sigma-hole located on the oposite
                # site of the N. So, we must ignore the Carbon.
                # Ref: https://doi.org/10.1021/jm501853m.
                if r_atm.atomic_num == 6 and ignore_carbon is True:
                    continue

                yr_vect = r_atm.coord - donor_grp.centroid
                rya_angle = im.angle(yr_vect, ya_vect)

                if (self.is_within_boundary(rya_angle,
                                            "min_rya_ang_ybond_inter", ge)):

                    # If no heavy atom is bonded to the acceptor, it means that
                    # only hydrogens may be bound to it. Then, we do not
                    # calculate the angles because hydrogens are too dynamic,
                    # i.e., the acceptor could be ionized or not at a specific
                    # moment in time and its hydrogens may be positioned in
                    # different ways.
                    if len(n_coords) == 0:
                        params = {"ya_dist_ybond_inter": ya_dist,
                                  "rya_ang_ybond_inter": rya_angle,
                                  "yan_ang_ybond_inter": -1}

                        inter = InteractionType(donor_grp,
                                                acceptor_grp,
                                                "Chalcogen bond",
                                                directional=True,
                                                params=params)
                        interactions.append(inter)
                    else:
                        # AY vector is always the same.
                        # Obs: the donor_grp is the Chalcogen (Y).
                        ay_vect = donor_grp.centroid - acceptor_grp.centroid

                        lowest_yan_angle = None
                        for n_coord in n_coords:
                            an_vect = n_coord - acceptor_grp.centroid
                            yan_angle = im.angle(ay_vect, an_vect)

                            # Update the XAR angle with the lowest value.
                            if (lowest_yan_angle is None
                                    or yan_angle < lowest_yan_angle):
                                lowest_yan_angle = yan_angle

                        # The angle will be None when any R (heavy atom)
                        # atom was found. In this case, the criteria must
                        # always fail.
                        if lowest_yan_angle is None:
                            lowest_yan_angle = -1

                        prop_name = "min_yan_ang_ybond_inter"
                        if self.is_within_boundary(lowest_yan_angle,
                                                   prop_name, ge):
                            params = {"ya_dist_ybond_inter": ya_dist,
                                      "rya_ang_ybond_inter": rya_angle,
                                      "yan_ang_ybond_inter": lowest_yan_angle}

                            inter = InteractionType(donor_grp,
                                                    acceptor_grp,
                                                    "Chalcogen bond",
                                                    directional=True,
                                                    params=params)
                            interactions.append(inter)
        return interactions

    @staticmethod
    def calc_chalc_bond_pi(self, params):
        """Default method to calculate chalcogen bonds between chalcogens and
        aromatic rings.

        Parameters
        ----------
        params : tuple of (:class:`~luna.mol.groups.AtomGroup`,\
                           :class:`~luna.mol.groups.AtomGroup`,\
                           :class:`~luna.mol.features.ChemicalFeature`,\
                           :class:`~luna.mol.features.ChemicalFeature`)
            The tuple follows the order (:math:`A`, :math:`B`, :math:`A_f`,
            :math:`B_f`), where :math:`A` and :math:`B` are two
            :class:`~luna.mol.groups.AtomGroup` objects, and :math:`A_f` and
            :math:`B_f` are their features
            (:class:`~luna.mol.features.ChemicalFeature` objects),
            respectively.

        Returns
        -------
         : list
        """
        if not self.add_non_cov:
            return []

        group1, group2, feat1, feat2 = params
        interactions = []

        if (feat1.name == "Aromatic" and feat2.name == "ChalcogenDonor"):
            donor_grp = group2
            ring_grp = group1
        else:
            donor_grp = group1
            ring_grp = group2

        # There are always just one donor/acceptor atom.
        donor_atm = donor_grp.atoms[0]

        # Interaction model: R-Y --- A, where A is the ring center
        # Defining the YA distance.
        ya_dist = im.euclidean_distance(donor_grp.centroid, ring_grp.centroid)

        if (self.is_within_boundary(ya_dist, "bsite_cutoff", le)
                and self.is_within_boundary(ya_dist,
                                            "max_yc_dist_ybond_inter", le)):

            ay_vect = donor_grp.centroid - ring_grp.centroid
            disp_angle = im.to_quad1(im.angle(ring_grp.normal, ay_vect))

            if (self.is_within_boundary(disp_angle,
                                        "max_disp_ang_ybond_inter", le)):
                # Interaction model: R-Y ---- A, where A is the ring center
                # YA vector is always the same
                ya_vect = ring_grp.centroid - donor_grp.centroid

                # Defining the angle RYA, where A is the ring center
                r_atms = [nbi for nbi in donor_atm.neighbors_info
                          if nbi.atomic_num != 1]
                r_elems = sorted([nbi.atomic_num
                                  for nbi in donor_atm.neighbors_info
                                  if nbi.atomic_num != 1])

                # Isothiazoles have only one sigma-hole located on the oposite
                # site of the N. Therefore, we should only evaluate this
                # sigma-hole by ignoring the angle formed with the carbon.
                # Beno et al (2015). DOI: https://doi.org/10.1021/jm501853m.
                ignore_carbon = False
                if len(r_elems) == 2 and r_elems[0] == 6 and r_elems[1] == 7:
                    ignore_carbon = True

                for r_atm in r_atms:
                    # Isothiazoles have only one sigma-hole located on the
                    # oposite site of the N. So, we must ignore the Carbon.
                    # Beno et al (2015).
                    #   DOI: https://doi.org/10.1021/jm501853m.
                    if r_atm.atomic_num == 6 and ignore_carbon is True:
                        continue

                    yr_vect = r_atm.coord - donor_grp.centroid
                    rya_angle = im.angle(yr_vect, ya_vect)

                    if self.is_within_boundary(rya_angle,
                                               "min_rya_ang_ybond_inter",
                                               ge):
                        params = {"yc_dist_ybond_inter": ya_dist,
                                  "disp_ang_ybond_inter": disp_angle,
                                  "rya_ang_ybond_inter": rya_angle}

                        inter = InteractionType(donor_grp,
                                                ring_grp,
                                                "Chalcogen-pi",
                                                directional=True,
                                                params=params)
                        interactions.append(inter)

        return interactions

    @staticmethod
    def calc_hbond(self, params):
        """Default method to calculate hydrogen bonds.

        Parameters
        ----------
        params : tuple of (:class:`~luna.mol.groups.AtomGroup`,\
                           :class:`~luna.mol.groups.AtomGroup`,\
                           :class:`~luna.mol.features.ChemicalFeature`,\
                           :class:`~luna.mol.features.ChemicalFeature`)
            The tuple follows the order (:math:`A`, :math:`B`, :math:`A_f`,
            :math:`B_f`), where :math:`A` and :math:`B` are two
            :class:`~luna.mol.groups.AtomGroup` objects, and :math:`A_f`
            and :math:`B_f` are their features
            (:class:`~luna.mol.features.ChemicalFeature` objects),
            respectively.

        Returns
        -------
         : list
        """
        if not self.add_non_cov:
            return []

        group1, group2, feat1, feat2 = params
        interactions = []

        if len(group1.atoms) != 1 or len(group2.atoms) != 1:
            logger.warning("One or more invalid atom groups were informed: "
                           "%s and %s. In hydrogen bonds, donor and acceptor "
                           "groups should always contain only one atom."
                           % (group1, group2))
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
        da_dist = im.euclidean_distance(donor_grp.centroid,
                                        acceptor_grp.centroid)

        if (self.is_within_boundary(da_dist, "bsite_cutoff", le)
                and self.is_within_boundary(da_dist,
                                            "max_da_dist_hb_inter", le)):

            # Interaction model: D-H ---- A-R.
            # Recover only hydrogen coordinates bonded to the donor.
            hydrog_coords = [nbi.coord
                             for nbi in donor_atm.neighbors_info
                             if nbi.atomic_num == 1]

            # Interaction model: D-H ---- A-R.
            # R coordinates, in which R is a heavy atom.
            r_coords = [nbi.coord for nbi in acceptor_atm.neighbors_info
                        if nbi.atomic_num != 1]

            # Firstly, it checks if it is not necessary to apply a strict
            # hbond rule, i.e., hydrogens must exist and all geometrical
            # criteria should be evaluated. Then it checks if no hydrogen is
            # bonded to the donor, or if the donor has hydrogens and only
            # hydrogens as neighbours (water, solvents, ammonia, SH2).
            # In the latter case, the hydrogens can be positioned in many
            # different ways, and each run of a tool like OpenBabel would vary
            # the hydrogen bond list when one applies this algorithm.
            #
            # If the user has also defined a list of lazy compounds, we can
            # skip the application of strict rules on them as well.
            # By default, the list contains only Water molecules.
            if ((self.strict_donor_rules is False and (len(hydrog_coords) == 0
                                                       or len(hydrog_coords) == len(donor_atm.neighbors_info)))
                    or (donor_atm.parent.resname in self.lazy_comps_list)):

                # When the position of the hydrogen cannot be defined, it
                # assumes the hydrogen to be located 1A away from the donor
                # in a line formed by the donor and the acceptor.
                ha_dist = da_dist - 1
                if self.is_within_boundary(ha_dist, "max_ha_dist_hb_inter",
                                           le):

                    # If no heavy atom is bonded to the acceptor, it means that
                    # only hydrogens may be bound to it. Then, we do not
                    # calculate the angles because hydrogens are too dynamic,
                    # i.e., the acceptor could be ionized or not at a specific
                    # moment in time and its hydrogens may be positioned in
                    # different ways.
                    if len(r_coords) == 0:
                        params = {"da_dist_hb_inter": da_dist,
                                  "ha_dist_hb_inter": -1,
                                  "dha_ang_hb_inter": -1,
                                  "har_ang_hb_inter": -1,
                                  "dar_ang_hb_inter": -1}

                        inter = InteractionType(donor_grp,
                                                acceptor_grp,
                                                "Hydrogen bond",
                                                directional=True,
                                                params=params)
                        interactions.append(inter)
                    else:
                        # AD vector is always the same.
                        ad_vect = donor_grp.centroid - acceptor_grp.centroid

                        lowest_dar_angle = None
                        for r_coord in r_coords:
                            ar_vect = r_coord - acceptor_grp.centroid
                            dar_angle = im.angle(ad_vect, ar_vect)

                            # Update the DAR angle with the lowest value.
                            if (lowest_dar_angle is None
                                    or dar_angle < lowest_dar_angle):
                                lowest_dar_angle = dar_angle

                        # The angle will be None when any R (heavy atom) atom
                        # was found. In this case, the criteria must always
                        # fail.
                        if lowest_dar_angle is None:
                            lowest_dar_angle = -1

                        if self.is_within_boundary(lowest_dar_angle,
                                                   "min_dar_ang_hb_inter",
                                                   ge):
                            params = {"da_dist_hb_inter": da_dist,
                                      "ha_dist_hb_inter": -1,
                                      "dha_ang_hb_inter": -1,
                                      "har_ang_hb_inter": -1,
                                      "dar_ang_hb_inter": lowest_dar_angle}

                            inter = InteractionType(donor_grp,
                                                    acceptor_grp,
                                                    "Hydrogen bond",
                                                    directional=True,
                                                    params=params)
                            interactions.append(inter)
            else:
                # It may happen that D is covalently bound to more than one
                # hydrogen atom. In this case, it is necessary to check the
                # distances and angles for each atom. It will produce a
                # hydrogen bond for each valid hydrogen.
                for h_coord in hydrog_coords:
                    ha_dist = \
                        im.euclidean_distance(h_coord, acceptor_grp.centroid)

                    hd_vect = donor_grp.centroid - h_coord
                    ha_vect = acceptor_grp.centroid - h_coord
                    dha_angle = im.angle(hd_vect, ha_vect)

                    if (self.is_within_boundary(ha_dist,
                                                "max_ha_dist_hb_inter", le)
                            and self.is_within_boundary(dha_angle,
                                                        "min_dha_ang_hb_inter",
                                                        ge)):

                        # If no heavy atom is bonded to the acceptor, it means
                        # that only hydrogens may be bound to it. Then, we do
                        # not calculate the angles because hydrogens are too
                        # dynamic, i.e., the acceptor could be ionized or not
                        # at a specific moment in time and its hydrogens may be
                        # positioned in different ways.
                        if len(r_coords) == 0:
                            params = {"da_dist_hb_inter": da_dist,
                                      "ha_dist_hb_inter": ha_dist,
                                      "dha_ang_hb_inter": dha_angle,
                                      "har_ang_hb_inter": -1,
                                      "dar_ang_hb_inter": -1}

                            inter = InteractionType(donor_grp,
                                                    acceptor_grp,
                                                    "Hydrogen bond",
                                                    directional=True,
                                                    params=params)
                            interactions.append(inter)
                        else:
                            # Interaction model: D-H ---- A-R
                            # AH vector is always the same.
                            ah_vect = h_coord - acceptor_grp.centroid
                            # AD vector is always the same.
                            ad_vect = \
                                donor_grp.centroid - acceptor_grp.centroid

                            # Interaction model: D-H ---- A-R
                            # Check the angles formed at the acceptor. When
                            # A is covalently bonded to more than one R atom,
                            # it is necessary to evaluate all possible angles
                            # D-A-R and H-A-R. In this case, all angles should
                            # satisfy the angle criterion. To do so, we could
                            # analyze only the lowest D-A-R and H-A-R angles.
                            # It guarantees that all angles will satisfy the
                            # criteria. OBS: it may happen that each one of the
                            # angles would belong to a different R atom.
                            lowest_har_angle = None
                            lowest_dar_angle = None
                            for r_coord in r_coords:
                                ar_vect = r_coord - acceptor_grp.centroid
                                har_angle = im.angle(ah_vect, ar_vect)
                                dar_angle = im.angle(ad_vect, ar_vect)

                                # Update the HAR angle with the lowest value.
                                if (lowest_har_angle is None
                                        or har_angle < lowest_har_angle):
                                    lowest_har_angle = har_angle
                                # Update the DAR angle with the lowest value.
                                if (lowest_dar_angle is None
                                        or dar_angle < lowest_dar_angle):
                                    lowest_dar_angle = dar_angle

                            # The angles will be None when any R (heavy atom)
                            # atom was found. In this case, the criteria must
                            # always fail.
                            if (lowest_har_angle is None
                                    or lowest_dar_angle is None):
                                lowest_har_angle = -1
                                lowest_dar_angle = -1

                            prop_name1 = "min_har_ang_hb_inter"
                            prop_name2 = "min_dar_ang_hb_inter"
                            if (self.is_within_boundary(lowest_har_angle,
                                                        prop_name1, ge)
                                    and self.is_within_boundary(lowest_dar_angle,
                                                                prop_name2, ge)):

                                # Only the lowest D-A-R and H-A-R angles
                                # are provided.
                                params = {"da_dist_hb_inter": da_dist,
                                          "ha_dist_hb_inter": ha_dist,
                                          "dha_ang_hb_inter": dha_angle,
                                          "har_ang_hb_inter": lowest_har_angle,
                                          "dar_ang_hb_inter": lowest_dar_angle}

                                inter = InteractionType(donor_grp,
                                                        acceptor_grp,
                                                        "Hydrogen bond",
                                                        directional=True,
                                                        params=params)
                                interactions.append(inter)

        return interactions

    @staticmethod
    def calc_weak_hbond(self, params):
        """Default method to calculate weak hydrogen bonds.

        Parameters
        ----------
        params : tuple of (:class:`~luna.mol.groups.AtomGroup`,\
                           :class:`~luna.mol.groups.AtomGroup`,\
                           :class:`~luna.mol.features.ChemicalFeature`,\
                           :class:`~luna.mol.features.ChemicalFeature`)
            The tuple follows the order (:math:`A`, :math:`B`, :math:`A_f`,
            :math:`B_f`), where :math:`A` and :math:`B` are two
            :class:`~luna.mol.groups.AtomGroup` objects, and :math:`A_f` and
            :math:`B_f` are their features
            (:class:`~luna.mol.features.ChemicalFeature` objects),
            respectively.

        Returns
        -------
         : list
        """
        if not self.add_non_cov:
            return []

        group1, group2, feat1, feat2 = params
        interactions = []

        if len(group1.atoms) != 1 or len(group2.atoms) != 1:
            logger.warning("One or more invalid atom groups were informed: "
                           "%s and %s. In weak hydrogen bonds, weak donor and "
                           "(weak) acceptor groups should always contain only "
                           "one atom." % (group1, group2))
            return []

        if ((feat1.name == "Acceptor" or feat1.name == "WeakAcceptor")
                and feat2.name == "WeakDonor"):
            donor_grp = group2
            acceptor_grp = group1
        elif (feat1.name == "WeakDonor" and (feat2.name == "Acceptor"
                                             or feat2.name == "WeakAcceptor")):
            donor_grp = group1
            acceptor_grp = group2

        else:
            logger.warning("Weak hydrogen bond requires a weak donor and an "
                           "(weak) acceptor groups. However, the informed "
                           "groups have the features %s and %s."
                           % (group1.feature_names, group2.feature_names))
            return []

        donor_atm = donor_grp.atoms[0]
        acceptor_atm = acceptor_grp.atoms[0]

        da_dist = im.euclidean_distance(donor_grp.centroid,
                                        acceptor_grp.centroid)
        if (self.is_within_boundary(da_dist, "bsite_cutoff", le)
                and self.is_within_boundary(da_dist,
                                            "max_da_dist_whb_inter", le)):

            # Interaction model: D-H ---- A-R.
            # Recover only hydrogen coordinates bonded to the donor.
            hydrog_coords = [nbi.coord
                             for nbi in donor_atm.neighbors_info
                             if nbi.atomic_num == 1]

            # Interaction model: D-H ---- A-R.
            # R coordinates, in which R is a heavy atom.
            r_coords = [nbi.coord
                        for nbi in acceptor_atm.neighbors_info
                        if nbi.atomic_num != 1]

            # Firstly, it checks if it is not necessary to apply a strict rule,
            # i.e., hydrogens must exist and all geometrical criteria should be
            # evaluated. Then it checks if no hydrogen is bonded to the donor,
            # or if the donor has hydrogens and only hydrogens as neighbours.
            #
            # If the user has also defined a list of lazy compounds, we can
            # skip the application of strict rules on them as well. By default,
            # the list contains only water, ammonia, and ammonium ion.
            if ((self.strict_weak_donor_rules is False and (len(hydrog_coords) == 0
                                                            or len(hydrog_coords) == len(donor_atm.neighbors_info)))
                    or (donor_atm.parent.resname in self.lazy_comps_list)):

                # When the position of the hydrogen cannot be defined,
                # it assumes the hydrogen to be located 1A away from the donor
                # in a line formed by the donor and the acceptor.
                ha_dist = da_dist - 1
                if self.is_within_boundary(ha_dist,
                                           "max_ha_dist_whb_inter", le):

                    # If no heavy atom is bonded to the acceptor, it means that
                    # only hydrogens may be bound to it. Then, we do not
                    # calculate the angles because hydrogens are too dynamic,
                    # i.e., the acceptor could be ionized or not at a specific
                    # moment in time and its hydrogens may be positioned in
                    # different ways.
                    if len(r_coords) == 0:
                        params = {"da_dist_whb_inter": da_dist,
                                  "ha_dist_whb_inter": -1,
                                  "dha_ang_whb_inter": -1,
                                  "har_ang_whb_inter": -1,
                                  "dar_ang_whb_inter": -1}

                        inter = InteractionType(donor_grp,
                                                acceptor_grp,
                                                "Weak hydrogen bond",
                                                directional=True,
                                                params=params)
                        interactions.append(inter)
                    else:
                        # AD vector is always the same.
                        ad_vect = donor_grp.centroid - acceptor_grp.centroid

                        lowest_dar_angle = None
                        for r_coord in r_coords:
                            ar_vect = r_coord - acceptor_grp.centroid
                            dar_angle = im.angle(ad_vect, ar_vect)

                            # Update the DAR angle with the lowest value.
                            if (lowest_dar_angle is None
                                    or dar_angle < lowest_dar_angle):
                                lowest_dar_angle = dar_angle

                        # The angle will be None when any R (heavy atom) atom
                        # was found. In this case, the criteria must always
                        # fail.
                        if lowest_dar_angle is None:
                            lowest_dar_angle = -1

                        if self.is_within_boundary(lowest_dar_angle,
                                                   "min_dar_ang_whb_inter",
                                                   ge):
                            params = {"da_dist_whb_inter": da_dist,
                                      "ha_dist_whb_inter": -1,
                                      "dha_ang_whb_inter": -1,
                                      "har_ang_whb_inter": -1,
                                      "dar_ang_whb_inter": lowest_dar_angle}

                            inter = InteractionType(donor_grp,
                                                    acceptor_grp,
                                                    "Weak hydrogen bond",
                                                    directional=True,
                                                    params=params)
                            interactions.append(inter)
            else:
                # It may happen that D is covalently bound to more than one
                # hydrogen atom. In such cases, it's necessary to check the
                # distances and angles for each atom.
                for h_coord in hydrog_coords:
                    ha_dist = im.euclidean_distance(h_coord,
                                                    acceptor_grp.centroid)

                    hd_vect = donor_grp.centroid - h_coord
                    ha_vect = acceptor_grp.centroid - h_coord
                    dha_angle = im.angle(hd_vect, ha_vect)

                    prop_name1 = "max_ha_dist_whb_inter"
                    prop_name2 = "min_dha_ang_whb_inter"
                    if (self.is_within_boundary(ha_dist, prop_name1, le)
                            and self.is_within_boundary(dha_angle,
                                                        prop_name2, ge)):

                        # If no heavy atom is bonded to the acceptor, it means
                        # that only hydrogens may be bound to it. Then, we do
                        # not calculate the angles because hydrogens are too
                        # dynamic, i.e., the acceptor could be ionized or not
                        # at a specific moment in time and its hydrogens may
                        # be positioned in different ways.
                        if len(r_coords) == 0:
                            params = {"da_dist_whb_inter": da_dist,
                                      "ha_dist_whb_inter": ha_dist,
                                      "dha_ang_whb_inter": dha_angle,
                                      "har_ang_whb_inter": -1,
                                      "dar_ang_whb_inter": -1}

                            inter = InteractionType(donor_grp,
                                                    acceptor_grp,
                                                    "Weak hydrogen bond",
                                                    directional=True,
                                                    params=params)
                            interactions.append(inter)
                        else:
                            # Interaction model: D-H ---- A-R
                            # AH vector is always the same.
                            ah_vect = h_coord - acceptor_grp.centroid
                            # AD vector is always the same.
                            ad_vect = (donor_grp.centroid
                                       - acceptor_grp.centroid)

                            # Interaction model: D-H ---- A-R
                            # Check the angles formed at the acceptor.
                            # When A is covalently bonded to more than one R
                            # atom, it is necessary to evaluate all possible
                            # angles D-A-R and H-A-R. In this case, all angles
                            # should satisfy the angle criterion. To do so, we
                            # could analyze only the lowest D-A-R and H-A-R
                            # angles. It guarantees that all angles will
                            # satisfy the criteria.
                            #
                            # OBS: it may happen that each one of the angles
                            # would belong to a different R atom.
                            lowest_har_angle = None
                            lowest_dar_angle = None
                            for r_coord in r_coords:
                                ar_vect = r_coord - acceptor_grp.centroid
                                har_angle = im.angle(ah_vect, ar_vect)
                                dar_angle = im.angle(ad_vect, ar_vect)

                                # Update the HAR angle with the lowest value.
                                if (lowest_har_angle is None
                                        or har_angle < lowest_har_angle):
                                    lowest_har_angle = har_angle
                                # Update the DAR angle with the lowest value.
                                if (lowest_dar_angle is None
                                        or dar_angle < lowest_dar_angle):
                                    lowest_dar_angle = dar_angle

                            # The angles will be None when any R (heavy atom)
                            # atom was found. In this case, the criteria must
                            # always fail.
                            if (lowest_har_angle is None
                                    or lowest_dar_angle is None):
                                lowest_har_angle = -1
                                lowest_dar_angle = -1

                            prop_name1 = "min_har_ang_whb_inter"
                            prop_name2 = "min_dar_ang_whb_inter"
                            if (self.is_within_boundary(lowest_har_angle,
                                                        prop_name1, ge)
                                    and self.is_within_boundary(lowest_dar_angle,
                                                                prop_name2, ge)):

                                # Only the lowest D-A-R and H-A-R angles are
                                # provided.
                                params = {"da_dist_whb_inter": da_dist,
                                          "ha_dist_whb_inter": ha_dist,
                                          "dha_ang_whb_inter": dha_angle,
                                          "har_ang_whb_inter": lowest_har_angle,
                                          "dar_ang_whb_inter": lowest_dar_angle}

                                inter = InteractionType(donor_grp,
                                                        acceptor_grp,
                                                        "Weak hydrogen bond",
                                                        directional=True,
                                                        params=params)
                                interactions.append(inter)

        return interactions

    @staticmethod
    def calc_hbond_pi(self, params):
        """Default method to calculate hydrogen bonds between (weak) donors
        and aromatic rings.

        Parameters
        ----------
        params : tuple of (:class:`~luna.mol.groups.AtomGroup`,\
                           :class:`~luna.mol.groups.AtomGroup`,\
                           :class:`~luna.mol.features.ChemicalFeature`,\
                           :class:`~luna.mol.features.ChemicalFeature`)
            The tuple follows the order (:math:`A`, :math:`B`, :math:`A_f`,
            :math:`B_f`), where :math:`A` and :math:`B` are two
            :class:`~luna.mol.groups.AtomGroup` objects, and :math:`A_f` and
            :math:`B_f` are their features
            (:class:`~luna.mol.features.ChemicalFeature` objects),
            respectively.

        Returns
        -------
         : list
        """
        if not self.add_non_cov:
            return []

        group1, group2, feat1, feat2 = params
        interactions = []

        if (feat1.name == "Aromatic"
                and (feat2.name == "Donor" or feat2.name == "WeakDonor")):
            ring_grp = group1
            donor_grp = group2
        elif (feat2.name == "Aromatic"
                and (feat1.name == "Donor" or feat1.name == "WeakDonor")):
            ring_grp = group2
            donor_grp = group1

        else:
            logger.warning("Hydrogen bond involving pi-systems requires an "
                           "aromatic and donor (weak donor) groups. However, "
                           "the informed groups have the features %s and %s."
                           % (group1.feature_names, group2.feature_names))
            return []

        if len(donor_grp.atoms) != 1:
            logger.warning("Invalid (weak) donor group was informed: %s. "
                           "In hydrogen bonds involving pi-systems, (weak) "
                           "donor groups should always contain only one "
                           "atom." % donor_grp)
            return []

        # There are always just one donor/acceptor atom.
        donor_atm = donor_grp.atoms[0]

        # Interaction model: D-H ---- A, in which A is the ring center.
        da_dist = im.euclidean_distance(donor_grp.centroid, ring_grp.centroid)
        if (self.is_within_boundary(da_dist, "bsite_cutoff", le)
                and self.is_within_boundary(da_dist,
                                            "max_dc_dist_whb_inter", le)):

            # Interaction model: D-H ---- A, in which A is the ring center.
            # Recover only hydrogen coordinates bonded to the donor.
            hydrog_coords = [nbi.coord
                             for nbi in donor_atm.neighbors_info
                             if nbi.atomic_num == 1]

            # Firstly, it checks if it is not necessary to apply a strict
            # hbond rule, i.e., hydrogens must exist and all geometrical
            # criteria should be evaluated. Then it checks if no hydrogen
            # is bonded to the donor, or if the donor has hydrogens and
            # only hydrogens as neighbours (water, solvents, ammonia, SH2).
            # In the latter case, the hydrogens can be positioned in many
            # different set of ways, and each run of a tool like OpenBabel
            # would vary the hydrogen bond list when one applies this
            # algorithm.
            #
            # If the user has also defined a list of lazy compounds, we can
            # skip the application of strict rules on them as well.
            # By default, the list contains only water, ammonia, and
            # ammonium ion.
            if ((self.strict_weak_donor_rules is False
                    and (len(hydrog_coords) == 0
                         or len(hydrog_coords) == len(donor_atm.neighbors_info)))
                    or (donor_atm.parent.resname in self.lazy_comps_list)):

                # When the position of the hydrogen cannot be defined, it
                # assumes the hydrogen to be located 1A away from the donor
                # in a line formed by the donor and the acceptor.
                ha_dist = da_dist - 1
                if self.is_within_boundary(ha_dist,
                                           "max_hc_dist_whb_inter", le):

                    # Interaction model: D-H ---- A, in which A is the ring
                    # center. Calculate the displacement angle formed between
                    # the ring normal and the vector Donor-Centroid.
                    ad_vect = donor_grp.centroid - ring_grp.centroid
                    disp_angle = im.to_quad1(im.angle(ring_grp.normal,
                                                      ad_vect))

                    if (self.is_within_boundary(disp_angle,
                                                "max_disp_ang_whb_inter", le)):
                        params = {"dc_dist_whb_inter": da_dist,
                                  "hc_dist_whb_inter": -1,
                                  "dhc_ang_whb_inter": -1,
                                  "disp_ang_whb_inter": disp_angle}

                        inter = InteractionType(donor_grp,
                                                ring_grp,
                                                "Weak hydrogen bond",
                                                directional=True,
                                                params=params)
                        interactions.append(inter)
            else:
                # It may happen that D is covalently bound to more than one
                # hydrogen atom. In this case, it is necessary to check the
                # distances and angles for each atom. It will produce a
                # hydrogen bond for each valid hydrogen.
                for h_coord in hydrog_coords:
                    ha_dist = im.euclidean_distance(h_coord, ring_grp.centroid)

                    hd_vect = donor_grp.centroid - h_coord
                    ha_vect = ring_grp.centroid - h_coord
                    dha_angle = im.angle(hd_vect, ha_vect)

                    prop_name1 = "max_hc_dist_whb_inter"
                    prop_name2 = "min_dhc_ang_whb_inter"
                    if (self.is_within_boundary(ha_dist, prop_name1, le)
                            and self.is_within_boundary(dha_angle,
                                                        prop_name2, ge)):

                        # Interaction model: D-H ---- A, in which A is the ring
                        # center. Calculate the displacement angle formed
                        # between the ring normal and the vector
                        # Donor-Centroid.
                        ad_vect = donor_grp.centroid - ring_grp.centroid
                        disp_angle = im.to_quad1(im.angle(ring_grp.normal,
                                                          ad_vect))

                        if (self.is_within_boundary(disp_angle,
                                                    "max_disp_ang_whb_inter",
                                                    le)):
                            params = {"dc_dist_whb_inter": da_dist,
                                      "hc_dist_whb_inter": ha_dist,
                                      "dhc_ang_whb_inter": dha_angle,
                                      "disp_ang_whb_inter": disp_angle}

                            inter = InteractionType(donor_grp,
                                                    ring_grp,
                                                    "Weak hydrogen bond",
                                                    directional=True,
                                                    params=params)
                            interactions.append(inter)
        return interactions

    @staticmethod
    def calc_ionic(self, params):
        """Default method to calculate attractive ionic interactions.

        Parameters
        ----------
        params : tuple of (:class:`~luna.mol.groups.AtomGroup`,\
                           :class:`~luna.mol.groups.AtomGroup`,\
                           :class:`~luna.mol.features.ChemicalFeature`,\
                           :class:`~luna.mol.features.ChemicalFeature`)
            The tuple follows the order (:math:`A`, :math:`B`, :math:`A_f`,
            :math:`B_f`), where :math:`A` and :math:`B` are two
            :class:`~luna.mol.groups.AtomGroup` objects, and :math:`A_f` and
            :math:`B_f` are their features
            (:class:`~luna.mol.features.ChemicalFeature` objects),
            respectively.

        Returns
        -------
         : list
        """
        if not self.add_non_cov:
            return []

        group1, group2, feat1, feat2 = params
        interactions = []

        cc_dist = im.euclidean_distance(group1.centroid, group2.centroid)
        if (self.is_within_boundary(cc_dist, "bsite_cutoff", le)
                and self.is_within_boundary(cc_dist,
                                            "max_dist_attract_inter", le)):

            params = {"dist_attract_inter": cc_dist}
            inter = InteractionType(group1, group2, "Ionic", params=params)

            interactions.append(inter)
        return interactions

    @staticmethod
    def calc_repulsive(self, params):
        """Default method to calculate repulsive ionic interactions.

        Parameters
        ----------
        params : tuple of (:class:`~luna.mol.groups.AtomGroup`,\
                           :class:`~luna.mol.groups.AtomGroup`,\
                           :class:`~luna.mol.features.ChemicalFeature`,\
                           :class:`~luna.mol.features.ChemicalFeature`)
            The tuple follows the order (:math:`A`, :math:`B`, :math:`A_f`,
            :math:`B_f`), where :math:`A` and :math:`B` are two
            :class:`~luna.mol.groups.AtomGroup` objects, and :math:`A_f` and
            :math:`B_f` are their features
            (:class:`~luna.mol.features.ChemicalFeature` objects),
            respectively.

        Returns
        -------
         : list
        """
        if not self.add_non_cov:
            return []

        group1, group2, feat1, feat2 = params
        interactions = []

        cc_dist = im.euclidean_distance(group1.centroid, group2.centroid)
        if (self.is_within_boundary(cc_dist, "bsite_cutoff", le)
                and self.is_within_boundary(cc_dist,
                                            "max_dist_repuls_inter", le)):

            params = {"dist_repuls_inter": cc_dist}
            inter = InteractionType(group1, group2, "Repulsive", params=params)

            interactions.append(inter)
        return interactions

    @staticmethod
    def calc_proximal(self, params):
        """Default method to calculate proximal interactions.

        Parameters
        ----------
        params : tuple of (:class:`~luna.mol.groups.AtomGroup`,\
                           :class:`~luna.mol.groups.AtomGroup`,\
                           :class:`~luna.mol.features.ChemicalFeature`,\
                           :class:`~luna.mol.features.ChemicalFeature`)
            The tuple follows the order (:math:`A`, :math:`B`, :math:`A_f`,
            :math:`B_f`), where :math:`A` and :math:`B` are two
            :class:`~luna.mol.groups.AtomGroup` objects, and :math:`A_f` and
            :math:`B_f` are their features
            (:class:`~luna.mol.features.ChemicalFeature` objects),
            respectively.

        Returns
        -------
         : list
        """
        if not self.add_proximal:
            return []

        group1, group2, feat1, feat2 = params
        interactions = []

        cc_dist = im.euclidean_distance(group1.centroid, group2.centroid)
        if (self.is_within_boundary(cc_dist, "min_dist_proximal", ge)
                and self.is_within_boundary(cc_dist, "max_dist_proximal", le)):

            params = {"dist_proximal": cc_dist}
            inter = InteractionType(group1, group2, "Proximal", params=params)
            interactions.append(inter)

        return interactions

    @staticmethod
    def calc_metal_coord(self, params):
        # Ignore covalent bonds.
        if not self.add_cov:
            return []

        group1, group2, feat1, feat2 = params
        interactions = []

        if feat1.name == "Atom" and feat2.name == "Metal":
            other_grp = group1
            metal_grp = group2
        elif feat1.name == "Metal" and feat2.name == "Atom":
            other_grp = group2
            metal_grp = group1
        else:
            logger.warning("Metal complexes require an "
                           "other and a metal group. However, "
                           "the informed groups have the features %s and %s."
                           % (group1.feature_names, group2.feature_names))
            return []

        if other_grp.atoms[0].element not in ["O", "N", "S"]:
            return []

        # Interaction model: M ---- A.
        # MA distance
        ma_dist = im.euclidean_distance(metal_grp.centroid,
                                        other_grp.centroid)

        if self.is_within_boundary(ma_dist, "max_ma_dist_metal_coord", le):
            params = {"ma_dist_metal_coord": ma_dist}
            inter = InteractionType(other_grp,
                                    metal_grp,
                                    "Metal coordination",
                                    directional=True,
                                    params=params)
            interactions.append(inter)

        return interactions

    @staticmethod
    def calc_atom_atom(self, params):
        """Default method to calculate atom-atom interactions, which include
        covalent bonds, Van der Waals, Van der Waals clash, and atom overlap.

        Note that covalent bonds are controlled by the flag ``add_cov``, while
        the other three interactions are controlled by the flag
        ``add_atom_atom``.

        .. note::
            We opted to separate `Van der Waals` from other non-covalent
            interactions because LUNA may generate an unnecessary number of
            additional interactions that are usually already represented by
            other non-covalent interactions as weak hydrogen bonds,
            hydrophobic, or dipole-dipole interactions. Thus, to give users a
            fine-grain control over which interactions to calculate, we
            provided this additional flag to turn off the calculation of
            Van der Waals interactions.

        Parameters
        ----------
        params : tuple of (:class:`~luna.mol.groups.AtomGroup`,\
                           :class:`~luna.mol.groups.AtomGroup`,\
                           :class:`~luna.mol.features.ChemicalFeature`,\
                           :class:`~luna.mol.features.ChemicalFeature`)
            The tuple follows the order (:math:`A`, :math:`B`, :math:`A_f`,
            :math:`B_f`), where :math:`A` and :math:`B` are two
            :class:`~luna.mol.groups.AtomGroup` objects, and :math:`A_f` and
            :math:`B_f` are their features
            (:class:`~luna.mol.features.ChemicalFeature` objects),
            respectively.

        Returns
        -------
         : list
        """
        group1, group2, feat1, feat2 = params
        interactions = []

        atm1 = group1.atoms[0]
        atm2 = group2.atoms[0]

        cc_dist = im.euclidean_distance(group1.centroid, group2.centroid)
        params = {"dist_atom_atom": cc_dist}

        # It checks if the two atoms are neighbors, i.e., if they are
        # covalently bonded. The covalent bonds are detected by OpenBabel,
        # which besides other evaluations, states that two atoms are
        # covalently bonded if:
        #       0.4 <= d(a1, a2) <= cov_rad(a1) + cov_rad(a2) + 0.45
        if atm1.is_neighbor(atm2):
            # Ignore covalent bonds.
            if not self.add_cov:
                return []

            bond_type = atm1.get_neighbor_info(atm2).bond_type

            if bond_type in COV_BONDS_MAPPING:
                bond_name = COV_BONDS_MAPPING[bond_type]
            else:
                bond_name = "Other bond"
                logger.warning("An unexpected bond (%s) was found between the "
                               "atoms %s and %s. Therefore, a general bond "
                               "name (%s) will be used instead."
                               % (bond_type, atm1, atm2, bond_name))

            inter = InteractionType(group1, group2, bond_name, params=params)
            interactions.append(inter)
        else:
            if not self.add_atom_atom:
                return []

            cov1 = ob.GetCovalentRad(ob.GetAtomicNum(atm1.element))
            cov2 = ob.GetCovalentRad(ob.GetAtomicNum(atm2.element))

            if cc_dist <= cov1 + cov2:
                inter = InteractionType(group1,
                                        group2,
                                        "Atom overlap",
                                        params=params)
                interactions.append(inter)
            else:
                rdw1 = ob.GetVdwRad(ob.GetAtomicNum(atm1.element))
                rdw2 = ob.GetVdwRad(ob.GetAtomicNum(atm2.element))

                # r1 + r2 - d < 0 => no clash
                # r1 + r2 - d = 0 => in the limit, i.e., spheres are touching.
                # r1 + r2 - d > 0 => clash.
                if ((rdw1 + rdw2 - cc_dist)
                        >= self.inter_config.get("vdw_clash_tolerance", 0)):

                    # Ignore Van der Waals and clashes for atoms separated from
                    # each other by only N bonds. Covalent bonds keep atoms
                    # very tightly, producing distances lower than their sum of
                    # Van der Waals radius. As a consequence the algorithm will
                    # find a lot of false clashes and Van der Waals
                    # interactions.
                    #
                    # It is better to keep this function inside the IFs to
                    # avoid the Dijkstra processing for pairs of atoms that
                    # wouldn't enter inside the IF.
                    min_bond_sep = \
                        self.inter_config.get("min_bond_separation", 0)
                    shortest_path_length = \
                        group1.get_shortest_path_length(group2, min_bond_sep)

                    # If get_shortest_path_length() returns any value that is
                    # not infinite (INF), it means these two groups contain a
                    # path with at less than or equal to the cutoff
                    # 'min_bond_separation'. Therefore, ignore them.
                    if shortest_path_length != float('inf'):
                        return []

                    inter = InteractionType(group1,
                                            group2,
                                            "Van der Waals clash",
                                            params=params)
                    interactions.append(inter)

                elif (cc_dist
                        <= (rdw1 + rdw2
                            + self.inter_config.get("vdw_tolerance", 0))):

                    # Ignore Van der Waals and clashes for atoms separated from
                    # each other by only N bonds. Covalent bonds keep atoms
                    # very tightly, producing distances lower than their sum of
                    # Van der Waals radius. As a consequence the algorithm will
                    # find a lot of false clashes and Van der Waals
                    # interactions.
                    #
                    # It is better to keep this function inside the IFs to
                    # avoid the Dijkstra processing for pairs of atoms that
                    # wouldn't enter inside the IF.
                    min_bond_sep = \
                        self.inter_config.get("min_bond_separation", 0)
                    shortest_path_length = \
                        group1.get_shortest_path_length(group2, min_bond_sep)

                    # If get_shortest_path_length() returns any value that is
                    # not infinite (INF), it means these two groups contain a
                    # path with at less than or equal to the cutoff
                    # 'min_bond_separation'. Therefore, ignore them.
                    if shortest_path_length != float('inf'):
                        return []

                    inter = InteractionType(group1,
                                            group2,
                                            "Van der Waals",
                                            params=params)
                    interactions.append(inter)

        return interactions

    def _is_intramol_inter(self, grp1, grp2):
        comps1 = grp1.compounds
        comps2 = grp2.compounds
        return len(comps1) == 1 and len(comps2) == 1 and comps1 == comps2

    def is_within_boundary(self, value, key, func):
        """Check if a value is within the boundary defined for a given
        parameter.

        .. note::
            It will always return True if the parameter does not exist in
            ``inter_config``.

        Parameters
        ----------
        value : any
            The value to be evaluated.
        key :
            A parameter defined in ``inter_config``.
        func : callable
            The function that evaluates if ``value`` is within the
            boundaries defined for the parameter ``key``.

            Usually, the comparison functions (e.g., lt, le, ge, gt)
            available in the Python module :py:mod:`operator` are enough for
            number comparisons. If you need custom comparison functions,
            just provide them here.

        Returns
        -------
         : bool
        """
        if key not in self.inter_config:
            return True
        return func(value, self.inter_config[key])

    def is_feature_pair_valid(self, feat1, feat2):
        """Check if the provided pair of features is valid or not.

        It will be valid if the pair exists in ``funcs``, i.e., there is one or
        more functions to calculate interactions defined for that given pair of
        features.

        It also return False if non-covalent interactions is turned off
        (``add_non_cov = False``) and at least one of the features is not
        `Atom`. This is useful to save processing time as it skips pairs that
        have functions to calculate non-covalent interactions right away.

        Parameters
        ----------
        feat1, feat2: :class:`~luna.mol.features.ChemicalFeature`

        Returns
        -------
         : bool
        """
        if isinstance(feat1, ChemicalFeature):
            feat1 = feat1.name
        if isinstance(feat2, ChemicalFeature):
            feat2 = feat2.name

        if not self.add_non_cov and (feat1 != "Atom" or feat2 != "Atom"):
            return False

        funcs = self.funcs
        return (True if ((feat1, feat2) in funcs
                         or (feat2, feat1) in funcs) else False)

    def get_functions(self, feat1, feat2):
        """Get the functions to calculate interactions for the given features.

        Parameters
        ----------
        feat1, feat2: :class:`~luna.mol.features.ChemicalFeature`

        Returns
        -------
         : iterable of callable
        """

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

    def set_functions_to_pair(self, pair, funcs):
        """Set functions to calculate interaction for the given pair of
        features.

        Parameters
        ----------
        pair: tuple of (:class:`~luna.mol.features.ChemicalFeature`, \
                    :class:`~luna.mol.features.ChemicalFeature`)
        funcs : iterable of callable
        """
        self.funcs[pair] = funcs
