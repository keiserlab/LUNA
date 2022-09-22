from collections import defaultdict
from itertools import chain

import numpy as np

import networkx as nx
from networkx.algorithms.shortest_paths.weighted import single_source_dijkstra

from Bio.KDTree import KDTree

from luna.MyBio.selector import Selector, AtomSelector
from luna.MyBio.util import biopython_entity_to_mol
from luna.interaction.contact import get_proximal_compounds
from luna.interaction.type import InteractionType
from luna.mol.atom import ExtendedAtom, AtomData
from luna.mol.charge_model import OpenEyeModel
from luna.mol.features import ChemicalFeature
from luna.wrappers.base import MolWrapper
from luna.util.exceptions import MoleculeSizeError, IllegalArgumentError
from luna.util.default_values import ACCEPTED_MOL_OBJ_TYPES, COV_SEARCH_RADIUS
from luna.util import math as im
from luna.util.file import pickle_data, unpickle_data
from luna.version import __version__

import logging
logger = logging.getLogger()


SS_BOND_FEATURES = ["Atom", "Acceptor", "ChalcogenDonor", "Hydrophobic"]


class AtomGroupsManager():
    """Store and manage `AtomGroup` objects.

    Parameters
    ----------
    atm_grps : iterable of `AtomGroup`, optional
        An initial sequence of `AtomGroup` objects.
    entry : :class:`~luna.mol.entry.Entry`, optional
        The chain or molecule from where the atom groups were perceived.

    Attributes
    ----------
    entry : :class:`~luna.mol.entry.Entry`
        The chain or molecule from where the atom groups were perceived.
    graph : :py:class:`networkx.Graph`
        Represent ``entry`` as a graph and its vicinity.
    version : str
        The LUNA version when the object was created.
    """

    def __init__(self, atm_grps=None, entry=None):
        self._atm_grps = []
        self.entry = entry
        self._child_dict = {}

        self.graph = nx.Graph()

        self.version = __version__

        self._compounds = set()

        self.add_atm_grps(atm_grps)

    @property
    def atm_grps(self):
        """iterable of `AtomGroup`, read-only: The sequence of `AtomGroup` objects.\
        Additional objects should be added using the method :py:meth:`add_atm_grps`."""
        return self._atm_grps

    @property
    def compounds(self):
        """set of :class:`~luna.MyBio.PDB.Residue.Residue`, read-only:
        Compounds comprising the ``atm_grps``."""
        return self._compounds

    @property
    def child_dict(self):
        """dict, read-only: Mapping between atoms (`ExtendedAtom`) and atom groups (`AtomGroup`).

        The mapping is a dict of {tuple of `ExtendedAtom` instances : `AtomGroup`} and is automatically
        updated when :py:meth:`add_atm_grps` is called."""
        return self._child_dict

    @property
    def size(self):
        """int, read-only: The number of atom groups in ``atm_grps``."""
        return len(self._atm_grps)

    @property
    def summary(self):
        """dict, read-only: The number of physicochemical features in ``atm_grps``."""
        summary = defaultdict(int)
        for grp in self.atm_grps:
            for feature in grp.features:
                summary[feature] += 1

        return summary

    def find_atm_grp(self, atoms):
        """ Find the atom group that contains the sequence of atoms ``atoms``.

        Returns
        -------
         : `AtomGroup` or None
            An atom group object or None if ``atoms`` is not in the ``child_dict`` mapping.
        """
        return self.child_dict.get(tuple(sorted(atoms)), None)

    def get_all_interactions(self):
        """Return all interactions established by the atom groups in ``atm_grps``.

        Returns
        -------
         : set of :class:`~luna.interactions.type.InteractionType`
            All interactions.
        """
        return set(chain.from_iterable([atm_grp.interactions for atm_grp in self.atm_grps]))

    def apply_filter(self, func):
        """Apply a filtering function over the atom groups in ``atm_grps`.

        Parameters
        ----------
        func : callable
            A filtering function that returns True case an `AtomGroup` object is valid and False otherwise.

        Yields
        ------
        `AtomGroup`
            A valid `AtomGroup` object.
        """
        for atm_grp in self.atm_grps:
            if func(atm_grp):
                yield atm_grp

    def filter_by_types(self, types, must_contain_all=True):
        """Filter `AtomGroup` objects by their physicochemical features.

        Parameters
        ----------
        types : iterable of str
            A sequence of physicochemical features.
        must_contain_all : bool
            If True, an `AtomGroup` object should contain all physicochemical features in ``types`` to be accepted.
            Otherwise, it will be filtered out.

        Yields
        ------
        `AtomGroup`
            A valid `AtomGroup` object.
        """
        for atm_grp in self.atm_grps:
            if must_contain_all:
                if set(types).issubset(set(atm_grp.feature_names)):
                    yield atm_grp
            else:
                if len(set(types) & set(atm_grp.feature_names)) > 0:
                    yield atm_grp

    def add_atm_grps(self, atm_grps):
        """Add one or more `AtomGroup` objects to ``atm_grps`` and automatically update ``child_dict``."""

        atm_grps = atm_grps or []

        self._atm_grps = list(set(self._atm_grps + list(atm_grps)))

        for atm_grp in atm_grps:
            self.child_dict[tuple(sorted(atm_grp.atoms))] = atm_grp
            self._compounds.update(atm_grp.compounds)
            atm_grp.manager = self

    def remove_atm_grps(self, atm_grps):
        """Remove one or more `AtomGroup` objects from ``atm_grps`` and automatically update ``child_dict``.

        Any recursive references to the removed objects will also be cleared.
        """
        self._atm_grps = list(set(self._atm_grps) - set(atm_grps))

        for atm_grp in atm_grps:
            # ExtendedAtom objects keep a list of all AtomGroup objects to which they belong to. So, if we don't clear the references
            # directly, the AtomGroup objects will still exist in the ExtendedAtom list even when they were already removed from an
            # instance of AtomGroupsManager().
            atm_grp.clear_refs()

            # Remove the atom group from the dict.
            key = tuple(sorted(atm_grp.atoms))
            if key in self.child_dict:
                del self.child_dict[key]

    def new_atm_grp(self, atoms, features=None, interactions=None):
        """Create a new `AtomGroup` object for ``atoms`` if one does not exist yet.
        Otherwise, return the existing `AtomGroup` object.

        Parameters
        ----------
        atoms : iterable of :class:`~luna.mol.atom.ExtendedAtom`
            A sequence of atoms.
        features : iterable of :class:`~luna.mol.features.ChemicalFeature`, optional
            If provided, add ``features`` to a new or an already existing `AtomGroup` object.
        interactions : iterable of :class:`~luna.interaction.type.InteractionType`, optional
            If provided, add ``interactions`` to a new or an already existing `AtomGroup` object.

        Returns
        -------
         : AtomGroup
            A new or an already existing `AtomGroup` object.
        """
        key = tuple(sorted(atoms))
        features = features or []
        interactions = interactions or []

        if key in self.child_dict:
            atm_grp = self.child_dict[key]
        else:
            atm_grp = AtomGroup(atoms)
            self.add_atm_grps([atm_grp])

        if features:
            atm_grp.add_features(features)

        if interactions:
            atm_grp.add_interactions(interactions)

        return atm_grp

    def merge_hydrophobic_atoms(self, interactions_mngr):
        """Create hydrophobic islands by merging covalently bonded hydrophobic atoms in ``atm_grps``.
        Hydrophobic islands are atom groups having the feature `Hydrophobe`.

        Atom-atom hydrophobic interactions in ``interactions_mngr`` are also converted to
        island-island interactions.

        Parameters
        ----------
        interactions_mngr : :class:`~luna.interaction.calc.InteractionsManager`
            An :class:`~luna.interaction.calc.InteractionsManager` object from where
            hydrophobic interactions are selected and convert from atom-atom to island-island interactions.
        """

        # Only hydrophobic atom groups.
        hydrop_atm_grps = list(self.filter_by_types(["Hydrophobic"]))

        # Hydrophobic islands dictionary. Keys are integer values and items are defined by a set of atom groups.
        hydrop_islands = defaultdict(set)

        # It stores a mapping of an atom (represented by its full id) and a hydrophobic island (defined by its keys).
        atm_mapping = {}

        island_id = 0
        for atm_grp in hydrop_atm_grps:
            # Hydrophobic atoms are defined always as only one atom.
            atm = atm_grp.atoms[0]

            # Recover the groups of all neighbors of this atom (it will merge all existing islands).
            nb_grps = set([atm_mapping[nbi.full_id] for nbi in atm.neighbors_info if nbi.full_id in atm_mapping])

            # Already there are hydrophobic islands formed by the neighbors of this atom.
            if nb_grps:
                # Merge all groups of the neighbors of this atom.
                new_island = set(chain.from_iterable([hydrop_islands.pop(nb_grp_id) for nb_grp_id in nb_grps]))
                # Include this atom to the merged group.
                new_island.add(atm)

                for k in atm_mapping:
                    if atm_mapping[k] in nb_grps:
                        atm_mapping[k] = island_id

                hydrop_islands[island_id] = new_island
                atm_mapping[atm.get_full_id()] = island_id
            else:
                atm_mapping[atm.get_full_id()] = island_id
                hydrop_islands[island_id].add(atm)

            island_id += 1

        # Create AtomGroup objects for the hydrophobic islands
        for island_id in hydrop_islands:
            # It will update an existing atom group or create a new one with the informed parameters.
            hydrophobe = self.new_atm_grp(hydrop_islands[island_id], [ChemicalFeature("Hydrophobe")])
            # Update the island information
            hydrop_islands[island_id] = hydrophobe

        hydrop_interactions = list(interactions_mngr.filter_by_types(["Hydrophobic"]))
        island_island_inter = defaultdict(set)
        for inter in hydrop_interactions:
            src_atm = inter.src_grp.atoms[0]
            trgt_atm = inter.trgt_grp.atoms[0]

            # The two island ids are used as key.
            key = tuple(sorted([atm_mapping[src_atm.get_full_id()], atm_mapping[trgt_atm.get_full_id()]]))

            island_island_inter[key].add(inter)

        interactions = set()
        for k in island_island_inter:
            island_atms = defaultdict(set)
            for inter in island_island_inter[k]:
                src_atm = inter.src_grp.atoms[0]
                trgt_atm = inter.trgt_grp.atoms[0]

                island_atms[atm_mapping[src_atm.get_full_id()]].add(src_atm)
                island_atms[atm_mapping[trgt_atm.get_full_id()]].add(trgt_atm)

            centroid1 = im.centroid(im.atom_coordinates(island_atms[k[0]]))
            centroid2 = im.centroid(im.atom_coordinates(island_atms[k[1]]))
            cc_dist = im.euclidean_distance(centroid1, centroid2)

            params = {"dist_hydrop_inter": cc_dist}

            inter = InteractionType(hydrop_islands[k[0]], hydrop_islands[k[1]], "Hydrophobic",
                                    src_interacting_atms=island_atms[k[0]], trgt_interacting_atms=island_atms[k[1]], params=params)
            interactions.add(inter)

        # Update the list of interactions with the new island-island interactions.
        interactions_mngr.add_interactions(interactions)
        # Remove atom-atom hydrophobic interactions.
        interactions_mngr.remove_interactions(hydrop_interactions)

        for atm_grp in hydrop_atm_grps:
            features = [f for f in atm_grp.features if f.name != "Hydrophobic"]

            # It may happen that a atom group ends up having no feature after the remotion of the feature "Hydrophobic".
            # This is unlikely to occur as all atoms (by default) will have at least the feature 'Atom'. But, depending one the
            # pharmacophore rules definition, it can occur.
            atm_grp.features = features

    def get_shortest_path_length(self, src_grp, trgt_grp, cutoff=None):
        """Compute the shortest path length between two atom groups ``src_grp`` and ``trgt_grp``.

        The shortest path between two atom groups is defined as the shortest path between any of their atoms,
        which are calculated using Dijkstra’s algorithm and the graph ``graph``.

        If there is not any path between ``src_grp`` and ``trgt_grp``, infinite is returned.

        Parameters
        ----------
        src_grp, trgt_grp : `AtomGroup`
            Two atom groups to calculate the shortest path.
        cutoff : int
            Only paths of length <= ``cutoff`` are returned. If None, all path lengths are considered.

        Returns
        -------
         : int or float('inf'):
            The shortest path.
        """
        shortest_path_size = float('inf')
        for src_atm in src_grp.atoms:
            for trgt_atm in trgt_grp.atoms:
                try:
                    dist, path = single_source_dijkstra(self.graph, src_atm, trgt_atm, cutoff=cutoff)
                    if dist < shortest_path_size:
                        shortest_path_size = dist
                except Exception:
                    pass
        return shortest_path_size

    def save(self, output_file, compressed=True):
        """Write the pickled representation of the `AtomGroupsManager` object to the file ``output_file``.

        Parameters
        ----------
        output_file : str
            The output file.
        compressed : bool, optional
            If True (the default), compress the pickled representation as a gzip file (.gz).

        Raises
        -------
        FileNotCreated
            If the file could not be created.
        """
        pickle_data(self, output_file, compressed)

    @staticmethod
    def load(input_file):
        """Load the pickled representation of an `AtomGroupsManager` object saved at the file ``input_file``.

        Returns
        ----------
         : `AtomGroupsManager`
            The reconstituted `AtomGroupsManager` object, including its set of atom groups and interactions.

        Raises
        -------
        PKLNotReadError
            If the file could not be loaded.
        """
        return unpickle_data(input_file)

    def __len__(self):
        # Number of atom groups.
        return self.size

    def __iter__(self):
        """Iterate over children."""
        for atm_grp in self.atm_grps:
            yield atm_grp


class AtomGroup():
    """ Represent single atoms, chemical functional groups, or simply an arrangement of atoms as in hydrophobes.

    Parameters
    ----------
    atoms : iterable of :class:`~luna.mol.atom.ExtendedAtom`
        A sequence of atoms.
    features : iterable of :class:`~luna.mol.features.ChemicalFeature`, optional
        A sequence of chemical features.
    interactions : iterable of :class:`~luna.interaction.type.InteractionType`, optional
        A sequence of interactions established by an atom group.
    recursive : bool
        If True, add the new atom group to the list of atom groups of each atom in ``atoms``.
    manager : `AtomGroupsManager`, optional
        The `AtomGroupsManager` object that contains this `AtomGroup` object.
    """

    def __init__(self, atoms, features=None, interactions=None, recursive=True, manager=None):
        self._atoms = sorted(atoms)

        # Atom properties
        self._coords = im.atom_coordinates(atoms)
        self._centroid = im.centroid(self.coords)
        self._normal = None

        features = features or []
        self._features = sorted(features)

        self._interactions = interactions or []
        self._hash_cache = None

        self._manager = manager

        self._recursive = recursive

        if recursive:
            for atm in self.atoms:
                atm.add_atm_grps([self])

    @property
    def atoms(self):
        """iterable of :class:`~luna.mol.atom.ExtendedAtom`, read-only: The sequence of atoms that belong to an atom group."""
        return self._atoms

    @property
    def compounds(self):
        """set of :class:`~luna.MyBio.PDB.Residue.Residue`, read-only: The set of unique compounds that contain the atoms in ``atoms``.

        As an atom group can be formed by the union of two or more compounds (e.g., amide of peptide bonds), it may return
        more than one compound.
        """
        return set([a.parent for a in self._atoms])

    @property
    def coords(self):
        """ array-like of floats : Atomic coordinates (x, y, z) of each atom in ``atoms``."""
        return self._coords

    @property
    def centroid(self):
        """ array-like of floats, read-only: The centroid (x, y, z) of the atom group.

        If ``atoms`` contains only one atom, then ``centroid`` returns the same as ``coords``.
        """
        return self._centroid

    @property
    def normal(self):
        """array-like of floats, read-only: The normal vector (x, y, z) of the points given by ``coords``."""
        if self._normal is None:
            self._normal = im.calc_normal(self.coords)
        return self._normal

    @property
    def features(self):
        """iterable of :class:`~luna.mol.features.ChemicalFeature`: A sequence of chemical features.

        To add or remove a feature use :py:meth:`add_features`
        or :py:meth:`remove_features`, respectively."""
        return self._features

    @features.setter
    def features(self, features):
        self._features = sorted(features)
        # Reset hash.
        self._hash_cache = None

    @property
    def feature_names(self):
        """iterable of str: The name of each chemical feature in ``features``."""
        return [f.name for f in self.features]

    @property
    def interactions(self):
        """iterable of :class:`~luna.interaction.type.InteractionType`: The sequence of interactions
        established by an atom group.

        To add or remove an interaction use :py:meth:`add_interactions`
        or :py:meth:`remove_interactions`, respectively."""
        return self._interactions

    @interactions.setter
    def interactions(self, interactions):
        self._interactions = interactions

    @property
    def manager(self):
        """`AtomGroupsManager`: The `AtomGroupsManager` object that contains an `AtomGroup` object."""
        return self._manager

    @manager.setter
    def manager(self, manager):
        if isinstance(manager, AtomGroupsManager):
            self._manager = manager
        else:
            raise IllegalArgumentError("The informed atom group manager must be an instance of '%s'." % AtomGroupsManager)

    @property
    def size(self):
        """int: The number of atoms comprising an atom group."""
        return len(self.atoms)

    def has_atom(self, atom):
        """Check if an atom group contains a given atom ``atom``.

        Parameters
        ----------
        atom : :class:`~luna.mol.atom.ExtendedAtom`

        Returns
        -------
         : bool
            If the atom group contains or not ``atom``.
        """
        return atom in self.atoms

    def contain_group(self, atm_grp):
        """Check if the atom group ``atm_grp`` is a subset of this atom group.

        For example, consider the benzene molecule.
        Its aromatic ring itself forms an `AtomGroup` object composed of all of its six atoms.
        Consider now any subset of carbons in the benzene molecule.
        This subset forms an `AtomGroup` object that is part of the group formed by the aromatic ring.
        Therefore, in this example, :meth:`contain_group` will return True because the aromatic ring
        contains the subset of hydrophobic atoms.

        Parameters
        ----------
        atm_grp : :class:`~luna.mol.groups.AtomGroup`

        Returns
        -------
         : bool
            If one atom group contains another atom group.
        """
        return set(atm_grp.atoms).issubset(set(self.atoms))

    def get_serial_numbers(self):
        """Get the serial number of each atom in an atom group."""
        return [a.get_serial_number() for a in self.atoms]

    def get_chains(self):
        """Get all unique chains in an atom group."""
        return sorted(set([a.get_parent_by_level("C").id for a in self.atoms]))

    def get_interactions_with(self, atm_grp):
        """Get all interactions that an atom group establishes with another atom group ``atm_grp``.

        Returns
        -------
         : iterable of :class:`~luna.interactions.type.InteractionType`
           All interactions established with the atom group ``atm_grp``.
        """
        target_interactions = []

        for inter in self.interactions:
            if inter.src_grp == atm_grp or inter.trgt_grp == atm_grp:
                target_interactions.append(inter)
        return target_interactions

    def get_shortest_path_length(self, trgt_grp, cutoff=None):
        """Compute the shortest path length between this atom group to another atom group ``trgt_grp``.

        The shortest path between two atom groups is defined as the shortest path between any of their atoms,
        which are calculated using Dijkstra’s algorithm.

        If ``manager`` is not provided, None is returned.

        If there is not any path between ``src_grp`` and ``trgt_grp``, infinite is returned.

        Parameters
        ----------
        trgt_grp : `AtomGroup`
            The target atom group to calculate the shortest path.
        cutoff : int, optional
            Only paths of length <= ``cutoff`` are returned. If None, all path lengths are considered.

        Returns
        -------
         : int, float('inf'), or None:
            The shortest path.
        """
        if self.manager is not None:
            return self.manager.get_shortest_path_length(self, trgt_grp, cutoff)
        return None

    def add_features(self, features):
        """ Add :class:`~luna.mol.features.ChemicalFeature` objects to ``features``."""
        self._features = sorted(set(self.features + list(features)))
        # Reset hash.
        self._hash_cache = None

    def remove_features(self, features):
        """ Remove :class:`~luna.mol.features.ChemicalFeature` objects from ``features``."""
        self._features = sorted(set(self.features) - set(features))
        # Reset hash.
        self._hash_cache = None

    def add_interactions(self, interactions):
        """ Add :class:`~luna.interaction.type.InteractionType` objects to ``interactions``."""
        self._interactions = list(set(self.interactions + list(interactions)))

    def remove_interactions(self, interactions):
        """ Remove :class:`~luna.interaction.type.InteractionType` objects from ``interactions``."""
        self._interactions = list(set(self.interactions) - set(interactions))

    def is_water(self):
        """Return True if all atoms in the atom group belong to water molecules."""
        return all([a.parent.is_water() for a in self.atoms])

    def is_hetatm(self):
        """Return True if all atoms in the atom group belong to hetero group, i.e., non-standard residues of proteins,
        DNAs, or RNAs, as well as atoms in other kinds of groups, such as carbohydrates,
        substrates, ligands, solvent, and metal ions.

        Hetero groups are designated by the flag HETATM in the PDB format."""
        return all([a.parent.is_hetatm() for a in self.atoms])

    def is_residue(self):
        """Return True if all atoms in the atom group belong to standard residues of proteins."""
        return all([a.parent.is_residue() for a in self.atoms])

    def is_nucleotide(self):
        """Return True if all atoms in the atom group belong to nucleotides."""
        return all([a.parent.is_nucleotide() for a in self.atoms])

    def is_mixed(self):
        """Return True if the atoms in the atom group belong to different compound classes
        (water, hetero group, residue, or nucleotide)."""
        return len(set([a.parent.get_class() for a in self.atoms])) > 1

    def has_water(self):
        """Return True if at least one atom in the atom group belongs to a water molecule."""
        return any([a.parent.is_water() for a in self.atoms])

    def has_hetatm(self):
        """Return True if at least one atom in the atom group belongs to a hetero group, i.e., non-standard residues of proteins,
        DNAs, or RNAs, as well as atoms in other kinds of groups, such as carbohydrates,
        substrates, ligands, solvent, and metal ions."""
        return any([a.parent.is_hetatm() for a in self.atoms])

    def has_residue(self):
        """Return True if at least one atom in the atom group belongs to a standard residue of proteins."""
        return any([a.parent.is_residue() for a in self.atoms])

    def has_nucleotide(self):
        """Return True if at least one atom in the atom group belongs to a nucleotide."""
        return any([a.parent.is_nucleotide() for a in self.atoms])

    def has_target(self):
        """Return True if at least one compound is the target of LUNA's analysis"""
        return any([a.parent.is_target() for a in self.atoms])

    def as_json(self):
        """Represent the atom group as a dict containing the atoms, compounds,
        features, and compound classes (water, hetero group, residue, or nucleotide).

        The dict is defined as follows:

            * ``atoms`` (iterable of :class:`~luna.mol.atom.ExtendedAtom`): the list of atoms comprising the atom group;
            * ``compounds`` (iterable of :class:`~luna.MyBio.PDB.Residue.Residue`): the list of unique compounds that contain the atoms;
            * ``classes`` (iterable of str): the list of compound classes;
            * ``features`` (iterable of :class:`~luna.mol.features.ChemicalFeature`): the atom group's list of chemical features.
        """
        grp_obj = {}
        grp_obj["atoms"] = [atm.as_json() for atm in self.atoms]
        grp_obj["compounds"] = [comp.as_json() for comp in self.compounds]
        grp_obj["features"] = [feat.name for feat in self.features if feat.name != "Atom"]
        grp_obj["classes"] = [comps.get_class() for comps in self.compounds]

        return grp_obj

    def clear_refs(self):
        """References to this `AtomGroup` instance will be removed from the list of atom groups of
        each atom in ``atoms``."""
        if self._recursive:
            for atm in self.atoms:
                atm.remove_atm_grps([self])

    def __repr__(self):
        return '<AtomGroup: [%s]>' % ', '.join([str(x) for x in self.atoms])

    def __eq__(self, other):
        """Overrides the default implementation"""
        if type(self) == type(other):
            return self.atoms == other.atoms and self.features == other.features
        return False

    def __ne__(self, other):
        """Overrides the default implementation"""
        return not self.__eq__(other)

    def __lt__(self, other):
        atms1 = tuple(sorted(self.atoms))
        atms2 = tuple(sorted(other.atoms))
        return atms1 < atms2

    def __len__(self):
        # Number of atoms.
        return self.size

    def __hash__(self):
        """Overrides the default implementation"""
        if self._hash_cache is None:
            # Transform atoms and features list into an imutable data structure.
            # The lists are sorted in order to avoid dependence on appending order.
            atoms_tuple = tuple(self.atoms)
            feat_tuple = tuple(self.features)
            self._hash_cache = hash((atoms_tuple, feat_tuple, self.__class__))
        return self._hash_cache


class PseudoAtomGroup(AtomGroup):
    """Represent only the atoms from an `AtomGroup` object that are involved in an interaction.

    Currently, this class is only used during the generation of LUNA's fingerprints.

    Parameters
    ----------
    parent_grp : `AtomGroup`
        The atom group that contains the subset of atoms ``atoms``.
    atoms : iterable of :class:`~luna.mol.atom.ExtendedAtom`
        A sequence of atoms.
    features : iterable of :class:`~luna.mol.features.ChemicalFeature`, optional
        A sequence of chemical features.
    interactions : iterable of :class:`~luna.interaction.type.InteractionType`, optional
        A sequence of interactions established by an atom group.
    """

    def __init__(self, parent_grp, atoms, features=None, interactions=None):
        self.parent_grp = parent_grp

        super().__init__(atoms, features, interactions, recursive=False)

    def __repr__(self):
        return '<PseudoAtomGroup: [%s]>' % ', '.join([str(x) for x in self.atoms])


class AtomGroupPerceiver():

    """Perceive and create atom groups for molecules.

    Parameters
    ----------
    feature_extractor : :class:`~luna.mol.features.FeatureExtractor`
        Perceive pharmacophoric properties from molecules.
    add_h : bool
        If True, add hydrogen to the molecules.
    ph : float, optional
        If not None, add hydrogens appropriate for pH ``ph``.
    amend_mol : bool
        If True, apply validation and standardization of molecules read from a PDB file.
    charge_model : class:`~luna.mol.charge_model.ChargeModel`
        A charge model object. By default, the implementation of OpenEye charge model is used.
    mol_obj_type : {"rdkit", "openbabel"}
        If "rdkit", parse the converted molecule with RDKit and return an instance of :class:`rdkit.Chem.rdchem.Mol`.
        If "openbabel", parse the converted molecule with Open Babel and return an instance of
        :class:`openbabel.pybel.Molecule`.
    expand_selection : bool
        If True (the default), perceive features for a given molecule considering all nearby molecules.
        The goal is to identify any covalently bonded molecules that may alter the
        pharmacophoric properties or chemical functional groups.

        For instance, consider an amide of a peptide bond.
        If ``expand_selection`` is False, the residues forming the peptide bond will be analyzed separately,
        which will make the oxygen and the nitrogen of the amide to be perceived as carbonyl oxygen
        and amine, respectively.
        On the other hand, if ``expand_selection`` is True, the covalent bond between the residues will be
        identified and the amide will be correctly perceived.
    radius : float
        If ``expand_selection`` is True, select all molecules up to a maximum of ``radius`` away (measured in Å).
        The default value is 2.2, which comprises covalent bond distances.
    tmp_path : str, optional
        A temporary directory to where temporary files will be saved.
        If not provided, the system's default temporary directory will be used instead.
    critical : bool
        If False, ignore any errors during the processing a molecule and continue to the next one.
        The default value is True, which implies that any errors will raise an exception.

    Raises
    ------
    IllegalArgumentError
        If ``mol_obj_type`` is not either 'rdkit' nor 'openbabel'.
    """

    def __init__(self, feature_extractor, add_h=False, ph=None, amend_mol=True,
                 charge_model=OpenEyeModel(), mol_obj_type="rdkit",
                 expand_selection=True, radius=COV_SEARCH_RADIUS,
                 cache=None, tmp_path=None, critical=True):

        if mol_obj_type not in ACCEPTED_MOL_OBJ_TYPES:
            acc_mol_obj_types = ", ".join(ACCEPTED_MOL_OBJ_TYPES)
            raise IllegalArgumentError("Objects of type '%s' are not currently"
                                       "accepted. The available options are: "
                                       "%s." % (mol_obj_type,
                                                acc_mol_obj_types))

        self.feature_extractor = feature_extractor

        self.add_h = add_h
        self.ph = ph
        # If the user decided not to add hydrogens,
        # it will try to use the existing ones.
        self.keep_hydrog = not self.add_h

        self.amend_mol = amend_mol
        self.charge_model = charge_model
        self.mol_obj_type = mol_obj_type
        self.expand_selection = expand_selection
        self.radius = radius

        self.cache = cache
        self.tmp_path = tmp_path

        # If the pharmacophoric perception is critical, any exception during
        # the processing of a molecule will stop the whole processing.
        self.critical = critical

    def perceive_atom_groups(self, compounds, mol_objs_dict=None):
        """Perceive and create atom groups for each molecule in ``compounds``.

        Parameters
        ----------
        compounds : iterable of :class:`~luna.MyBio.PDB.Residue.Residue`
            A sequence of molecules.
        mol_objs_dict : dict
            Map a compound, represented by its id, to a molecular object (:class:`~luna.wrappers.base.MolWrapper`, :class:`rdkit.Chem.rdchem.Mol`,
            or :class:`openbabel.pybel.Molecule`).

            This parameter can be used in cases where the ligand is read from a molecular file
            and no standardization or validation is required.

        Returns
        -------
         : `AtomGroupsManager`
            An `AtomGroupsManager` object containing all atom groups perceived for the molecules in ``compounds``.
        """

        mol_objs_dict = mol_objs_dict or {}

        self.atm_grps_mngr = AtomGroupsManager()
        self.atm_mapping = {}

        compounds = set(compounds)
        init_comps_set = set(compounds)

        if self.cache:
            cached_compounds = set([c for c in compounds
                                    if self.cache.is_compound_cached(c)])
            compounds = compounds - cached_compounds

        # Controls which compounds are allowed to expand
        # when 'expand_selection' is set ON.
        pruned_comps = set()
        # Create a queue of compounds to be processed.
        comp_queue = set(compounds)

        while comp_queue:
            try:
                comp = comp_queue.pop()

                logger.debug("It will try to perceive the features of the "
                             "compound '%s' as no predefined properties "
                             "were provided." % comp.full_name)

                mol_obj = mol_objs_dict.get(comp.id, None)

                # If no OBMol object was defined and expand_selection
                # was set ON and it is not a border compound, it will
                # get all compounds around the target and create a new
                # OBMol object with them.
                if mol_obj is None and self.expand_selection:
                    # Remove the ligand from the list when
                    # it was provided as an OBMol object.
                    comp_list = [c for c in get_proximal_compounds(comp)
                                 if c.id not in mol_objs_dict]

                    # Expands the queue of compounds
                    # with any new border compound.
                    if comp not in pruned_comps:
                        border_comps = set(comp_list) - compounds
                        pruned_comps |= border_comps
                        comp_queue |= border_comps

                # Otherwise, the compound list will be composed
                # only by the target compound.
                else:
                    comp_list = [comp]

                # Recover all atoms from the previous selected compounds.
                selector = Selector(keep_altloc=False,
                                    keep_hydrog=self.keep_hydrog)
                atoms = tuple([a for r in comp_list
                               for a in r.get_unpacked_list()
                               if selector.accept_atom(a)])

                self._assign_properties(comp, atoms, mol_obj)
            except Exception:
                logger.debug("Features for the compound '%s' were "
                             "not correctly perceived." % comp)

                if self.critical:
                    raise

        if self.cache and cached_compounds:
            self._apply_cache(cached_compounds, comp.get_parent_by_level("M"))

        # Remove atom groups not comprising the provided
        # compound list (parameter 'compounds').
        remove_atm_grps = []
        for atm_grp in self.atm_grps_mngr:
            if any([c in init_comps_set for c in atm_grp.compounds]) is False:
                remove_atm_grps.append(atm_grp)
        self.atm_grps_mngr.remove_atm_grps(remove_atm_grps)

        return self.atm_grps_mngr

    def _apply_cache(self, compounds, model):

        def copy_atom(ref_atm, force_update=False):
            ref_res = ref_atm.parent
            ref_chain = ref_res.parent

            # Recover current structure entities.
            chain = model[ref_chain.id]
            res = chain[ref_res.id]
            atm = res[ref_atm.id]

            # If an ExtendedAtom already exists, return it.
            if not force_update and atm in self.atm_mapping:
                return self.atm_mapping[atm]

            # Define a new ExtendedAtom.
            ext_atm = self._new_extended_atom(atm)

            # Set its atomic invariants.
            ext_atm.invariants = ref_atm.invariants

            # Add neighbors' information to the new ExtendedAtom.
            for nbi in ref_atm.neighbors_info:
                atom_info = AtomData(nbi.atomic_num,
                                     nbi.coord, nbi.bond_type,
                                     nbi.full_id)
                ext_atm.add_nb_info([atom_info])

            return ext_atm

        comp_names = set([c.full_name for c in compounds])

        for atm_grp in sorted(self.cache.atm_grps_mngr):
            if any([c.full_name in comp_names for c in atm_grp.compounds]):
                atoms = []
                for atm in atm_grp.atoms:
                    # Copy extended atom, i.e., create a new extended
                    # atom or recover an existing one and set invariants
                    # and neighbors' information when necessary.
                    ext_atm = copy_atom(atm, force_update=True)
                    atoms.append(ext_atm)

                # Create the new atom group.
                features = atm_grp.features
                self.atm_grps_mngr.new_atm_grp(atoms, features)

        for edge in sorted(self.cache.atm_grps_mngr.graph.edges):
            # Skip already existing edges.
            if edge in self.atm_grps_mngr.graph.edges:
                continue

            atoms = []
            for atm in edge:
                # Copy extended atom, i.e., create a new
                # extended atom or recover an existing one
                # and set invariants and neighbors' information
                # if necessary.
                ext_atm = copy_atom(atm)
                atoms.append(ext_atm)

            # Update the AtomGroupsManager object with a new edge.
            self.atm_grps_mngr.graph.add_edge(atoms[0], atoms[1], weight=1)

    def _assign_properties(self, target_compound, target_atoms, mol_obj=None):

        # If no OBMol was defined, create a new one with the compound list.
        ignored_atoms = []
        if mol_obj is None:
            mol_obj, ignored_atoms = self._get_mol_from_entity(target_compound, target_atoms)
        # Create a new MolWrapper object.
        mol_obj = MolWrapper(mol_obj)

        # If the add_h property is set to False, the code will not remove any existing hydrogens from the PDB structure.
        # In these situations, the list of atoms may contain hydrogens. But, we do not need to attribute properties to hydrogens.
        # We just need them to correctly set properties to heavy atoms. So let's just ignore them.
        target_atoms = [atm for atm in target_atoms if atm.element != "H" and atm not in ignored_atoms]

        if mol_obj.get_num_heavy_atoms() != len(target_atoms):
            raise MoleculeSizeError("The number of heavy atoms in the PDB selection and in the MOL file are different.")

        # Ignore hydrogen atoms.
        atm_obj_list = [atm for atm in mol_obj.get_atoms() if atm.get_atomic_num() != 1]

        atm_map = {}
        trgt_atms = {}
        for i, atm_obj in enumerate(atm_obj_list):
            atm_map[atm_obj.get_idx()] = target_atoms[i].get_full_id()

            trgt_atms[target_atoms[i].get_full_id()] = self._new_extended_atom(target_atoms[i])

            # Invariants are still None, so it means a new ExtendedAtom has just been created.
            if trgt_atms[target_atoms[i].get_full_id()].invariants is None:
                # Update atomic invariants for new ExtendedAtoms created.
                trgt_atms[target_atoms[i].get_full_id()].invariants = atm_obj.get_atomic_invariants()

            # The current atom already has its invariants updated, which means the atom was created previously
            # by another target compound.
            else:
                # In this case, if the current atom belongs to the target compound, we need to prioritize the current target compound
                # and overwrite the atom's invariants. By doing so, we always guarantee the most accurate information of an atom.
                #
                # It is done because some atoms are updated by centroids (target compound) not corresponding to its real compound.
                # When it happens, the information of covalently bonded atoms may be lost due to the short COV_SEARCH_RADIUS used when
                # selecting nearby compounds. As a consequence, these atoms will have their properties and topology incorrectly
                # perceived.
                #
                # However, we need to highlight that the main goal of selecting atoms using COV_SEARCH_RADIUS is to find atoms
                # covalently bonded to the current compound. Therefore, atoms in the border of nearby molecules are not important
                # right now and they will be updated in their own time.
                if trgt_atms[target_atoms[i].get_full_id()].parent == target_compound:
                    # Update the atomic invariants for an ExtendedAtom.
                    trgt_atms[target_atoms[i].get_full_id()].invariants = atm_obj.get_atomic_invariants()

        # Set all neighbors, i.e., covalently bonded atoms.
        for bond_obj in mol_obj.get_bonds():
            bgn_atm_obj = bond_obj.get_begin_atom()
            end_atm_obj = bond_obj.get_end_atom()

            # At least one of the atoms must be a non-hydrogen atom.
            if bgn_atm_obj.get_atomic_num() != 1 or end_atm_obj.get_atomic_num() != 1:
                # If the atom 1 is not a hydrogen, add atom 2 to its neighbor list.
                if bgn_atm_obj.get_atomic_num() != 1:
                    # If the current bgn_atm_obj consists of an atom from the current target compound, we can update its information.
                    # Other compounds are ignored for now as they will have their own time to update its information.
                    if trgt_atms[atm_map[bgn_atm_obj.get_idx()]].parent == target_compound:
                        full_id = atm_map.get(end_atm_obj.get_idx())
                        coord = mol_obj.get_atom_coord_by_id(end_atm_obj.get_id())
                        atom_info = AtomData(end_atm_obj.get_atomic_num(), coord, bond_obj.get_bond_type(), full_id)
                        trgt_atms[atm_map[bgn_atm_obj.get_idx()]].add_nb_info([atom_info])

                # If the atom 2 is not a hydrogen, add atom 1 to its neighbor list.
                if end_atm_obj.get_atomic_num() != 1:
                    # If the current end_atm_obj consists of an atom from the current target compound, we can update its information.
                    # Other compounds are ignored for now as they will have their own time to update its information.
                    if trgt_atms[atm_map[end_atm_obj.get_idx()]].parent == target_compound:
                        full_id = atm_map.get(bgn_atm_obj.get_idx())
                        coord = mol_obj.get_atom_coord_by_id(bgn_atm_obj.get_id())
                        atom_info = AtomData(bgn_atm_obj.get_atomic_num(), coord, bond_obj.get_bond_type(), full_id)
                        trgt_atms[atm_map[end_atm_obj.get_idx()]].add_nb_info([atom_info])

        group_features = self.feature_extractor.get_features_by_groups(mol_obj, atm_map)
        for key in group_features:
            grp_obj = group_features[key]

            # It only accepts groups containing at least one atom from the target compound.
            if any([trgt_atms[i].parent == target_compound for i in grp_obj["atm_ids"]]) is False:
                continue

            atoms = [trgt_atms[i] for i in grp_obj["atm_ids"]]
            self.atm_grps_mngr.new_atm_grp(atoms, grp_obj["features"])

        # Update the graph in the AtomGroupsManager object with the current compound information.
        for mybio_atm in target_compound.get_atoms():
            if mybio_atm.get_full_id() in trgt_atms and mybio_atm.parent == target_compound:
                atm = trgt_atms[mybio_atm.get_full_id()]

                for nb_info in atm.neighbors_info:
                    if nb_info.full_id in trgt_atms:
                        self.atm_grps_mngr.graph.add_edge(atm, trgt_atms[nb_info.full_id], weight=1)
        return True

    def _new_extended_atom(self, atm, invariants=None):
        if atm not in self.atm_mapping:
            self.atm_mapping[atm] = ExtendedAtom(atm, invariants=invariants)

        return self.atm_mapping[atm]

    def _get_mol_from_entity(self, target_compound, target_atoms):
        validate_mol = self.amend_mol
        standardize_mol = self.amend_mol and target_compound.is_residue()

        atom_selector = AtomSelector(target_atoms, keep_altloc=False, keep_hydrog=self.keep_hydrog)

        return biopython_entity_to_mol(target_compound.get_parent_by_level('M'), atom_selector,
                                       validate_mol=validate_mol, standardize_mol=standardize_mol, add_h=self.add_h,
                                       ph=self.ph, mol_obj_type=self.mol_obj_type, tmp_path=self.tmp_path)


class AtomGroupNeighborhood:
    """ Class for fast neighbor atom groups searching.

    ``AtomGroupNeighborhood`` makes use of a KD Tree implemented in C, so it's fast.

    Parameters
    ----------
    atm_grps : iterable of `AtomGroup`, optional
        A sequence of `AtomGroup` objects, which is used in the queries.
        It can contain atom groups from different molecules.
    bucket_size : int
        Bucket size of KD tree.
        You can play around with this to optimize speed if you feel like it.
        The default value is 10.

    """

    def __init__(self, atm_grps, bucket_size=10):
        self.atm_grps = list(atm_grps)

        # get the coordinates
        coord_list = [ga.centroid for ga in self.atm_grps]

        # to Nx3 array of type float
        self.coords = np.array(coord_list).astype("f")
        assert(bucket_size > 1)
        assert(self.coords.shape[1] == 3)
        self.kdt = KDTree(3, bucket_size)
        self.kdt.set_coords(self.coords)

    def search(self, center, radius):
        """Return all atom groups in ``atm_grps`` that is up to a maximum of ``radius`` away (measured in Å) of ``center``.

        For atom groups with more than one atom, their centroid is used as a reference.
        """
        self.kdt.search(center, radius)
        indices = self.kdt.get_indices()
        n_grps_list = []
        atm_grps = self.atm_grps
        for i in indices:
            a = atm_grps[i]
            n_grps_list.append(a)

        return n_grps_list
