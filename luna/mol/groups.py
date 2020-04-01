from collections import defaultdict
from itertools import chain

import numpy as np

from Bio.KDTree import KDTree

# Open Babel
from openbabel import openbabel as ob
from openbabel.pybel import readfile
from openbabel.pybel import Molecule as PybelWrapper
# RDKit
from rdkit.Chem import MolFromMolBlock, MolFromMolFile, SanitizeFlags, SanitizeMol

from luna.MyBio.selector import ResidueSelector, AtomSelector
from luna.MyBio.util import save_to_file
from luna.graph.bellman_ford import bellman_ford
from luna.mol.interaction.contact import get_contacts_for_entity, get_cov_contacts_for_entity
from luna.mol.atom import ExtendedAtom, AtomData
from luna.mol.charge_model import OpenEyeModel
from luna.mol.validator import MolValidator
from luna.mol.standardiser import ResiduesStandardiser
from luna.mol.wrappers.obabel import convert_molecule
from luna.mol.wrappers.base import MolWrapper
from luna.mol.features import ChemicalFeature
from luna.util.file import get_unique_filename, remove_files
from luna.util.exceptions import MoleculeSizeError, IllegalArgumentError, MoleculeObjectError
from luna.util.default_values import ACCEPTED_MOL_OBJ_TYPES, COV_SEARCH_RADIUS, OPENBABEL
from luna.mol.interaction.type import InteractionType
from luna.mol.interaction import math as im
from luna.util.file import pickle_data, unpickle_data

import logging
logger = logging.getLogger()


SS_BOND_FEATURES = ["Atom", "Acceptor", "ChalcogenDonor", "Hydrophobic"]


class AtomGroupsManager():

    def __init__(self, atm_grps=None):
        self._atm_grps = []
        self._child_dict = {}
        self.add_atm_grps(atm_grps)

    @property
    def atm_grps(self):
        return self._atm_grps

    @property
    def child_dict(self):
        return self._child_dict

    @property
    def size(self):
        return len(self._atm_grps)

    @property
    def summary(self):
        summary = defaultdict(int)
        for grp in self.atm_grps:
            for feature in grp.features:
                summary[feature] += 1

        return summary

    def find_atm_grp(self, atoms):
        return self.child_dict.get(tuple(sorted(atoms)), None)

    def get_all_interactions(self):
        return set(chain.from_iterable([atm_grp.interactions for atm_grp in self.atm_grps]))

    def filter_by_types(self, types):
        for atm_grp in self.atm_grps:
            if set(types).issubset(set(atm_grp.feature_names)):
                yield atm_grp

    def add_atm_grps(self, atm_grps):
        atm_grps = atm_grps or []

        self._atm_grps = list(set(self._atm_grps + list(atm_grps)))

        for atm_grp in atm_grps:
            self.child_dict[tuple(sorted(atm_grp.atoms))] = atm_grp

    def remove_atm_grps(self, atm_grps):
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

        # Only hydrophobic atom groups.
        hydrop_atm_grps = list(self.filter_by_types(["Hydrophobic"]))

        # Hydrophobic islands dictionary. Keys are integer values and items are defined by a set of atom groups.
        hydrop_islands = defaultdict(set)

        # It stores a mapping of an atom (represented by its serial number) and a hydrophobic island (defined by its keys).
        atm_mapping = {}

        island_id = 0
        for atm_grp in hydrop_atm_grps:
            # Hydrophobic atoms are defined always as only one atom.
            atm = atm_grp.atoms[0]

            # Recover the groups of all neighbors of this atom (it will merge all existing islands).
            nb_grps = set([atm_mapping[nb] for nb in atm.neighborhood[atm.serial_number].keys() if nb in atm_mapping])

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
                atm_mapping[atm.serial_number] = island_id
            else:
                atm_mapping[atm.serial_number] = island_id
                hydrop_islands[island_id].add(atm)
            island_id += 1

        # Create AtomGroup objects for the hydrophobe islands
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
            key = tuple(sorted([atm_mapping[src_atm.serial_number], atm_mapping[trgt_atm.serial_number]]))

            island_island_inter[key].add(inter)

        interactions = set()
        for k in island_island_inter:
            island_atms = defaultdict(set)
            for inter in island_island_inter[k]:
                src_atm = inter.src_grp.atoms[0]
                trgt_atm = inter.trgt_grp.atoms[0]

                island_atms[atm_mapping[src_atm.serial_number]].add(src_atm)
                island_atms[atm_mapping[trgt_atm.serial_number]].add(trgt_atm)

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

    def save(self, output_file, compressed=True):
        pickle_data(self, output_file, compressed)

    @staticmethod
    def load(input_file):
        return unpickle_data(input_file)

    def __len__(self):
        # Number of atom groups.
        return self.size

    def __iter__(self):
        """Iterate over children."""
        for atm_grp in self.atm_grps:
            yield atm_grp


class AtomGroup():

    def __init__(self, atoms, features=None, interactions=None, recursive=True):
        self._atoms = sorted(atoms)

        # Atom properties
        self._coords = im.atom_coordinates(atoms)
        self._centroid = im.centroid(self.coords)
        self._normal = None

        features = features or []
        self._features = sorted(features)

        self._interactions = interactions or []
        self._hash_cache = None

        self._recursive = recursive

        if recursive:
            for atm in self.atoms:
                atm.add_atm_grps([self])

    @property
    def atoms(self):
        return self._atoms

    @property
    def compounds(self):
        return set([a.parent for a in self._atoms])

    @property
    def coords(self):
        return self._coords

    @property
    def centroid(self):
        return self._centroid

    @property
    def normal(self):
        if self._normal is None:
            self._normal = im.calc_normal(self.coords)
        return self._normal

    @property
    def features(self):
        return self._features

    @features.setter
    def features(self, features):
        self._features = sorted(features)
        # Reset hash.
        self._hash_cache = None

    @property
    def feature_names(self):
        return [f.name for f in self.features]

    @property
    def interactions(self):
        return self._interactions

    @interactions.setter
    def interactions(self, interactions):
        self._interactions = interactions

    @property
    def size(self):
        return len(self.atoms)

    def has_atom(self, atom):
        return atom in self.atoms

    def contain_group(self, atm_grp):
        return set(atm_grp.atoms).issubset(set(self.atoms))

    def get_serial_numbers(self):
        return [a.get_serial_number() for a in self.atoms]

    def get_chains(self):
        return sorted(set([a.get_parent_by_level("C").id for a in self.atoms]))

    def get_interactions_with(self, atm_grp):
        target_interactions = []

        for inter in self.interactions:
            if inter.src_grp == atm_grp or inter.trgt_grp == atm_grp:
                target_interactions.append(inter)

        return target_interactions

    def get_shortest_path_size(self, trgt_grp, merge_neighborhoods=False):
        if merge_neighborhoods:
            nb_graph = {**self.atoms[0].neighborhood, **trgt_grp.atoms[0].neighborhood}
        else:
            nb_graph = self.atoms[0].neighborhood

        shortest_path_size = float('inf')
        for src_atm in self.atoms:
            # d stores the path size from the source to the target, and p stores the predecessors from each target.
            d, p = bellman_ford(nb_graph, src_atm.serial_number)
            for trgt_atm in trgt_grp.atoms:
                if trgt_atm.serial_number in d:
                    if d[trgt_atm.serial_number] < shortest_path_size:
                        shortest_path_size = d[trgt_atm.serial_number]
        return shortest_path_size

    def add_features(self, features):
        self._features = sorted(set(self.features + list(features)))
        # Reset hash.
        self._hash_cache = None

    def remove_features(self, features):
        self._features = sorted(set(self.features) - set(features))
        # Reset hash.
        self._hash_cache = None

    def add_interactions(self, interactions):
        self._interactions = list(set(self.interactions + list(interactions)))

    def remove_interactions(self, interactions):
        self._interactions = list(set(self.interactions) - set(interactions))

    def is_water(self):
        """Return 1 if all compounds are water molecules."""
        return all([a.parent.is_water() for a in self.atoms])

    def is_hetatm(self):
        """Return 1 if all compounds are hetero groups."""
        return all([a.parent.is_hetatm() for a in self.atoms])

    def is_residue(self):
        """Return 1 if all compounds are amino acids."""
        return all([a.parent.is_residue() for a in self.atoms])

    def is_nucleotide(self):
        """Return 1 if all compounds are nucleotides."""
        return all([a.parent.is_nucleotide() for a in self.atoms])

    def is_mixed(self):
        """Return 1 if the compounds are from different classes."""
        return len(set([a.parent.get_class() for a in self.atoms])) > 1

    def has_water(self):
        """Return 1 if at least one compound is a water molecule."""
        return any([a.parent.is_water() for a in self.atoms])

    def has_hetatm(self):
        """Return 1 if at least one compound is an hetero groups."""
        return any([a.parent.is_hetatm() for a in self.atoms])

    def has_residue(self):
        """Return 1 if at least one compound is an amino acids."""
        return any([a.parent.is_residue() for a in self.atoms])

    def has_nucleotide(self):
        """Return 1 if at least one compound is a nucleotides."""
        return any([a.parent.is_nucleotide() for a in self.atoms])

    def has_target(self):
        """Return 1 if at least one compound is the target."""
        return any([a.parent.is_target() for a in self.atoms])

    def clear_refs(self):
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

    def __init__(self, parent_grp, atoms, features=None, interactions=None):
        self.parent_grp = parent_grp

        super().__init__(atoms, features, interactions, recursive=False)

    def __repr__(self):
        return '<PseudoAtomGroup: [%s]>' % ', '.join([str(x) for x in self.atoms])


class AtomGroupPerceiver():

    def __init__(self, feature_extractor, add_h=False, ph=None, amend_mol=True, mol_obj_type="rdkit",
                 charge_model=OpenEyeModel(), tmp_path=None, expand_selection=True, default_properties=None,
                 radius=COV_SEARCH_RADIUS):

        if mol_obj_type not in ACCEPTED_MOL_OBJ_TYPES:
            raise IllegalArgumentError("Objects of type '%s' are not currently accepted. "
                                       "The available options are: %s." % (mol_obj_type, ", ".join(ACCEPTED_MOL_OBJ_TYPES)))

        self.feature_extractor = feature_extractor

        self.add_h = add_h
        self.ph = ph
        # If the user decided not to add hydrogens it will try to use the existing ones.
        self.keep_hydrog = not self.add_h

        self.amend_mol = amend_mol
        self.mol_obj_type = mol_obj_type
        self.charge_model = charge_model
        self.tmp_path = tmp_path
        self.expand_selection = expand_selection
        self.default_properties = default_properties
        self.radius = radius

    def perceive_atom_groups(self, compounds, mol_objs_dict=None):
        mol_objs_dict = mol_objs_dict or {}

        self.atm_grps_mngr = AtomGroupsManager()

        self.atm_mapping = {}

        compounds = set(compounds)
        # Controls which compounds are allowed to expand when this option was set ON.
        pruned_comps = set()
        # Create a queue of compounds to be processed.
        comp_queue = set(compounds)

        while comp_queue:
            comp = comp_queue.pop()

            props_set = False

            if self.default_properties is not None:
                props_set = self._set_default_properties(comp)

            if props_set is False:
                logger.info("It will try to perceive the features of the compound %s as no predefined properties "
                            "was provided." % comp)

                mol_obj = mol_objs_dict.get(comp.id, None)

                # If no OBMol object was defined and expand_selection was set ON and it is not a border compound,
                # it will get all compounds around the target and create a new OBMol object with them.
                if mol_obj is None and self.expand_selection:
                    comp_list = self._get_proximal_compounds(comp)

                    # Expands the queue of compounds with any new border compound.
                    if comp not in pruned_comps:
                        border_comps = set(comp_list) - compounds
                        pruned_comps |= border_comps
                        comp_queue |= border_comps

                # Otherwise, the compound list will be composed only by the target compound.
                else:
                    comp_list = [comp]

                # Select compounds based on the compounds list defined previously.
                comp_sel = ResidueSelector(comp_list, keep_altloc=False, keep_hydrog=self.keep_hydrog)

                # Recover all atoms from the previous selected compounds.
                atoms = tuple([a for r in comp_list for a in r.get_unpacked_list() if comp_sel.accept_atom(a)])

                props_set = self._calculate_properties(comp, atoms, comp_sel, mol_obj)

            if not props_set:
                logger.warning("Features for the compound '%s' were not correctly perceived." % comp)

        # Remove atom groups not comprising the provided compound list (parameter 'compounds').
        remove_atm_grps = []
        for atm_grp in self.atm_grps_mngr:
            if any([c in compounds for c in atm_grp.compounds]) is False:
                remove_atm_grps.append(atm_grp)
        self.atm_grps_mngr.remove_atm_grps(remove_atm_grps)

        return self.atm_grps_mngr

    def _get_proximal_compounds(self, target_compound):
        model = target_compound.get_parent_by_level('M')
        proximal = get_contacts_for_entity(entity=model, source=target_compound, radius=self.radius, level='R')

        # Sorted by the compound order as in the PDB.
        return sorted(list(set([p[1] for p in proximal])), key=lambda r: (r.parent.parent.id, r.parent.id, r.idx))

    def _set_default_properties(self, target_compound):
        #
        # TODO: Limitations in this method, although not so impactful:
        #           - How to detect if the compound is the first one? In this case the N should be +.
        #           - What to do when there is a covalently bonded ligand? It may change some atom properties.
        #           - What to do with coordinated metals (similar to the previous problem)? It may change some atom properties.
        #           - What to do with hydrogens if the users asked to add it?
        #
        # Although this cases are not so important they may slightly change the results of the analysis.
        #
        # This function only works if the user informed default properties to the target compound.
        if target_compound.resname in self.default_properties:
            try:
                atm_sel = AtomSelector(keep_altloc=False, keep_hydrog=self.keep_hydrog)

                # If the add_h property is set to False, the code will not remove any existing hydrogens from the PDB structure.
                # In these situations, the list of atoms may contain hydrogens. But, we do not need to attribute properties to hydrogens.
                # We just need them to correctly set properties to heavy atoms. So let's just ignore them.
                atm_map = {atm.name: self._new_extended_atom(atm) for atm in target_compound.get_atoms() if (atm_sel.accept_atom(atm) and
                                                                                                             atm.element != "H")}
                atms_in_ss_bonds = set()

                # It stores compounds covalently bound to the informed compound.
                neighbors = set()

                if target_compound.is_water() is False:
                    model = target_compound.get_parent_by_level('M')
                    cov_atoms = get_cov_contacts_for_entity(entity=model, source=target_compound)

                    for atm1, atm2 in cov_atoms:
                        # Validate the atoms using the AtomSelector created previously.
                        if atm_sel.accept_atom(atm1) and atm_sel.accept_atom(atm2):
                            # Add all compounds covalently bonded to the target compound into the neighbors list.
                            if atm1.parent != target_compound:
                                neighbors.add(atm1.parent)
                            if atm2.parent != target_compound:
                                neighbors.add(atm2.parent)

                            if atm1.parent != atm2.parent and atm1.parent.resname == "CYS" and atm2.parent.resname == "CYS":
                                atms_in_ss_bonds.add(atm_map[atm1.name])
                                atms_in_ss_bonds.add(atm_map[atm2.name])

                            # If the atom 1 belongs to the target and is not a hydrogen, add atom 2 to its neighbor's list.
                            if atm1.parent == target_compound and atm1.element != "H":
                                atom_info = AtomData(ob.GetAtomicNum(atm2.element), atm2.coord, atm2.serial_number)
                                atm_map[atm1.name].add_nb_info([atom_info])
                            # If the atom 2 belongs to the target and is not a hydrogen, add atom 1 to its neighbor's list.
                            if atm2.parent == target_compound and atm2.element != "H":
                                atom_info = AtomData(ob.GetAtomicNum(atm1.element), atm1.coord, atm1.serial_number)
                                atm_map[atm2.name].add_nb_info([atom_info])

                for atms_str in self.default_properties[target_compound.resname]:
                    atoms = []
                    missing_atoms = []
                    for atm_name in atms_str.split(","):
                        if atm_name in atm_map:
                            atoms.append(atm_map[atm_name])
                        else:
                            missing_atoms.append(atm_name)

                    if missing_atoms:
                        # If the only missing atom is the OXT.
                        if len(missing_atoms) == 1 and missing_atoms[0] == "OXT":
                            # If there is no successor compound, it may be an indication that there is a missing atom, the OXT in this case.
                            # But, sometimes not having the successor compound may be caused by missing compounds instead of missing atoms.
                            # As it is not so important, we will only print a warning.
                            if "next" not in neighbors:
                                logger.warning("The group composed by the atoms '%s' will be ignored because the OXT atom is missing "
                                               "in the compound %s. If the atom is not the C-terminal just ignore this message. "
                                               % (atms_str.replace(",", ", "), target_compound))
                        else:
                            logger.warning("The group composed by the atoms '%s' will be ignored because some of them were not found in "
                                           "the compound %s. The missing atoms are: %s." % (atms_str.replace(",", ", "), target_compound,
                                                                                            ", ".join(missing_atoms)))
                    else:
                        # It checks if a cysteine atom is establishing a disulfide bond. If it does, it will read a predefined set of
                        # properties (SS_BOND_FEATURES). In this case, the SG is hydrophobic and is not a donor anymore.
                        if target_compound.is_residue() and target_compound.resname == "CYS" and atoms[0] in atms_in_ss_bonds:
                            features = [ChemicalFeature(f) for f in SS_BOND_FEATURES]
                        else:
                            features = [ChemicalFeature(f) for f in self.default_properties[target_compound.resname][atms_str].split(",")]

                        # It considers all first compounds in the chain as the N-terminal and won't evaluate if there are missing compounds
                        # in the N-terminal. Since not all PDBs will have a 'REMARK 465 MISSING RESIDUES' field and, also, given that this
                        # field is not reliable enough, the best way to deal with the N-terminal would be by aligning the structure with
                        # the protein sequence. However, it would force one to provide the sequence and it would also increase the
                        # processing time only to identify if the compound is an N-terminal or not. In most of the cases, the first
                        # compound will be the N-terminal. Moreover, in general (maybe always), binding sites will not be at the ends
                        # of the protein. Therefore, I opted to always attribute a 'PositivelyIonizable' property to the N atom from the
                        # first compound in the chain.
                        if atms_str == "N" and target_compound.idx == 0:
                            features.append(ChemicalFeature("PositivelyIonizable"))

                        self.atm_grps_mngr.new_atm_grp(atoms, list(set(features)))

                # Check for amides only when the target is an amino acid.
                if target_compound.is_residue():
                    for nb in neighbors:
                        # Ignore compounds that are not amino acids.
                        if not nb.is_residue():
                            continue
                        # Recover valid atoms by applying an atom selection.
                        nb_atms = {atm.name: atm for atm in nb.get_atoms() if atm_sel.accept_atom(atm)}
                        if nb.idx < target_compound.idx:
                            # If this compound is neighbor of the target compound N, then they form an amide (N-C=O).
                            if "N" in atm_map and "O" in nb_atms and "C" in nb_atms and atm_map["N"].is_neighbor(nb_atms["C"]):
                                # Define an amide group and add it to the compound groups.
                                # In the atm_map, the atoms were already transformed into ExtendedAtoms().
                                amide = [atm_map["N"], self._new_extended_atom(nb_atms["C"]), self._new_extended_atom(nb_atms["O"])]

                                self.atm_grps_mngr.new_atm_grp(amide, [ChemicalFeature("Amide")])

                        elif nb.idx > target_compound.idx:
                            # If this compound is neighbor of the target compound C, then they form an amide (N-C=O).
                            if "N" in nb_atms and "O" in atm_map and "C" in atm_map and atm_map["C"].is_neighbor(nb_atms["N"]):
                                # Define an amide group and add it to the compound groups.
                                # In the atm_map, the atoms were already transformed into ExtendedAtoms().
                                amide = [self._new_extended_atom(nb_atms["N"]), atm_map["C"], atm_map["O"]]

                                self.atm_grps_mngr.new_atm_grp(amide, [ChemicalFeature("Amide")])

                # Define a new graph as a dict object: each key is a serial number and the values are dictionaries with the serial
                # number of each neighbor of an atom. The algorithm of Bellman Ford requires weights, so we set all weights to 1.
                nb_graph = defaultdict(dict)

                for atm in atm_map.values():
                    # Add the neighbors of this atom. The number '1' defined below represents the edge weight (mandatory).
                    for nb_info in atm.neighbors_info:
                        # Ignore hydrogens: we do not need to include them into the neighborhood graph.
                        if nb_info.atomic_num != 1:
                            nb_graph[atm.serial_number][nb_info.serial_number] = 1
                            nb_graph[nb_info.serial_number][atm.serial_number] = 1

                    # Set the graph dictionary to each ExtendedAtom object. Since, dictionaries are passed as references, we can change
                    # the variable nb_graph and the changes will be updated in each ExtendedAtom object automatically.
                    atm.set_neighborhood(nb_graph)
                return True
            except Exception as e:
                logger.exception(e)
                return False
        return False

    def _calculate_properties(self, target_compound, target_atoms, compound_selector, mol_obj=None):

        try:
            # If no OBMol was defined, create a new one with the compound list.
            ignored_atoms = []
            if mol_obj is None:
                mol_obj, ignored_atoms = self._get_mol_from_entity(target_compound.get_parent_by_level('M'),
                                                                   target_compound, target_atoms, compound_selector)
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
                atm_map[atm_obj.get_idx()] = target_atoms[i].serial_number

                trgt_atms[target_atoms[i].serial_number] = self._new_extended_atom(target_atoms[i])

                # Invariants are still None, so it means a new ExtendedAtom has just been created.
                if trgt_atms[target_atoms[i].serial_number].invariants is None:
                    # Update atomic invariants for new ExtendedAtoms created.
                    trgt_atms[target_atoms[i].serial_number].invariants = atm_obj.get_atomic_invariants()

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
                    if trgt_atms[target_atoms[i].serial_number].parent == target_compound:

                        # Update the atomic invariants for an ExtendedAtom.
                        trgt_atms[target_atoms[i].serial_number].invariants = atm_obj.get_atomic_invariants()

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
                            serial_number = atm_map.get(end_atm_obj.get_idx())
                            coord = mol_obj.get_atom_coord_by_id(end_atm_obj.get_id())
                            atom_info = AtomData(end_atm_obj.get_atomic_num(), coord, serial_number)
                            trgt_atms[atm_map[bgn_atm_obj.get_idx()]].add_nb_info([atom_info])

                    # If the atom 2 is not a hydrogen, add atom 1 to its neighbor list.
                    if end_atm_obj.get_atomic_num() != 1:
                        # If the current end_atm_obj consists of an atom from the current target compound, we can update its information.
                        # Other compounds are ignored for now as they will have their own time to update its information.
                        if trgt_atms[atm_map[end_atm_obj.get_idx()]].parent == target_compound:
                            serial_number = atm_map.get(bgn_atm_obj.get_idx())
                            coord = mol_obj.get_atom_coord_by_id(bgn_atm_obj.get_id())
                            atom_info = AtomData(bgn_atm_obj.get_atomic_num(), coord, serial_number)
                            trgt_atms[atm_map[end_atm_obj.get_idx()]].add_nb_info([atom_info])

            group_features = self.feature_extractor.get_features_by_groups(mol_obj, atm_map)
            for key in group_features:
                grp_obj = group_features[key]

                # It only accepts groups containing at least one atom from the target compound.
                if any([trgt_atms[i].parent == target_compound for i in grp_obj["atm_ids"]]) is False:
                    continue

                atoms = [trgt_atms[i] for i in grp_obj["atm_ids"]]

                self.atm_grps_mngr.new_atm_grp(atoms, grp_obj["features"])

            # Define a new graph as a dict object: each key is a serial number and the values are dictionaries with the serial
            # number of each neighbor of an atom. The algorithm of Bellman Ford requires weights, so we set all weights to 1.
            nb_graph = {atm.serial_number: None for atm in trgt_atms.values()}
            for atm in trgt_atms.values():
                # Add the neighbors of this atom. The number '1' defined below represents the edge weight.
                nb_graph[atm.serial_number] = {i.serial_number: 1 for i in atm.neighbors_info if i.serial_number in nb_graph}
                # Set the graph dictionary to each ExtendedAtom object. Since, dictionaries are passed as references, we can change
                # the variable nb_graph and the changes will be updated in each ExtendedAtom object automatically.
                atm.set_neighborhood(nb_graph)

            return True
        except Exception as e:
            logger.exception(e)
            return False

    def _new_extended_atom(self, atm, invariants=None):
        if atm not in self.atm_mapping:
            self.atm_mapping[atm] = ExtendedAtom(atm, invariants=invariants)

        return self.atm_mapping[atm]

    def _get_mol_from_entity(self, entity, target_compound, target_atoms, compound_selector):
        # First it saves the selection into a PDB file and then it converts the file to .mol.
        # I had to do it because the OpenBabel 2.4.1 had a problem with some molecules containing aromatic rings.
        # In such cases, the aromatic ring was being wrongly perceived and some atoms received more double bonds than it
        # was expected. The version 2.3.2 works better. Therefore, I defined this version manually (openbabel property).
        filename = get_unique_filename(self.tmp_path)
        pdb_file = '%s_pdb-file.pdb' % filename
        logger.info("Saving the PDB object as a PDB file named '%s'." % pdb_file)
        save_to_file(entity, pdb_file, compound_selector)

        mol_file = '%s_mol-file.mol' % filename
        ob_opt = {"error-level": 5}
        logger.info("Converting the PDB file '%s' to a Mol file named '%s' using Open Babel." % (pdb_file, mol_file))
        if self.add_h:
            logger.info("Hydrogens will be added to the molecule.")
            if self.ph is not None:
                ob_opt["p"] = self.ph
            else:
                ob_opt["h"] = ""
        convert_molecule(pdb_file, mol_file, opt=ob_opt, openbabel=OPENBABEL)

        # Currently, ignored atoms are only metals.
        ignored_atoms = []

        mol_obj = None
        if self.amend_mol:
            logger.info("A validation will be performed and it will try to fix some errors.")

            try:
                mol_obj = next(readfile("mol", mol_file))
            except Exception as e:
                logger.exception(e)
                raise MoleculeObjectError("An error occurred while parsing the file '%s' with Open Babel and the molecule "
                                          "object could not be created. Check the logs for more information." % mol_file)

            if target_compound.is_residue():
                # If the add_h property is set to False, the code will not remove any existing hydrogens from the PDB structure.
                # In these situations, the list of atoms may contain hydrogens. But, we do not need to attribute properties to hydrogens.
                # We just need them to correctly set properties to heavy atoms. So let's just ignore them.
                target_atoms = [atm for atm in target_atoms if atm.element != "H"]

                mol_obj = MolWrapper(mol_obj)
                if mol_obj.get_num_heavy_atoms() != len(target_atoms):
                    raise MoleculeSizeError("The number of heavy atoms in the PDB selection and in the MOL file are different.")

                # Ignore hydrogen atoms.
                atm_obj_list = [atm for atm in mol_obj.get_atoms() if atm.get_atomic_num() != 1]

                atom_pairs = []
                for i, atm_obj in enumerate(atm_obj_list):
                    if target_atoms[i].parent.is_residue():
                        atom_pairs.append((atm_obj, target_atoms[i]))

                rs = ResiduesStandardiser(break_metal_bonds=False)
                mol_obj.unwrap().DeleteHydrogens()
                rs.standardise(atom_pairs)

                for i, atm_obj in enumerate(atm_obj_list):
                    if atm_obj.get_id() in rs.removed_atoms:
                        ignored_atoms.append(target_atoms[i])

                # After standardizing residues, we need to recreate the Mol file, otherwise implicit hydrogens will not be included in the
                # MOL object and, therefore, their coordinates could not be accessed. If you try to generate coordinates directly from the
                # object, hydrogens will be incorrectly placed.
                mol_obj = PybelWrapper(mol_obj.unwrap())
                new_mol_file = '%s_tmp-mol-file.mol' % filename
                mol_obj.write("mol", new_mol_file, overwrite=True)
                # Overwrite mol_file by converting the new molecular file using the user specified parameters.
                # Note that right now it will add explicit hydrogens to the molecules according to the provided pH.
                convert_molecule(new_mol_file, mol_file, opt=ob_opt, openbabel=OPENBABEL)

                # Let's finally read the correct and standardized molecular file.
                try:
                    mol_obj = next(readfile("mol", mol_file))
                except Exception as e:
                    logger.exception(e)
                    raise MoleculeObjectError("An error occurred while parsing the file '%s' with Open Babel and the molecule "
                                              "object could not be created. Check the logs for more information." % mol_file)

                # Remove temporary files.
                remove_files([new_mol_file])

            mv = MolValidator()
            is_valid = mv.validate_mol(mol_obj)
            logger.info('Validation finished!!!')

            if not is_valid:
                logger.warning("The molecular file '%s' contain invalid atoms. Check the logs for more information." % mol_file)

            if self.mol_obj_type == "rdkit":
                try:
                    # The sanitization is set off. We will apply it in the next step.
                    mol_obj = MolFromMolBlock(mol_obj.write('mol'), sanitize=False, removeHs=False)
                    # Sanitize molecule is applied now, so we will be able to catch the exceptions raised by RDKit,
                    # otherwise it would not be possible.
                    SanitizeMol(mol_obj, SanitizeFlags.SANITIZE_ALL)
                except Exception as e:
                    logger.exception(e)
                    raise MoleculeObjectError("An error occurred while parsing the molecular block generated by Open Babel with "
                                              "RDKit. Check the logs for more information." % mol_file)
        else:
            try:
                # Create a new Mol object.
                if self.mol_obj_type == "openbabel":
                    mol_obj = next(readfile("mol", mol_file))
                else:
                    # The sanitization is set off. We will apply it in the next statement.
                    mol_obj = MolFromMolFile(mol_file, sanitize=False, removeHs=False)
                    # Sanitize molecule is applied now, so we will be able to catch the exceptions raised by RDKit,
                    # otherwise it would not be possible.
                    SanitizeMol(mol_obj, SanitizeFlags.SANITIZE_ALL)
            except Exception as e:
                logger.exception(e)
                tool = "Open Babel" if self.mol_obj_type == "openbabel" else "RDKit"
                raise MoleculeObjectError("An error occurred while parsing the file '%s' with %s and the molecule "
                                          "object could not be created. Check the logs for more information." % (mol_file, tool))

        # Remove temporary files.
        remove_files([pdb_file, mol_file])

        return mol_obj, ignored_atoms


class AtomGroupNeighborhood:

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
        self.kdt.search(center, radius)
        indices = self.kdt.get_indices()
        n_grps_list = []
        atm_grps = self.atm_grps
        for i in indices:
            a = atm_grps[i]
            n_grps_list.append(a)

        return n_grps_list