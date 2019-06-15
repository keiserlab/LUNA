from collections import defaultdict
from itertools import chain

from rdkit.Chem import MolFromMolBlock, MolFromMolFile, SanitizeFlags, SanitizeMol
from pybel import readfile
from openbabel import etab

from MyBio.selector import ResidueSelector, AtomSelector
from MyBio.util import save_to_file

from graph.bellman_ford import bellman_ford
from mol.interaction.contact import get_contacts_for_entity, get_cov_contacts_for_entity
from mol.neighborhood import NbAtom, NbAtomData
from mol.charge_model import OpenEyeModel
from mol.validator import MolValidator
from mol.wrappers.obabel import convert_molecule
from mol.wrappers.base import MolWrapper
from mol.features import ChemicalFeature
from util.file import get_unique_filename, remove_files
from util.exceptions import MoleculeSizeError, IllegalArgumentError, MoleculeObjectError
from util.default_values import ACCEPTED_MOL_OBJ_TYPES, COV_SEARCH_RADIUS, OPENBABEL

import mol.interaction.math as im

import logging
logger = logging.getLogger()


SS_BOND_FEATURES = ["Atom", "Acceptor", "ChalcogenDonor", "Hydrophobic"]


class CompoundGroups():

    def __init__(self, target_residue, atm_grps=None):
        self._residue = target_residue
        self._atm_grps = atm_grps or []

    @property
    def residue(self):
        return self._residue

    @property
    def atm_grps(self):
        return self._atm_grps

    @property
    def summary(self):
        summary = defaultdict(int)
        for grp in self.atm_grps:
            for feature in grp.features:
                summary[feature] += 1

        return summary

    def add_atm_grps(self, atm_grps):
        self._atm_grps = list(set(self._atm_grps + atm_grps))

    def remove_atm_grps(self, atm_grps):
        self._atm_grps = list(set(self._atm_grps) - set(atm_grps))

        # NbAtom objects keep a list of all AtomGroup objects to which they belong to. So, if we don't clear the references directly,
        # the AtomGroup objects will still exist in the NbAtom list even when they were already removed from an instance of CompoundGroups.
        for atm_grp in atm_grps:
            atm_grp.clear_refs()


class AtomGroup():

    def __init__(self, atoms, features, interactions=None, recursive=True):

        self._atoms = list(atoms)

        # Atom properties
        self._coords = im.atom_coordinates(atoms)
        self._centroid = im.centroid(self.coords)
        self._normal = None

        self.features = features
        self._interactions = interactions or []
        self._recursive = recursive
        self._hash_cache = None

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
    def feature_names(self):
        return [f.name for f in self.features]

    @property
    def interactions(self):
        return self._interactions

    def has_atom(self, atom):
        return atom in self.atoms

    def contain_group(self, atm_grp):
        return set(atm_grp.atoms).issubset(set(self.atoms))

    def get_serial_numbers(self):
        return [a.get_serial_number() for a in self.atoms]

    def get_interactions_with(self, atm_grp):
        target_interactions = []

        for i in self.interactions:
            if i.atm_grp1 == atm_grp or i.atm_grp2 == atm_grp:
                target_interactions.append(i)

        return target_interactions

    def get_shortest_path_size(self, trgt_grp):
        shortest_path_size = float('inf')
        for src_atm in self.atoms:
            # d stores the path size from the source to the target, and p stores the predecessors from each target.
            d, p = bellman_ford(src_atm.neighborhood, src_atm.serial_number)
            for trgt_atm in trgt_grp.atoms:
                if trgt_atm.serial_number in d:
                    if d[trgt_atm.serial_number] < shortest_path_size:
                        shortest_path_size = d[trgt_atm.serial_number]
        return shortest_path_size

    def add_interactions(self, interactions):
        self._interactions = list(set(self.interactions + interactions))

    def remove_interactions(self, interactions):
        self._interactions = list(set(self.interactions) - set(interactions))

    def is_water(self):
        """Return 1 if all compounds are water molecules."""
        return all([a.parent.is_water() for a in self.atoms])

    def is_hetatm(self):
        """Return 1 if all compounds are hetero groups."""
        return all([a.parent.is_hetatm() for a in self.atoms])

    def is_aminoacid(self):
        """Return 1 if all compounds are amino acids."""
        return all([a.parent.is_aminoacid() for a in self.atoms])

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

    def has_aminoacid(self):
        """Return 1 if at least one compound is an amino acids."""
        return any([a.parent.is_aminoacid() for a in self.atoms])

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
        if isinstance(self, other.__class__):
            return self.atoms == other.atoms and self.features == other.features
        return False

    def __ne__(self, other):
        """Overrides the default implementation"""
        return not self.__eq__(other)

    def __hash__(self):
        """Overrides the default implementation"""
        if self._hash_cache is None:
            # Transform atoms and features list into an imutable data structure.
            # The lists are sorted in order to avoid dependence on appending order.
            atoms_tuple = tuple(sorted(self.atoms, key=hash))
            feat_tuple = tuple(sorted(self.features, key=hash))

            self._hash_cache = hash((atoms_tuple, feat_tuple))
        return self._hash_cache


class CompoundGroupPerceiver():

    def __init__(self, feature_extractor, add_h=False, ph=None, amend_mol=True, mol_obj_type="rdkit",
                 charge_model=OpenEyeModel(), tmp_path=None, expand_selection=True, default_properties=None,
                 radius=COV_SEARCH_RADIUS, group_hydrophobes=True):

        if mol_obj_type not in ACCEPTED_MOL_OBJ_TYPES:
            raise IllegalArgumentError("Objects of type '%s' are not currently accepted. "
                                       "The available options are: %s." % (mol_obj_type, ", ".join(ACCEPTED_MOL_OBJ_TYPES)))

        self.feature_extractor = feature_extractor
        self.add_h = add_h
        self.ph = ph
        self.amend_mol = amend_mol
        self.mol_obj_type = mol_obj_type
        self.charge_model = charge_model
        self.tmp_path = tmp_path
        self.expand_selection = expand_selection
        self.default_properties = default_properties
        self.radius = radius
        self.group_hydrophobes = group_hydrophobes

    def perceive_compound_groups(self, target_residue, mol_obj=None):
        comp_grps = None
        if self.default_properties is not None:
            comp_grps = self._get_default_properties(target_residue)

        if comp_grps is None:
            logger.info("It will calculate the properties ")

            # If no OBMol object was defined and expand_selection was set ON, it will get all residues around the
            # target and create a new OBMol object with them.
            if mol_obj is None and self.expand_selection:
                res_list = self._get_all_around_residue(target_residue)
            # Otherwise, the residue list will be composed only by the target residue.
            else:
                res_list = [target_residue]

            # Sorted by residue order as in the PDB.
            res_list = sorted(res_list, key=lambda r: r.get_full_id())

            # Select residues based on the residue list defined previously.
            res_sel = ResidueSelector(res_list, keep_altloc=False, keep_hydrog=False)
            # Recover all atoms from the previous selected residues.
            atoms = tuple([a for r in res_list for a in r.get_unpacked_list() if res_sel.accept_atom(a)])

            comp_grps = self._get_calculated_properties(target_residue, atoms, res_sel, mol_obj)

        if self.group_hydrophobes:
            self._merge_hydrophobes(comp_grps)

        return comp_grps

    def _get_all_around_residue(self, target_residue):
        model = target_residue.get_parent_by_level('M')
        proximal = get_contacts_for_entity(entity=model, source=target_residue, radius=self.radius, level='R')

        return sorted(list(set([p[1] for p in proximal])), key=lambda r: r.id)

    def _get_default_properties(self, target_residue):
        #
        # TODO: Limitations in this method, although not so impactful:
        #           - How to detect if the residue is the first one? In this case the N should be +.
        #           - What to do when there is a covalently bonded ligand? It may change some atom properties.
        #           - What to do with coordinated metals (similar to the previous problem)? It may change some atom properties.
        #           - What to do with hydrogens if the users asked to add it?
        #
        # Although this cases are not so important they may slightly change the results of the analysis.
        #
        # This function only works if the user informed default properties to the target residue.
        if target_residue.resname in self.default_properties:

            atm_sel = AtomSelector(keep_altloc=False, keep_hydrog=False)
            atm_map = {atm.name: NbAtom(atm) for atm in target_residue.get_atoms() if atm_sel.accept_atom(atm)}

            atms_in_ss_bonds = set()

            if target_residue.is_water() is False:
                model = target_residue.get_parent_by_level('M')
                cov_atoms = get_cov_contacts_for_entity(entity=model, source=target_residue)

                neighbors = set()
                for atm1, atm2 in cov_atoms:
                    # Validate the atoms using the AtomSelector created previously.
                    if atm_sel.accept_atom(atm1) and atm_sel.accept_atom(atm2):
                        # Add all residues covalently bonded to the target residue into the neighbors list.
                        if atm1.parent != target_residue:
                            neighbors.add(atm1.parent)
                        if atm2.parent != target_residue:
                            neighbors.add(atm2.parent)

                        if atm1.parent != atm2.parent and atm1.parent.resname == "CYS" and atm2.parent.resname == "CYS":
                            atms_in_ss_bonds.add(atm_map[atm1.name])
                            atms_in_ss_bonds.add(atm_map[atm2.name])

                        # If the atom 1 belongs to the target and is not a hydrogen, add atom 2 to its neighbor list.
                        if atm1.parent == target_residue and atm1.element != "H":
                            atom_info = NbAtomData(etab.GetAtomicNum(atm2.element), atm2.coord, atm2.serial_number)
                            atm_map[atm1.name].add_nb_info([atom_info])
                        # If the atom 2 belongs to the target and is not a hydrogen, add atom 1 to its neighbor list.
                        if atm2.parent == target_residue and atm2.element != "H":
                            atom_info = NbAtomData(etab.GetAtomicNum(atm1.element), atm1.coord, atm1.serial_number)
                            atm_map[atm2.name].add_nb_info([atom_info])

            comp_grps = CompoundGroups(target_residue)

            for atms_str in self.default_properties[target_residue.resname]:
                nb_atoms = []
                missing_atoms = []
                for atm_name in atms_str.split(","):
                    if atm_name in atm_map:
                        nb_atoms.append(atm_map[atm_name])
                    else:
                        missing_atoms.append(atm_name)

                if missing_atoms:
                    # If the only missing atom is the OXT.
                    if len(missing_atoms) == 1 and missing_atoms[0] == "OXT":
                        # If there is no successor residue, it may be an indication that there is a missing atom, the OXT in this case.
                        # But, sometimes not having the successor residue may be caused by missing residues instead of missing atoms.
                        # As it is not so important, we will always print the alert.
                        if "next" not in neighbors:
                            logger.warning("The group composed by the atoms '%s' will be ignored because the OXT atom is missing in the "
                                           "residue %s. If the atom is not the C-terminal just ignore this message. "
                                           % (atms_str.replace(",", ", "), target_residue))
                    else:
                        logger.warning("The group composed by the atoms '%s' will be ignored because some of them were not found in "
                                       "the residue %s. The missing atoms are: %s." % (atms_str.replace(",", ", "), target_residue,
                                                                                       ", ".join(missing_atoms)))
                else:
                    # It checks if a cysteine atom is establishing a disulfide bond. If it does, it will read a predefined set of
                    # properties (SS_BOND_FEATURES). In this case, the SG is hydrophobic and is not a donor anymore.
                    if target_residue.is_aminoacid() and target_residue.resname == "CYS" and nb_atoms[0] in atms_in_ss_bonds:
                        features = [ChemicalFeature(f) for f in SS_BOND_FEATURES]
                    else:
                        features = [ChemicalFeature(f) for f in self.default_properties[target_residue.resname][atms_str].split(",")]

                    # It considers all first residues in the chain as the N-terminal and won't evaluate if there are missing residues
                    # in the N-terminal. Since not all PDBs will have a 'REMARK 465 MISSING RESIDUES' field and, also, given that this
                    # field is not reliable enough, the best way to deal with the N-terminal would be by aligning the structure with
                    # the protein sequence. However, it would force one to provide the sequence and it would also increase the processing
                    # time only to identify if the residue is an N-terminal or not. In most of the cases, the first residue will be the
                    # N-terminal. Moreover, in general (maybe always), the binding sites will not be at the ends of the protein.
                    # Therefore, I opted to always attribute a 'PositivelyIonizable' property to the N atom from the first
                    # residue in the chain.
                    if atms_str == "N" and target_residue.idx == 0:
                        features.append(ChemicalFeature("PositivelyIonizable"))

                    atm_grp = AtomGroup(nb_atoms, list(set(features)), recursive=True)
                    comp_grps.add_atm_grps([atm_grp])

            # Check for amides only when the target is an amino acid.
            if target_residue.is_aminoacid():
                for nb in neighbors:
                    # Ignore residues that are not amino acids.
                    if not nb.is_aminoacid():
                        continue
                    # Recover valid atoms by applying an atom selection.
                    nb_atms = {atm.name: atm for atm in nb.get_atoms() if atm_sel.accept_atom(atm)}
                    if nb.idx < target_residue.idx:
                        # If this residue is neighbor of the target residue N, then they form an amide (N-C=O).
                        if "N" in atm_map and "O" in nb_atms and "C" in nb_atms and atm_map["N"].is_neighbor(nb_atms["C"]):
                            # Define an amide group and add it to the compound groups.
                            # In the atm_map, the atoms were already transformed into NbAtoms().
                            amide = [atm_map["N"], NbAtom(nb_atms["C"]), NbAtom(nb_atms["O"])]
                            atm_grp = AtomGroup(amide, [ChemicalFeature("Amide")], recursive=True)
                            comp_grps.add_atm_grps([atm_grp])
                    elif nb.idx > target_residue.idx:
                        # If this residue is neighbor of the target residue C, then they form an amide (N-C=O).
                        if "N" in nb_atms and "O" in atm_map and "C" in atm_map and atm_map["C"].is_neighbor(nb_atms["N"]):
                            # Define an amide group and add it to the compound groups.
                            # In the atm_map, the atoms were already transformed into NbAtoms().
                            amide = [NbAtom(nb_atms["N"]), atm_map["C"], atm_map["O"]]
                            atm_grp = AtomGroup(amide, [ChemicalFeature("Amide")], recursive=True)
                            comp_grps.add_atm_grps([atm_grp])

            # Define a new graph as a dict object: each key is a serial number and the values are dictionaries with the serial
            # number of each neighbor of an atom. The algorithm of Bellman Ford requires weights, so we set all weights to 1.
            nb_graph = {atm.serial_number: None for atm in atm_map.values()}
            for atm in atm_map.values():
                # Add the neighbors of this atom. The number '1' defined below represents the edge weight. But, right now it's not used.
                nb_graph[atm.serial_number] = {i.serial_number: 1 for i in atm.neighbors_info if i.serial_number in nb_graph}
                # Set the graph dictionary to each NbAtom object. Since, dictionaries are passed as references, we can change
                # the variable nb_graph and the changes will be updated in each NbAtom object automatically.
                atm.set_neighborhood(nb_graph)

            return comp_grps

    def _get_calculated_properties(self, target_residue, target_atoms, residue_selector, mol_obj=None):

        # If no OBMol was defined, create a new one through the residue list.
        mol_obj = mol_obj or self._get_mol_from_entity(target_residue.get_parent_by_level('M'), residue_selector)
        # Wrap the molecule object.
        mol_obj = MolWrapper(mol_obj)

        if mol_obj.get_num_heavy_atoms() != len(target_atoms):
            raise MoleculeSizeError("The number of heavy atoms in the PDB selection and in the MOL file are different.")

        # Ignore hydrogen atoms
        atm_obj_list = [atm for atm in mol_obj.get_atoms(wrapped=True) if atm.get_atomic_num() != 1]

        atm_map = {}
        trgt_atms = {}
        for i, atm_obj in enumerate(atm_obj_list):
            atm_map[atm_obj.get_idx()] = target_atoms[i].serial_number
            trgt_atms[target_atoms[i].serial_number] = NbAtom(target_atoms[i])
        # Set all neighbors, i.e., covalently bonded atoms.
        for bond_obj in mol_obj.get_bonds(wrapped=True):
            bgn_atm_obj = bond_obj.get_begin_atom(wrapped=True)
            end_atm_obj = bond_obj.get_end_atom(wrapped=True)

            # At least one of the atoms must be a non-hydrogen atom.
            if bgn_atm_obj.get_atomic_num() != 1 or end_atm_obj.get_atomic_num() != 1:
                # If the atom 1 is not a hydrogen, add atom 2 to its neighbor list.
                if bgn_atm_obj.get_atomic_num() != 1:
                    serial_number = atm_map.get(end_atm_obj.get_idx())
                    coord = mol_obj.get_atom_coord_by_id(end_atm_obj.get_id())
                    atom_info = NbAtomData(end_atm_obj.get_atomic_num(), coord, serial_number)
                    trgt_atms[atm_map[bgn_atm_obj.get_idx()]].add_nb_info([atom_info])
                # If the atom 2 is not a hydrogen, add atom 1 to its neighbor list.
                if end_atm_obj.get_atomic_num() != 1:
                    serial_number = atm_map.get(bgn_atm_obj.get_idx())
                    coord = mol_obj.get_atom_coord_by_id(bgn_atm_obj.get_id())
                    atom_info = NbAtomData(bgn_atm_obj.get_atomic_num(), coord, serial_number)
                    trgt_atms[atm_map[end_atm_obj.get_idx()]].add_nb_info([atom_info])

        group_features = self.feature_extractor.get_features_by_groups(mol_obj, atm_map)

        comp_grps = CompoundGroups(target_residue)
        for key in group_features:
            grp_obj = group_features[key]

            # It only accepts groups containing at least one atom from the target residue.
            if any([trgt_atms[i].parent == target_residue for i in grp_obj["atm_ids"]]) is False:
                continue

            atoms = [trgt_atms[i] for i in grp_obj["atm_ids"]]
            atm_grp = AtomGroup(atoms, grp_obj["features"], recursive=True)
            comp_grps.add_atm_grps([atm_grp])

        # Define a new graph as a dict object: each key is a serial number and the values are dictionaries with the serial
        # number of each neighbor of an atom. The algorithm of Bellman Ford requires weights, so we set all weights to 1.
        nb_graph = {atm.serial_number: None for atm in trgt_atms.values()}
        for atm in trgt_atms.values():
            # Add the neighbors of this atom. The number '1' defined below represents the edge weight. But, right now it's not used.
            nb_graph[atm.serial_number] = {i.serial_number: 1 for i in atm.neighbors_info if i.serial_number in nb_graph}
            # Set the graph dictionary to each NbAtom object. Since, dictionaries are passed as references, we can change
            # the variable nb_graph and the changes will be updated in each NbAtom object automatically.
            atm.set_neighborhood(nb_graph)

        return comp_grps

    def _get_mol_from_entity(self, entity, residue_selector):
        # First it saves the selection into a PDB file and then it converts the file to .mol.
        # I had to do it because the OpenBabel 2.4.1 had a problem with some molecules containing aromatic rings.
        # In such cases, the aromatic ring was being wrongly perceived and some atoms received more double bonds than it
        # was expected. The version 2.3.2 works better. Therefore, I defined this version manually (openbabel property).
        filename = get_unique_filename(self.tmp_path)
        pdb_file = '%s.pdb' % filename
        logger.info("Saving the PDB object as a PDB file named '%s'." % pdb_file)
        save_to_file(entity, pdb_file, residue_selector)

        mol_file = '%s.mol' % filename
        ob_opt = {"error-level": 5}
        logger.info("Converting the PDB file '%s' to a Mol file named '%s' using Open Babel." % (pdb_file, mol_file))
        if self.add_h:
            logger.info("Hydrogens will be added to the molecule.")
            if self.ph is not None:
                ob_opt["p"] = self.ph
            else:
                ob_opt["h"] = ""
        convert_molecule(pdb_file, mol_file, opt=ob_opt, openbabel=OPENBABEL)

        mol_obj = None
        if self.amend_mol:
            logger.info("A validation will be performed and it will try to fix some errors.")

            try:
                mol_obj = next(readfile("mol", mol_file))
            except Exception as e:
                logger.exception(e)
                raise MoleculeObjectError("An error occurred while parsing the file '%s' with Open Babel and the molecule "
                                          "object could not be created. Check the logs for more information." % mol_file)

            mv = MolValidator(fix_charges=True, fix_implicit_valence=True)
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

        return mol_obj

    def _merge_hydrophobes(self, comp_grps):
        # Only hydrophobic atom groups.
        hydrop_atm_grps = [g for g in comp_grps.atm_grps if "Hydrophobic" in g.feature_names]
        # Hydrophobic islands dictionary. Keys are integer values and items are defined by a set of atom groups.
        hydrop_islands = defaultdict(set)
        # It stores a mapping of an atom (represented by its serial number) and a hydrophobic island (defined by its keys).
        atm_mapping = {}

        grp_id = 0
        for atm_grp in hydrop_atm_grps:
            # Hydrophobic atoms are defined always as only one atom.
            atm = atm_grp.atoms[0]

            # Recover the groups of all neighbors of this atom (it will merge all existing groups).
            nb_grps = set([atm_mapping[nb] for nb in atm.neighborhood[atm.serial_number].keys() if nb in atm_mapping])

            # Already there are hydrophobe groups formed by the neighbors of this atom.
            if nb_grps:
                # Merge all groups of the neighbors of this atom.
                new_grp = set(chain.from_iterable([hydrop_islands.pop(nb_grp_id) for nb_grp_id in nb_grps]))
                # Include this atom to the merged group.
                new_grp.add(atm)

                for k in atm_mapping:
                    if atm_mapping[k] in nb_grps:
                        atm_mapping[k] = grp_id

                hydrop_islands[grp_id] = new_grp
                atm_mapping[atm.serial_number] = grp_id
            else:
                atm_mapping[atm.serial_number] = grp_id
                hydrop_islands[grp_id].add(atm)
            grp_id += 1

        for atms in hydrop_islands.values():
            atm_grp = AtomGroup(atms, [ChemicalFeature("Hydrophobe")], recursive=True)
            comp_grps.add_atm_grps([atm_grp])

        for atm_grp in hydrop_atm_grps:
            features = [f for f in atm_grp.features if f.name != "Hydrophobic"]

            # If the atom group does not have any feature, we can remove it from the CompoundGroup object.
            # This is unlikely to occur as all atoms (by default) will have at least the feature 'Atom'. So, I added this test only
            # to avoid cases in which one edits the feature definitions in a way that an atom ends up having no features at all after
            # removing its hydrophobic feature.
            if len(features) == 0:
                comp_grps.remove_atm_grps([atm_grp])
            else:
                atm_grp.features = features
