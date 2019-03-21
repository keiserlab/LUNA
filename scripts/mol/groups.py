from MyBio.selector import ResidueSelector
from MyBio.util import save_to_file

import mol.interaction.math as im

from mol.interaction.contact import get_contacts_for_entity
from mol.neighborhood import (NbAtom, NbAtomData)
from mol.charge_model import OpenEyeModel
from mol.validator import MolValidator
from mol.wrappers.obabel import convert_molecule
from mol.wrappers.base import MolWrapper
from util.file import get_unique_filename
from util.exceptions import MoleculeSizeError, IllegalArgumentError, MoleculeObjectError
from util.default_values import ACCEPTED_MOL_OBJ_TYPES

from rdkit.Chem import MolFromMolBlock, MolFromMolFile, SanitizeFlags, SanitizeMol
from pybel import readfile

from collections import defaultdict

import logging
logger = logging.getLogger()


COV_CUTOFF = 2


class AtomGroup():

    def __init__(self, atoms, features, interactions=None, recursive=True):
        self.features = features

        # Update the atoms and coordinate properties.
        self._atoms = atoms
        self._coords = im.atom_coordinates(atoms)
        self._centroid = im.centroid(self.coords)
        self._normal = None

        self.interactions = interactions or []

        self._hash_cache = None

        if recursive:
            for a in self.atoms:
                a.add_atm_grp(self)

    @property
    def atoms(self):
        return self._atoms

    @atoms.setter
    def atoms(self, new_atoms):
        self._atoms = new_atoms
        self._coords = im.atom_coordinates(new_atoms)
        self._centroid = im.centroid(self.coords)
        self._normal = None
        self._hash_cache = None

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

    def add_interaction(self, interaction):
        self.interactions = list(set(self.interactions + [interaction]))

    def remove_interaction(self, interaction):
        self.interactions.remove(interaction)

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

    def __repr__(self):
        return '<AtomGroup: [%s]' % ', '.join([str(x) for x in self.atoms])

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


class CompoundGroups():

    def __init__(self, target_residue, atm_grps=None, only_grps_with_target=True):
        self.residue = target_residue
        self.atm_grps = atm_grps or []
        self.only_grps_with_target = only_grps_with_target

    @property
    def summary(self):
        summary = defaultdict(int)
        for grp in self.atm_grps:
            for feature in grp.features:
                summary[feature] += 1

        return summary

    def add_group(self, group):
        atm_grps = set(self.atm_grps + [group])
        self.atm_grps = list(atm_grps)


class CompoundGroupPerceiver():

    def __init__(self, feature_extractor, add_h=False, ph=None, amend_mol=True, mol_obj_type="rdkit",
                 charge_model=OpenEyeModel(), tmp_path=None, expand_selection=True, radius=COV_CUTOFF):

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
        self.expand_selection = expand_selection,
        self.radius = radius

    def perceive_compound_groups(self, target_residue, mol_obj=None, only_grps_with_target=True):
        # If no OBMol object was defined and expand_selection was set ON, it will get all residues around the
        # target and create a new OBMol object with them.
        if mol_obj is None and self.expand_selection:
            res_list = self._get_all_around_residue(target_residue)
        # Otherwise, the residue list will be composed only by the target residue.
        else:
            res_list = [target_residue]

        # Select residues based on the residue list defined previously.
        res_sel = ResidueSelector(res_list, keep_altloc=False, keep_hydrog=False)
        # If no OBMol was defined, create a new one through the residue list.
        mol_obj = mol_obj or self._get_mol_from_entity(target_residue.get_parent_by_level('M'), res_sel)

        # Recover all atoms from the previous selected residues.
        atoms = [a for r in res_list for a in r.get_unpacked_list() if res_sel.accept_atom(a)]
        atoms = tuple(sorted(atoms, key=lambda a: a.serial_number))

        mol_obj = MolWrapper(mol_obj)

        if mol_obj.get_num_heavy_atoms() != len(atoms):
            raise MoleculeSizeError("The number of heavy atoms in the PDB selection and in the MOL file are different.")

        # Ignore hydrogen atoms
        atm_obj_list = [atm for atm in mol_obj.get_atoms(wrappered=True) if atm.get_atomic_num() != 1]

        atm_map = {}
        trgt_atms = {}
        for i, atm_obj in enumerate(atm_obj_list):
            atm_map[atm_obj.get_idx()] = atoms[i].serial_number
            trgt_atms[atoms[i].serial_number] = NbAtom(atoms[i])

        # Set all neighbors, i.e., covalently bonded atoms.
        for bond_obj in mol_obj.get_bonds(wrappered=True):
            bgn_atm_obj = bond_obj.get_begin_atom(wrappered=True)
            end_atm_obj = bond_obj.get_end_atom(wrappered=True)

            # At least one of the atoms must be a non-hydrogen atom.
            if bgn_atm_obj.get_atomic_num() != 1 or end_atm_obj.get_atomic_num() != 1:
                # If the atom 1 is not a hydrogen, add atom 2 to its neighbor list.
                if bgn_atm_obj.get_atomic_num() != 1:
                    serial_number = atm_map.get(end_atm_obj.get_idx())
                    coord = mol_obj.get_atom_coord_by_id(end_atm_obj.get_id())
                    atom_info = NbAtomData(end_atm_obj.get_atomic_num(), coord, serial_number)
                    trgt_atms[atm_map[bgn_atm_obj.get_idx()]].add_nb_atom(atom_info)
                # If the atom 2 is not a hydrogen, add atom 1 to its neighbor list.
                if end_atm_obj.get_atomic_num() != 1:
                    serial_number = atm_map.get(bgn_atm_obj.get_idx())
                    coord = mol_obj.get_atom_coord_by_id(bgn_atm_obj.get_id())
                    atom_info = NbAtomData(bgn_atm_obj.get_atomic_num(), coord, serial_number)
                    trgt_atms[atm_map[end_atm_obj.get_idx()]].add_nb_atom(atom_info)

        group_features = self.feature_extractor.get_features_by_groups(mol_obj, atm_map)

        comp_grps = CompoundGroups(target_residue, only_grps_with_target=only_grps_with_target)
        for key in group_features:
            grp_obj = group_features[key]

            # If has_target is set ON, it will only accept groups containing at least one atom from the target residue.
            if (only_grps_with_target and
                    any([trgt_atms[i].parent == target_residue for i in grp_obj["atm_ids"]]) is False):
                continue

            atoms = [trgt_atms[i] for i in grp_obj["atm_ids"]]
            atm_grp = AtomGroup(atoms, grp_obj["features"], recursive=True)
            comp_grps.add_group(atm_grp)

        return comp_grps

    def _get_all_around_residue(self, target_residue):
        model = target_residue.get_parent_by_level('M')
        proximal = get_contacts_for_entity(entity=model, source=target_residue, radius=self.radius, level='R')
        res_list = sorted(list(set([p[1] for p in proximal])), key=lambda r: r.id)

        return res_list

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
        convert_molecule(pdb_file, mol_file, opt=ob_opt, openbabel='/usr/bin/obabel')

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

            if self.mol_obj_type == "openbabel":
                return mol_obj
            else:
                try:
                    # The sanitization is set off. We will apply it in the next statement.
                    rdk_mol = MolFromMolBlock(mol_obj.write('mol'), sanitize=False, removeHs=False)
                    # Sanitize molecule is applied now, so we will be able to catch the exceptions raised by RDKit,
                    # otherwise it would not be possible.
                    SanitizeMol(rdk_mol, SanitizeFlags.SANITIZE_ALL)
                    return rdk_mol
                except Exception as e:
                    logger.exception(e)
                    raise MoleculeObjectError("An error occurred while parsing the molecular block generated by Open Babel with "
                                              "RDKit. Check the logs for more information." % mol_file)
        else:
            try:
                # Create a new Mol object.
                if self.mol_obj_type == "openbabel":
                    return next(readfile("mol", mol_file))
                else:
                    # The sanitization is set off. We will apply it in the next statement.
                    rdk_mol = MolFromMolFile(mol_file, sanitize=False, removeHs=False)
                    # Sanitize molecule is applied now, so we will be able to catch the exceptions raised by RDKit,
                    # otherwise it would not be possible.
                    SanitizeMol(rdk_mol, SanitizeFlags.SANITIZE_ALL)
                    return rdk_mol
            except Exception as e:
                logger.exception(e)
                tool = "Open Babel" if self.mol_obj_type == "openbabel" else "RDKit"
                raise MoleculeObjectError("An error occurred while parsing the file '%s' with %s and the molecule "
                                          "object could not be created. Check the logs for more information." % (mol_file, tool))
