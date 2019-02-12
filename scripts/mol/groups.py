from MyBio.selector import ResidueSelector
from MyBio.util import save_to_file

import mol.interaction.math as im
from mol.interaction.contact import get_contacts_for_entity
from mol.neighborhood import (NbAtom, NbAtomData)
from mol.charge_model import OpenEyeModel
from mol.validator import (MolValidator, RDKitValidator)
from mol.wrappers.obabel import convert_molecule
from util.file import get_unique_filename
from util.exceptions import (MoleculeSizeError, MoleculeInformationError, IllegalArgumentError)

from rdkit.Chem import MolFromMolBlock
from openbabel import (OBMolAtomIter, OBMolBondIter)
from pybel import readfile

from collections import defaultdict

import logging


logger = logging.getLogger(__name__)

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


# TODO:
# Option: use openbabel or rdkit
class CompoundGroupPerceiver():

    def __init__(self, feature_extractor, add_h=False, ph=None, error_level=4,
                 charge_model=OpenEyeModel(), tmp_path=None, expand_selection=True, radius=COV_CUTOFF):
        '''
            error_level:
                0 = permissive mode. No validation will be performed.
                1 = amending mode. A validation will be performed and it will try to fix some errors.
                2 = strict mode. A validation will be performed, and no errors will be accepted.
                3 = RDKit sanitization mode. A validation will be performed using RDKit sanitization functions, and no errors will be accepted.
                4 = Amending + RDKit sanitization mode. A validation will be performed using RDKit sanitization functions, but first it will try to fix errors.
        '''
        self.feature_extractor = feature_extractor
        self.add_h = add_h
        self.ph = ph
        self.error_level = error_level
        self.charge_model = charge_model
        self.tmp_path = tmp_path
        self.expand_selection = expand_selection,
        self.radius = radius

    def perceive_compound_groups(self, target_residue, ob_mol=None, only_grps_with_target=True):
        # If no OBMol object was defined and expand_selection was set ON, it will get all residues around the
        # target and create a new OBMol object with them.
        if ob_mol is None and self.expand_selection:
            res_list = self._get_all_around_residue(target_residue)
        # Otherwise, the residue list will be composed only by the target residue.
        else:
            res_list = [target_residue]

        # Select residues based on the residue list defined previously.
        res_sel = ResidueSelector(res_list, keep_altloc=False, keep_hydrog=False)
        # If no OBMol was defined, create a new one through the residue list.
        ob_mol = ob_mol or self._get_obmol_from_entity(target_residue.get_parent_by_level('M'), res_sel)

        # Validate ob_mol object according to the error_level informed when initializing this class.
        self._validate_obmol(ob_mol)

        # Recover all atoms from the previous selected residues.
        atoms = [a for r in res_list for a in r.get_unpacked_list() if res_sel.accept_atom(a)]
        atoms = tuple(sorted(atoms, key=lambda a: a.serial_number))

        if ob_mol.OBMol.NumHvyAtoms() != len(atoms):
            raise MoleculeSizeError("The number of heavy atoms in the PDB selection and in the MOL file are different.")

        # Ignore hydrogen atoms
        ob_atms = [atm for atm in OBMolAtomIter(ob_mol.OBMol) if atm.GetAtomicNum() != 1]

        atm_map = {}
        trgt_atms = {}
        for i, ob_atm in enumerate(ob_atms):
            atm_map[ob_atm.GetIdx()] = atoms[i].serial_number
            trgt_atms[atoms[i].serial_number] = NbAtom(atoms[i])

        # Set all neighbors, i.e., covalently bonded atoms.
        for ob_bond in OBMolBondIter(ob_mol.OBMol):
            bgn_ob_atm = ob_bond.GetBeginAtom()
            end_ob_atm = ob_bond.GetEndAtom()

            # At least one of the atoms must be a non-hydrogen atom.
            if bgn_ob_atm.GetAtomicNum() != 1 or end_ob_atm.GetAtomicNum() != 1:
                # If the atom 1 is not a hydrogen, add atom 2 to its neighbor list.
                if bgn_ob_atm.GetAtomicNum() != 1:
                    serial_number = atm_map.get(end_ob_atm.GetIdx())
                    coord = [end_ob_atm.GetX(), end_ob_atm.GetY(), end_ob_atm.GetZ()]
                    atom_info = NbAtomData(end_ob_atm.GetAtomicNum(), coord, serial_number)
                    trgt_atms[atm_map[bgn_ob_atm.GetIdx()]].add_nb_atom(atom_info)
                # If the atom 2 is not a hydrogen, add atom 1 to its neighbor list.
                if end_ob_atm.GetAtomicNum() != 1:
                    serial_number = atm_map.get(bgn_ob_atm.GetIdx())
                    coord = [bgn_ob_atm.GetX(), bgn_ob_atm.GetY(), bgn_ob_atm.GetZ()]
                    atom_info = NbAtomData(bgn_ob_atm.GetAtomicNum(), coord, serial_number)
                    trgt_atms[atm_map[end_ob_atm.GetIdx()]].add_nb_atom(atom_info)

        group_features = self.feature_extractor.get_features_by_groups(ob_mol, atm_map)

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

    def _get_obmol_from_entity(self, entity, residue_selector):
        # First it saves the selection into a PDB file and then it converts the file to .mol.
        # I had to do it because the OpenBabel 2.4.1 had a problem with some molecules containing aromatic rings.
        # In such cases, the aromatic ring was being wrongly perceived and some atoms received more double bonds than it
        # was exepcted. The version 2.3.2 works fine. Therefore, I defined this version manually (openbabel property).
        filename = get_unique_filename(self.tmp_path)
        pdb_file = '%s.pdb' % filename
        save_to_file(entity, pdb_file, residue_selector)

        mol_file = '%s.mol' % filename
        ob_opt = {"error-level": 5}
        if self.add_h:
            if self.ph is not None:
                ob_opt["p"] = self.ph
            else:
                ob_opt["h"] = ""
        convert_molecule(pdb_file, mol_file, opt=ob_opt, openbabel='/usr/bin/obabel')
        # Read molecule file and return an OBMol object.
        return next(readfile("mol", mol_file))

    def _validate_obmol(self, ob_mol):
        if self.error_level == 0:
            logger.warning("Error level set to 0. No validation will be performed. "
                           "Molecules containing errors may generate incorrect results.")
        elif self.error_level == 1:
            logger.warning("Error level set to 1. A validation will be performed and it will try to fix some errors.")
            mv = MolValidator(fix_charges=True, fix_implicit_valence=True)
            if not mv.is_mol_valid(ob_mol):
                raise MoleculeInformationError("Informed molecule '%s' is invalid. "
                                               "Check the logs for more information." % ob_mol.OBMol.GetTitle())
        elif self.error_level == 2:
            logger.warning("Error level set to 2. A validation will be performed and no errors will be accepted.")
            mv = MolValidator()
            if not mv.is_mol_valid(ob_mol):
                raise MoleculeInformationError("Informed molecule '%s' is invalid. "
                                               "Check the logs for more information." % ob_mol.OBMol.GetTitle())
        elif self.error_level == 3:
            logger.warning("Error level set to 3. A validation will be performed using RDKit sanitization functions. "
                           "No errors will be accepted.")
            mv = RDKitValidator()
            rdk_mol = MolFromMolBlock(ob_mol.write('mol'), removeHs=False, sanitize=False)
            if not mv.is_mol_valid(rdk_mol):
                raise MoleculeInformationError("Informed molecule '%s' is invalid. "
                                               "Check the logs for more information." % ob_mol.OBMol.GetTitle())
        elif self.error_level == 4:
            logger.warning("A validation will be performed using RDKit sanitization functions. "
                           "No errors will be accepted.")
            mv = MolValidator(fix_charges=True, fix_implicit_valence=True)
            if not mv.is_mol_valid(ob_mol):
                raise MoleculeInformationError("Informed molecule '%s' is invalid. "
                                               "Check the logs for more information." % ob_mol.OBMol.GetTitle())
            else:
                mv = RDKitValidator()
                rdk_mol = MolFromMolBlock(ob_mol.write('mol'), removeHs=False, sanitize=False)
                if not mv.is_mol_valid(rdk_mol):
                    raise MoleculeInformationError("Informed molecule '%s' is invalid. "
                                                   "Check the logs for more information." % ob_mol.OBMol.GetTitle())
        else:
            raise IllegalArgumentError("Error level must be either 0, 1, 2, 3, or 4.")
        logger.warning('Molecule validated!!!')
