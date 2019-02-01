from MyBio.selector import ResidueSelector
from MyBio.util import save_to_file

import mol.interaction.math as imath
from mol.interaction.contact import get_contacts_for_entity
from mol.neighborhood import (NbAtom, NbCoordinates)
from mol.coord import Coordinate
from mol.charge_model import OpenEyeModel
from mol.validator import (MolValidator, RDKitValidator)
from mol.wrappers.obabel import convert_molecule
from util.file import get_unique_filename
from util.exceptions import (MoleculeSizeError, MoleculeInformationError, IllegalArgumentError)

from rdkit.Chem import (MolFromMolFile, MolFromMolBlock)
from openbabel import (OBMolAtomIter, OBAtomAtomIter, OBMolBondIter)
from pybel import readfile

from collections import defaultdict

import logging


logger = logging.getLogger(__name__)


class AtomGroup():

    def __init__(self, atoms, features, interactions=None, recursive=True):
        self.features = features

        # Update the atoms and coordinate properties.
        self._atoms = atoms
        self._coords = imath.atom_coordinates(atoms)
        self._centroid = imath.centroid(self.coords)
        self._normal = None

        self.interactions = interactions or []

        self._hash_cache = None

        if recursive:
            for a in self.atoms:
                a.add_atom_group(self)

    @property
    def atoms(self):
        return self._atoms

    @atoms.setter
    def atoms(self, new_atoms):
        self._atoms = new_atoms
        self._coords = imath.atom_coordinates(new_atoms)
        self._centroid = imath.centroid(self.coords)
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
            self._normal = imath.calc_normal(self.coords)
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
        """Return 1 if all compounds are water molecules."""
        return any([a.parent.is_water() for a in self.atoms])

    def has_hetatm(self):
        """Return 1 if all compounds are hetero groups."""
        return any([a.parent.is_hetatm() for a in self.atoms])

    def has_aminoacid(self):
        """Return 1 if all compounds are amino acids."""
        return any([a.parent.is_aminoacid() for a in self.atoms])

    def has_nucleotide(self):
        """Return 1 if all compounds are nucleotides."""
        return any([a.parent.is_nucleotide() for a in self.atoms])

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

    def __init__(self, mybio_residue, atm_grps=None):
        self.residue = mybio_residue
        self.atm_grps = atm_grps or []

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
# Create a class for finding groups
# Option: use openbabel or rdkit

# Create a function for extracting openbabel smarts search

# Adapt code for accepting ligand files and use it
def find_compound_groups(mybio_residue, feature_extractor, output_path, add_h=False, ph=None,
                         error_level=2, charge_model=OpenEyeModel()):

    '''
        error_level:
            0 = permissive mode. No validation will be performed.
            1 = amending mode. A validation will be performed and it will try to fix some errors.
            2 = strict mode. A validation will be performed, and no errors will be accepted.
            3 = RDKit sanitization mode. A validation will be performed using RDKit sanitization functions, and no errors will be accepted.
            4 = Amending + RDKit sanitization mode. A validation will be performed using RDKit sanitization functions, but first it will try to fix errors.
    '''
    COV_CUTOFF = 2

    model = mybio_residue.get_parent_by_level('M')
    proximal = get_contacts_for_entity(entity=model, source=mybio_residue, radius=COV_CUTOFF, level='R')
    res_set = set([p[1] for p in proximal])
    res_list = sorted(list(res_set), key=lambda r: r.id)
    res_sel = ResidueSelector(res_list, keep_altloc=False, keep_hydrog=False)

    # First it saves the selection into a PDB file and then it converts the file to .mol.
    # I had to do it because the OpenBabel 2.4.1 had a problem with some molecules containing aromatic rings.
    # In such cases, the aromatic ring was being wrongly perceived and some atoms received more double bonds than
    # it was exepcted. The version 2.3.2 works fine. Therefore, I defined this version manually (openbabel property).
    # But why to convert to MOL file? Because RDKit does not perceive bond order from PDB files.
    filename = get_unique_filename(output_path)
    pdb_file = '%s.pdb' % filename
    save_to_file(model, pdb_file, res_sel)

    mol_file = '%s.mol' % filename
    ob_opt = {"error-level": 5}
    if add_h:
        if ph is not None:
            ob_opt["p"] = ph
        else:
            ob_opt["h"] = None
    convert_molecule(pdb_file, mol_file, opt=ob_opt, openbabel='/usr/bin/obabel')

    # Recover all atoms from the previous selected residues.
    atoms = [a for r in res_list for a in r.get_unpacked_list() if res_sel.accept_atom(a)]
    atoms = tuple(sorted(atoms, key=lambda a: a.serial_number))
    ob_mol = next(readfile("mol", mol_file))

    if error_level == 0:
        logger.warning("Error level set to 0. No validation will be performed. Molecules containing errors may generate incorrect results.")
    elif error_level == 1:
        logger.warning("Error level set to 1. A validation will be performed and it will try to fix some errors.")
        mv = MolValidator(fix_charges=True, fix_implicit_valence=True)
        if not mv.is_mol_valid(ob_mol):
            logger.warning("Molecule loaded from '%s' is invalid. Check the logs for more information." % mol_file)
            raise MoleculeInformationError("Molecule loaded from '%s' is invalid. Check the logs for more information." % mol_file)
    elif error_level == 2:
        logger.warning("Error level set to 2. A validation will be performed and no errors will be accepted.")
        mv = MolValidator()
        if not mv.is_mol_valid(ob_mol):
            logger.warning("Molecule loaded from '%s' is invalid. Check the logs for more information." % mol_file)
            raise MoleculeInformationError("Molecule loaded from '%s' is invalid. Check the logs for more information." % mol_file)
    elif error_level == 3:
        logger.warning("Error level set to 3. A validation will be performed using RDKit sanitization functions. No errors will be accepted.")
        mv = RDKitValidator()
        rdk_mol = MolFromMolFile(mol_file, removeHs=False, sanitize=False)
        if not mv.is_mol_valid(rdk_mol):
            logger.warning("Molecule loaded from '%s' is invalid. Check the logs for more information." % mol_file)
            raise MoleculeInformationError("Molecule loaded from '%s' is invalid. Check the logs for more information." % mol_file)
    elif error_level == 4:
        logger.warning("A validation will be performed using RDKit sanitization functions. No errors will be accepted.")
        mv = MolValidator(fix_charges=True, fix_implicit_valence=True)
        if not mv.is_mol_valid(ob_mol):
            logger.warning("Molecule loaded from '%s' is invalid. Check the logs for more information." % mol_file)
            raise MoleculeInformationError("Molecule loaded from '%s' is invalid. Check the logs for more information." % mol_file)
        else:
            mv = RDKitValidator()
            mol_block = ob_mol.write('mol')
            rdk_mol = MolFromMolBlock(mol_block, removeHs=False, sanitize=False)
            if not mv.is_mol_valid(rdk_mol):
                logger.warning("Molecule loaded from '%s' is invalid. Check the logs for more information." % mol_file)
                raise MoleculeInformationError("Molecule loaded from '%s' is invalid. Check the logs for more information." % mol_file)
    else:
        raise IllegalArgumentError("Error level must be either 0, 1, 2, 3, or 4.")
    logger.warning('Molecule validated!!!')

    if ob_mol.OBMol.NumHvyAtoms() != len(atoms):
        raise MoleculeSizeError("The number of heavy atoms in the PDB selection and in the MOL file are different.")

    atm_map = {}
    trgt_atms = {}
    for i, ob_atm in enumerate(OBMolAtomIter(ob_mol.OBMol)):
        # Ignore hydrogen atoms
        if (ob_atm.GetAtomicNum() != 1):
            atm_map[ob_atm.GetIdx()] = atoms[i].serial_number
            trgt_atms[atoms[i].serial_number] = NbAtom(atoms[i])

    # Set all neighbors, i.e., covalently bonded atoms.
    for ob_bond in OBMolBondIter(ob_mol.OBMol):
        bgn_ob_atm = ob_bond.GetBeginAtom()
        end_ob_atm = ob_bond.GetEndAtom()

        if bgn_ob_atm.GetAtomicNum() != 1 and end_ob_atm.GetAtomicNum() != 1:
            bgn_nb_atm = trgt_atms[atm_map[bgn_ob_atm.GetIdx()]]
            end_nb_atm = trgt_atms[atm_map[end_ob_atm.GetIdx()]]

            bgn_nb_atm.add_nb_atom(end_nb_atm)
            end_nb_atm.add_nb_atom(bgn_nb_atm)

    group_features = feature_extractor.get_features_by_groups(ob_mol, atm_map)

    comp_grps = CompoundGroups(mybio_residue)
    for key in group_features:
        grp_obj = group_features[key]

        # Only get groups containing atoms from the informed residue (mybio_residue parameter)
        if any([trgt_atms[i].parent == mybio_residue for i in grp_obj["atm_ids"]]):
            atoms = [trgt_atms[i] for i in grp_obj["atm_ids"]]
            atm_grp = AtomGroup(atoms, grp_obj["features"], recursive=True)
            comp_grps.add_group(atm_grp)

    return comp_grps
