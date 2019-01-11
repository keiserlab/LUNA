from MyBio.PDB.PDBIO import PDBIO
from MyBio.selector import ResidueSelector

import mol.interaction.math as imath
from mol.interaction.contact import get_contacts_for_entity
from mol.neighborhood import (NbAtom, NbCoordinates)
from mol.coord import Coordinate

from rdkit.Chem import MolFromMolBlock
from openbabel import (OBMolAtomIter, OBAtomAtomIter)

from collections import defaultdict
from pybel import readstring
from io import StringIO


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
                a.add_atm_grp(self)

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


def find_compound_groups(mybio_residue, feature_extractor, ph=None, has_explicit_hydrogen=False):

    COV_CUTOFF = 1.99999999999

    model = mybio_residue.get_parent_by_level('M')
    proximal = get_contacts_for_entity(entity=model, source=mybio_residue,
                                       radius=COV_CUTOFF, level='R')

    cov_residues = set()
    for pair in proximal:
        if (pair[1] != mybio_residue):
            cov_residues.add(pair[1])

    res_list = list(cov_residues) + [mybio_residue]
    res_sel = ResidueSelector(res_list)

    fh = StringIO()
    io = PDBIO()
    io.set_structure(model)
    io.save(fh, select=res_sel, preserve_atom_numbering=True, write_conects=True)
    fh.seek(0)
    pdb_block = ''.join(fh.readlines())

    obMol = readstring("pdb", pdb_block)
    if (ph is not None):
        obMol.OBMol.AddHydrogens(True, True, ph)

    atm_map = {}
    nb_coords_by_atm = {}

    filtered_obAtom = [a for a in OBMolAtomIter(obMol.OBMol) if a.GetAtomicNum() != 1]

    for idx, obAtom in enumerate(filtered_obAtom):
        # Ignore hydrogen atoms
        if (obAtom.GetAtomicNum() != 1):
            obResidue = obAtom.GetResidue()
            serialNumber = obResidue.GetSerialNum(obAtom)
            atm_map[idx] = serialNumber

            nb_coords = []
            for nb_obAtom in OBAtomAtomIter(obAtom):
                coords = Coordinate(nb_obAtom.GetX(), nb_obAtom.GetY(), nb_obAtom.GetZ(),
                                    atomic_num=nb_obAtom.GetAtomicNum())

                nb_coords.append(coords)

            nb_coords_by_atm[serialNumber] = NbCoordinates(nb_coords)

    mol_block = obMol.write('mol')
    rdMol = MolFromMolBlock(mol_block)

    group_features = feature_extractor.get_features_by_groups(rdMol, atm_map)

    trgt_atms = {}
    for mybio_res in res_list:
        for atm in mybio_res.get_atoms():
            # Ignore hydrogen atoms
            if (atm.element != "H"):
                nb_coords = nb_coords_by_atm[atm.get_serial_number()]
                nb_atom = NbAtom(atm, nb_coords)
                trgt_atms[atm.get_serial_number()] = nb_atom

    comp_grps = CompoundGroups(mybio_residue)
    for key in group_features:
        grp_obj = group_features[key]

        # Only get groups containing atoms from the informed residue (mybio_residue parameter)
        if any([trgt_atms[i].parent == mybio_residue for i in grp_obj["atm_ids"]]):
            atoms = [trgt_atms[i] for i in grp_obj["atm_ids"]]
            atm_grp = AtomGroup(atoms, grp_obj["features"], recursive=True)
            comp_grps.add_group(atm_grp)

    return comp_grps
