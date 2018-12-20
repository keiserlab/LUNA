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
        self.atoms = atoms
        self.features = features
        self.coords = imath.atom_coordinates(atoms)
        self._centroid = imath.centroid(self.coords)
        self._normal = None

        self.interactions = interactions or []

        if recursive:
            for a in self.atoms:
                a.add_atom_group(self)

    @property
    def compound(self):
        return self.atoms[0].get_parent()

    @property
    def centroid(self):
        return self._centroid

    @property
    def normal(self):
        if self._normal is None:
            self._calc_normal()

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

    def _calc_normal(self):
        if self._normal is None:
            self._normal = imath.calc_normal(self.coords)

    def __repr__(self):
        return '<AtomGroup: [%s]' % ', '.join([str(x) for x in self.atoms])


class CompoundGroups():

    def __init__(self, mybio_residue, atm_grps=None):
        self.residue = mybio_residue

        if atm_grps is None:
            atm_grps = []

        self.atm_grps = atm_grps

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

    res_sel = ResidueSelector(list(cov_residues) + [mybio_residue])

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
            if not mybio_residue.is_water():
                for nb_obAtom in OBAtomAtomIter(obAtom):
                    coords = Coordinate(nb_obAtom.GetX(), nb_obAtom.GetY(), nb_obAtom.GetZ(),
                                        atomic_num=nb_obAtom.GetAtomicNum())

                    nb_coords.append(coords)

            nb_coords_by_atm[serialNumber] = NbCoordinates(nb_coords)

    mol_block = obMol.write('mol')
    rdMol = MolFromMolBlock(mol_block)

    group_features = feature_extractor.get_features_by_groups(rdMol, atm_map)

    trgt_atms = {}
    for atm in mybio_residue.get_atoms():
        # Ignore hydrogen atoms
        if (atm.element != "H"):
            nb_coords = nb_coords_by_atm[atm.get_serial_number()]
            nb_atom = NbAtom(atm, nb_coords)
            trgt_atms[atm.get_serial_number()] = nb_atom

    comp_grps = CompoundGroups(mybio_residue)
    for key in group_features:
        grp_obj = group_features[key]

        is_grp_valid = True
        atoms = []
        for atm_idx in grp_obj["atomIds"]:
            if (atm_idx not in trgt_atms):
                is_grp_valid = False
                break
            else:
                atoms.append(trgt_atms[atm_idx])

        # This group has atoms that do not belong to the target residue.
        if (is_grp_valid is False):
            continue

        atm_grp = AtomGroup(atoms, grp_obj["features"], recursive=True)
        comp_grps.add_group(atm_grp)

    return comp_grps
