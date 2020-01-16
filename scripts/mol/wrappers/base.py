import numpy as np
from rdkit.Chem import Atom as RDAtom
from rdkit.Chem import Mol as RDMol
from rdkit.Chem import Bond as RDBond
from rdkit.Chem import MolToSmiles, MolToPDBBlock, GetPeriodicTable
from openbabel import OBMol, OBAtom, OBBond, OBMolAtomIter, OBAtomAtomIter, OBAtomBondIter, OBMolBondIter
from openbabel import etab
from pybel import Molecule as PybelMol


from util.exceptions import AtomObjectTypeError, MoleculeObjectTypeError


import logging

logger = logging.getLogger()


class AtomWrapper:

    def __init__(self, atm_obj):
        if isinstance(atm_obj, self.__class__):
            atm_obj = atm_obj.unwrap()

        if not isinstance(atm_obj, RDAtom) and not isinstance(atm_obj, OBAtom):
            logger.exception("Objects of type '%s' are not currently accepted." % atm_obj.__class__)
            raise AtomObjectTypeError("Objects of type '%s' are not currently accepted." % atm_obj.__class__)

        self._atm_obj = atm_obj

    @property
    def atm_obj(self):
        return self._atm_obj

    @atm_obj.setter
    def atm_obj(self, atm_obj):
        if not isinstance(atm_obj, RDAtom) and not isinstance(atm_obj, OBAtom):
            logger.exception("Objects of type '%s' are not currently accepted." % atm_obj.__class__)
            raise AtomObjectTypeError("Objects of type '%s' are not currently accepted." % atm_obj.__class__)
        else:
            self._atm_obj = atm_obj

    def get_idx(self):
        # Both RDKit and Openbabel have the same function name.
        return self._atm_obj.GetIdx()

    def get_id(self):
        if self.is_rdkit_obj():
            return self._atm_obj.GetIdx()
        elif self.is_openbabel_obj():
            return self._atm_obj.GetId()

    def get_atomic_num(self):
        # Both RDKit and Openbabel have the same function name.
        return self._atm_obj.GetAtomicNum()

    def get_symbol(self):
        if self.is_rdkit_obj():
            return self._atm_obj.GetSymbol()
        elif self.is_openbabel_obj():
            return etab.GetSymbol(self._atm_obj.GetAtomicNum())

    def get_charge(self):
        # Both RDKit and Openbabel have the same function name.
        return self._atm_obj.GetFormalCharge()

    def get_isotope(self):
        if self.is_rdkit_obj():
            return int(round(self._atm_obj.GetMass()))
        elif self.is_openbabel_obj():
            return int(round(self._atm_obj.GetExactMass()))

    def get_mass(self):
        if self.is_rdkit_obj():
            return self._atm_obj.GetMass()
        elif self.is_openbabel_obj():
            return self._atm_obj.GetExactMass()

    def get_atomic_mass(self):
        if self.is_rdkit_obj():
            return GetPeriodicTable().GetAtomicWeight(self.get_atomic_num())
        elif self.is_openbabel_obj():
            return self._atm_obj.GetAtomicMass()

    def get_neighbors_number(self, only_heavy_atoms=False):
        # Subtract the number of hydrogens from the degree if only heavy atoms must be considered.
        penalty = 0 if only_heavy_atoms is False else self.get_h_count()
        return self.get_degree() - penalty

    def get_neighbors(self):
        # TODO: wrap option
        if self.is_rdkit_obj():
            return self._atm_obj.GetNeighbors()
        elif self.is_openbabel_obj():
            return OBAtomAtomIter(self._atm_obj)

    def get_degree(self):
        if self.is_rdkit_obj():
            # In RDKit, GetDegree() returns the degree of the atom.
            return self._atm_obj.GetTotalDegree()
        elif self.is_openbabel_obj():
            # In OpenBabel, GetImplicitValence() returns the maximum number of connections expected (degree).
            return self._atm_obj.GetImplicitValence()

    def get_valence(self):
        if self.is_rdkit_obj():
            return self._atm_obj.GetExplicitValence() + self._atm_obj.GetImplicitValence()
        elif self.is_openbabel_obj():
            bonds = [b.GetBondOrder() for b in OBAtomBondIter(self._atm_obj)]
            return sum(bonds) + self._atm_obj.ImplicitHydrogenCount()

    def get_h_count(self):
        if self.is_rdkit_obj():
            return self._atm_obj.GetTotalNumHs(includeNeighbors=True)
        elif self.is_openbabel_obj():
            return self._atm_obj.ImplicitHydrogenCount() + self._atm_obj.ExplicitHydrogenCount()

    def get_bonds(self, wrapped=True):
        bonds = []
        if self.is_rdkit_obj():
            # TODO: IMPLEMENT
            bonds = []
        elif self.is_openbabel_obj():
            bonds = OBAtomBondIter(self._atm_obj)

        if bonds and wrapped:
            return [BondWrapper(bond) for bond in bonds]
        return bonds

    def get_atomic_invariants(self):
        return [self.get_neighbors_number(only_heavy_atoms=True),       # Number of heavy atoms
                (self.get_valence() - self.get_h_count()),              # Valence - Num. Hs
                self.get_atomic_num(),                                  # Atomic number
                self.get_isotope(),                                     # Isotope number
                self.get_charge(),                                      # Formal charge
                self.get_h_count(),                                     # Num. Hs
                int(self.is_in_ring())]                                 # If the atom belongs to a ring or not

    def is_in_ring(self):
        # Both RDKit and Openbabel have the same function name.
        return self._atm_obj.IsInRing()

    def set_charge(self, charge):
        # Both RDKit and Openbabel have the same function name.
        self._atm_obj.SetFormalCharge(charge)

    def set_implicit_valence(self, valence):
        if self.is_rdkit_obj():
            logger.warning("RDKit does not provide a setter for implicit valence.")
        elif self.is_openbabel_obj():
            self._atm_obj.SetImplicitValence(valence)

    def unwrap(self):
        return self._atm_obj

    def is_rdkit_obj(self):
        if isinstance(self._atm_obj, RDAtom):
            return True
        return False

    def is_openbabel_obj(self):
        if isinstance(self._atm_obj, OBAtom):
            return True
        return False

    def __getattr__(self, attr):
        return getattr(self._atm_obj, attr)


class BondWrapper:
    def __init__(self, bond_obj):
        if isinstance(bond_obj, self.__class__):
            bond_obj = bond_obj.unwrap()

        if not isinstance(bond_obj, RDBond) and not isinstance(bond_obj, OBBond):
            logger.exception("Objects of type '%s' are not currently accepted." % bond_obj.__class__)
            raise MoleculeObjectTypeError("Objects of type '%s' are not currently accepted." % bond_obj.__class__)

        self._bond_obj = bond_obj

    @property
    def bond_obj(self):
        return self._bond_obj

    @bond_obj.setter
    def bond_obj(self, bond_obj):
        if not isinstance(bond_obj, RDBond) and not isinstance(bond_obj, OBBond):
            logger.exception("Objects of type '%s' are not currently accepted." % bond_obj.__class__)
            raise AtomObjectTypeError("Objects of type '%s' are not currently accepted." % bond_obj.__class__)
        else:
            self._bond_obj = bond_obj

    def get_partner_atom(self, atm, wrapped=True):
        if isinstance(atm, AtomWrapper):
            atm = atm.unwrap()

        partner = None
        if atm == self._bond_obj.GetBeginAtom():
            partner = self._bond_obj.GetEndAtom()
        elif atm == self._bond_obj.GetEndAtom():
            partner = self._bond_obj.GetBeginAtom()

        if partner and wrapped:
            return AtomWrapper(partner)
        return partner

    def get_begin_atom(self, wrapped=True):
        # Both RDKit and Openbabel have the same function name.
        if wrapped:
            return AtomWrapper(self._bond_obj.GetBeginAtom())
        return self._bond_obj.GetBeginAtom()

    def get_end_atom(self, wrapped=True):
        # Both RDKit and Openbabel have the same function name.
        if wrapped:
            return AtomWrapper(self._bond_obj.GetEndAtom())
        return self._bond_obj.GetEndAtom()

    def get_bond_type(self):
        if self.is_rdkit_obj():
            return self._bond_obj.GetBondType()
        elif self.is_openbabel_obj():
            return self._bond_obj.GetBondOrder()

    def set_bond_type(self, bond_type):
        if self.is_rdkit_obj():
            pass
            # return self._bond_obj.GetBondType()
        elif self.is_openbabel_obj():
            return self._bond_obj.SetBondOrder(bond_type)

    def unwrap(self):
        return self._bond_obj

    def is_rdkit_obj(self):
        if isinstance(self._bond_obj, RDBond):
            return True
        return False

    def is_openbabel_obj(self):
        if isinstance(self._bond_obj, OBBond):
            return True
        return False

    def __getattr__(self, attr):
        return getattr(self._bond_obj, attr)


class MolWrapper:
    def __init__(self, mol_obj):
        if isinstance(mol_obj, self.__class__):
            mol_obj = mol_obj.unwrap()
        elif isinstance(mol_obj, PybelMol):
            mol_obj = mol_obj.OBMol

        if not isinstance(mol_obj, RDMol) and not isinstance(mol_obj, OBMol):
            logger.exception("Objects of type '%s' are not currently accepted." % mol_obj.__class__)
            raise MoleculeObjectTypeError("Objects of type '%s' are not currently accepted." % mol_obj.__class__)

        self._mol_obj = mol_obj

    @classmethod
    def from_mol_file(self, mol_file, mol_file_ext=None, mol_obj_type="rdkit"):
        # TODO: implement it
        pass

    @property
    def mol_obj(self):
        return self._mol_obj

    @mol_obj.setter
    def mol_obj(self, mol_obj):
        if isinstance(mol_obj, PybelMol):
            mol_obj = mol_obj.OBMol

        if not isinstance(mol_obj, RDMol) and not isinstance(mol_obj, OBMol):
            logger.exception("Objects of type '%s' are not currently accepted." % mol_obj.__class__)
            raise MoleculeObjectTypeError("Objects of type '%s' are not currently accepted." % mol_obj.__class__)

        self._mol_obj = mol_obj

    def get_atoms(self, wrapped=True):

        atoms = []
        if self.is_rdkit_obj():
            atoms = self._mol_obj.GetAtoms()
        elif self.is_openbabel_obj():
            atoms = OBMolAtomIter(self._mol_obj)

        if atoms and wrapped:
            return [AtomWrapper(atm) for atm in atoms]
        return atoms

    def get_bonds(self, wrapped=True):

        bonds = []
        if self.is_rdkit_obj():
            bonds = self._mol_obj.GetBonds()
        elif self.is_openbabel_obj():
            bonds = OBMolBondIter(self._mol_obj)

        if bonds and wrapped:
            return [BondWrapper(bond) for bond in bonds]
        return bonds

    def get_num_heavy_atoms(self):
        if self.is_rdkit_obj():
            return self._mol_obj.GetNumHeavyAtoms()
        elif self.is_openbabel_obj():
            return self._mol_obj.NumHvyAtoms()

    def get_atom_coord_by_id(self, atm_id):
        if self.is_rdkit_obj():
            return list(self._mol_obj.GetConformer().GetAtomPosition(atm_id))
        elif self.is_openbabel_obj():
            atm = self._mol_obj.GetAtomById(atm_id)
            return [atm.GetX(), atm.GetY(), atm.GetZ()]

    def to_smiles(self):
        if self.is_rdkit_obj():
            return MolToSmiles(self._mol_obj)
        elif self.is_openbabel_obj():
            return PybelMol(self._mol_obj).write("smi").split("\t")[0]

    def to_pdb_block(self):
        if self.is_rdkit_obj():
            return MolToPDBBlock(self._mol_obj)
        elif self.is_openbabel_obj():
            return PybelMol(self._mol_obj).write("pdb")

    def unwrap(self):
        return self._mol_obj

    def is_rdkit_obj(self):
        if isinstance(self._mol_obj, RDMol):
            return True
        return False

    def is_openbabel_obj(self):
        if isinstance(self._mol_obj, OBMol):
            return True
        elif isinstance(self._mol_obj, PybelMol):
            return True
        return False

    def is_pybel_obj(self):
        if isinstance(self._mol_obj, PybelMol):
            return True
        return False

    def __getattr__(self, attr):
        return getattr(self._mol_obj, attr)
