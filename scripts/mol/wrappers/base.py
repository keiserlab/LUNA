from util.exceptions import AtomObjectTypeError, MoleculeObjectTypeError
from rdkit.Chem import Atom as RDAtom
from rdkit.Chem import Mol as RDMol
from rdkit.Chem import Bond as RDBond
from rdkit.Chem import MolToSmiles
from openbabel import OBMol, OBAtom, OBBond, OBMolAtomIter, OBAtomAtomIter, OBAtomBondIter, OBMolBondIter
from openbabel import etab
from pybel import Molecule as PybelMol

import logging

logger = logging.getLogger()


class AtomWrapper:

    def __init__(self, atm_obj):
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

    def get_valence(self):
        if self.is_rdkit_obj():
            return self._atm_obj.GetExplicitValence() + self._atm_obj.GetImplicitValence()
        elif self.is_openbabel_obj():
            bonds = [b.GetBondOrder() for b in OBAtomBondIter(self._atm_obj)]
            return sum(bonds) + self._atm_obj.ImplicitHydrogenCount()

    def get_implicit_valence(self):
        # Both RDKit and Openbabel have the same function name.
        return self._atm_obj.GetImplicitValence()

    def get_neighbors(self):
        if self.is_rdkit_obj():
            return self._atm_obj.GetNeighbors()
        elif self.is_openbabel_obj():
            return OBAtomAtomIter(self._atm_obj)

    def get_degree(self):
        if self.is_rdkit_obj():
            # In RDKit, GetDegree() returns the degree of the atom.
            return self._atm_obj.GetDegree()
        elif self.is_openbabel_obj():
            # In OpenBabel, GetValence() returns the current number of explicit connections (degree).
            return self._atm_obj.GetValence()

    def get_implicit_h_count(self):
        if self.is_rdkit_obj():
            return self._atm_obj.GetNumImplicitHs()
        elif self.is_openbabel_obj():
            return self._atm_obj.ImplicitHydrogenCount()

    def get_explicit_h_count(self):
        if self.is_rdkit_obj():
            return self._atm_obj.GetNumExplicitHs()
        elif self.is_openbabel_obj():
            return self._atm_obj.ExplicitHydrogenCount()

    def get_charge(self):
        # Both RDKit and Openbabel have the same function name.
        return self._atm_obj.GetFormalCharge()

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

    def get_begin_atom(self, wrapped=False):
        # Both RDKit and Openbabel have the same function name.
        if wrapped:
            return AtomWrapper(self._bond_obj.GetBeginAtom())
        return self._bond_obj.GetBeginAtom()

    def get_end_atom(self, wrapped=False):
        # Both RDKit and Openbabel have the same function name.
        if wrapped:
            return AtomWrapper(self._bond_obj.GetEndAtom())
        return self._bond_obj.GetEndAtom()

    def get_bond_type(self):
        if self.is_rdkit_obj():
            return self._bond_obj.GetBondType()
        elif self.is_openbabel_obj():
            return self._bond_obj.GetBondOrder()

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
        if isinstance(mol_obj, PybelMol):
            mol_obj = mol_obj.OBMol

        if not isinstance(mol_obj, RDMol) and not isinstance(mol_obj, OBMol):
            logger.exception("Objects of type '%s' are not currently accepted." % mol_obj.__class__)
            raise MoleculeObjectTypeError("Objects of type '%s' are not currently accepted." % mol_obj.__class__)

        self._mol_obj = mol_obj

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

    def get_atoms(self, wrapped=False):
        if self.is_rdkit_obj():
            atoms = self._mol_obj.GetAtoms()
        elif self.is_openbabel_obj():
            atoms = OBMolAtomIter(self._mol_obj)

        if atoms:
            if wrapped:
                return [AtomWrapper(atm) for atm in atoms]
            return atoms
        return []

    def get_bonds(self, wrapped=False):
        if self.is_rdkit_obj():
            bonds = self._mol_obj.GetBonds()
        elif self.is_openbabel_obj():
            bonds = OBMolBondIter(self._mol_obj)

        if bonds:
            if wrapped:
                return [BondWrapper(bond) for bond in bonds]
            return bonds
        return []

    def get_num_heavy_atoms(self):
        if self.is_rdkit_obj():
            return self._mol_obj.GetNumHeavyAtoms()
        elif self.is_openbabel_obj():
            return self._mol_obj.NumHvyAtoms()

    def get_atom_coord_by_id(self, atom_id):
        if self.is_rdkit_obj():
            return list(self._mol_obj.GetConformer().GetAtomPosition(atom_id))
        elif self.is_openbabel_obj():
            atm = self._mol_obj.GetAtomById(atom_id)
            return [atm.GetX(), atm.GetY(), atm.GetZ()]

    def to_smiles(self):
        if self.is_rdkit_obj():
            return MolToSmiles(self._mol_obj)
        elif self.is_openbabel_obj():
            return PybelMol(self._mol_obj).write("smi").split("\t")[0]

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
