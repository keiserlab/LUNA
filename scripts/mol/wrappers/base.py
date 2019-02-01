from util.exceptions import (AtomObjectTypeError, MoleculeObjectTypeError)
from rdkit.Chem import Atom as RDAtom
from rdkit.Chem import Mol as RDMol
from openbabel import (OBMol, OBAtom, OBMolAtomIter, OBAtomAtomIter, OBAtomBondIter)
from openbabel import etab
from pybel import Molecule as PybelMol

import logging


logger = logging.getLogger(__name__)


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

    def get_id(self):
        # Both RDKit and Openbabel have the same function name.
        return self.atm_obj.GetIdx()

    def get_atomic_num(self):
        # Both RDKit and Openbabel have the same function name.
        return self.atm_obj.GetAtomicNum()

    def get_symbol(self):
        if self.is_rdkit_obj():
            return self.atm_obj.GetSymbol()
        elif self.is_openbabel_obj():
            return etab.GetSymbol(self.atm_obj.GetAtomicNum())

    def get_valence(self):
        if self.is_rdkit_obj():
            return self.atm_obj.GetExplicitValence() + self.atm_obj.GetImplicitValence()
        elif self.is_openbabel_obj():
            bonds = [b.GetBondOrder() for b in OBAtomBondIter(self.atm_obj)]
            return sum(bonds) + self.atm_obj.ImplicitHydrogenCount()

    def get_implicit_valence(self):
        # Both RDKit and Openbabel have the same function name.
        return self.atm_obj.GetImplicitValence()

    def get_neighbors(self):
        if self.is_rdkit_obj():
            return self.atm_obj.GetNeighbors()
        elif self.is_openbabel_obj():
            return OBAtomAtomIter(self.atm_obj)

    def get_degree(self):
        if self.is_rdkit_obj():
            # In RDKit, GetDegree() returns the degree of the atom.
            return self.atm_obj.GetDegree()
        elif self.is_openbabel_obj():
            # In OpenBabel, GetValence() returns the current number of explicit connections (degree).
            return self.atm_obj.GetValence()

    def get_implicit_h_count(self):
        if self.is_rdkit_obj():
            return self.atm_obj.GetNumImplicitHs()
        elif self.is_openbabel_obj():
            return self.atm_obj.ImplicitHydrogenCount()

    def get_explicit_h_count(self):
        if self.is_rdkit_obj():
            return self.atm_obj.GetNumExplicitHs()
        elif self.is_openbabel_obj():
            return self.atm_obj.ExplicitHydrogenCount()

    def get_charge(self):
        # Both RDKit and Openbabel have the same function name.
        return self.atm_obj.GetFormalCharge()

    def set_charge(self, charge):
        # Both RDKit and Openbabel have the same function name.
        self.atm_obj.SetFormalCharge(charge)

    def set_implicit_valence(self, valence):
        if self.is_rdkit_obj():
            logger.warning("RDKit does not provide a setter for implicit valence.")
        elif self.is_openbabel_obj():
            self.atm_obj.SetImplicitValence(valence)

    def is_rdkit_obj(self):
        if isinstance(self.atm_obj, RDAtom):
            return True
        return False

    def is_openbabel_obj(self):
        if isinstance(self.atm_obj, OBAtom):
            return True
        return False

    def __getattr__(self, attr):
        return getattr(self.atm_obj, attr)


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

    def get_atoms(self):
        if self.is_rdkit_obj():
            return self.mol_obj.GetAtoms()
        elif self.is_openbabel_obj():
            return OBMolAtomIter(self.mol_obj)

    def is_rdkit_obj(self):
        if isinstance(self.mol_obj, RDMol):
            return True
        return False

    def is_openbabel_obj(self):
        if isinstance(self.mol_obj, OBMol):
            return True
        elif isinstance(self.mol_obj, PybelMol):
            return True
        return False

    def is_pybel_obj(self):
        if isinstance(self.mol_obj, PybelMol):
            return True
        return False
