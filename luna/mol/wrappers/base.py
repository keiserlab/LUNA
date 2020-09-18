from enum import Enum, unique
from rdkit.Chem import Atom as RDAtom
from rdkit.Chem import Mol as RDMol
from rdkit.Chem import Bond as RDBond
from rdkit.Chem import BondType as RDBondType
from rdkit.Chem import MolFromSmiles, MolToSmiles, MolToPDBBlock, MolToMolBlock, GetPeriodicTable
from openbabel import openbabel as ob
from openbabel.pybel import Molecule as PybelMol
from openbabel.pybel import readstring, readfile, readstring

from luna.util.file import get_file_format, get_filename
from luna.mol.wrappers.rdkit import read_multimol_file, read_mol_from_file, new_mol_from_block
from luna.util.exceptions import (AtomObjectTypeError, MoleculeObjectTypeError, IllegalArgumentError, MoleculeObjectError,
                                  MoleculeNotFoundError)

import logging

logger = logging.getLogger()


@unique
class BondType(Enum):
    # Same values and bond types available at RDKit.
    UNSPECIFIED = 0
    SINGLE = 1
    DOUBLE = 2
    TRIPLE = 3
    QUADRUPLE = 4
    QUINTUPLE = 5
    HEXTUPLE = 6
    ONEANDAHALF = 7
    TWOANDAHALF = 8
    THREEANDAHALF = 9
    FOURANDAHALF = 10
    FIVEANDAHALF = 11
    AROMATIC = 12
    IONIC = 13
    HYDROGEN = 14
    THREECENTER = 15
    DATIVEONE = 16
    DATIVE = 17
    DATIVEL = 18
    DATIVER = 19
    OTHER = 20
    ZERO = 21


class OBBondType(Enum):
    # Same values and bond types available at OpenBabel.
    SINGLE = 1
    DOUBLE = 2
    TRIPLE = 3
    AROMATIC = 5


class AtomWrapper:

    def __init__(self, atm_obj, mol_obj=None):
        if isinstance(atm_obj, self.__class__):
            atm_obj = atm_obj.unwrap()

        if not isinstance(atm_obj, RDAtom) and not isinstance(atm_obj, ob.OBAtom):
            logger.exception("Objects of type '%s' are not currently accepted." % atm_obj.__class__)
            raise AtomObjectTypeError("Objects of type '%s' are not currently accepted." % atm_obj.__class__)

        self._atm_obj = atm_obj
        self._parent = None

        if mol_obj is None:
            mol_obj = self.get_parent()
        self._parent = MolWrapper(mol_obj)

    @property
    def atm_obj(self):
        return self._atm_obj

    @atm_obj.setter
    def atm_obj(self, atm_obj):
        if not isinstance(atm_obj, RDAtom) and not isinstance(atm_obj, ob.OBAtom):
            logger.exception("Objects of type '%s' are not currently accepted." % atm_obj.__class__)
            raise AtomObjectTypeError("Objects of type '%s' are not currently accepted." % atm_obj.__class__)
        else:
            self._atm_obj = atm_obj

    @property
    def parent(self):
        return self.get_parent()

    @parent.setter
    def parent(self, mol_obj):
        self._parent = MolWrapper(mol_obj)

    def get_idx(self):
        # Both RDKit and Open Babel have the same function name.
        return self._atm_obj.GetIdx()

    def get_id(self):
        if self.is_rdkit_obj():
            return self._atm_obj.GetIdx()
        elif self.is_openbabel_obj():
            return self._atm_obj.GetId()

    def get_parent(self):
        if self._parent is None:
            if self.is_rdkit_obj():
                mol_obj = self._atm_obj.GetOwningMol()
            elif self.is_openbabel_obj():
                mol_obj = self._atm_obj.GetParent()

            self._parent = mol_obj
        return self._parent

    def get_atomic_num(self):
        # Both RDKit and Open Babel have the same function name.
        return self._atm_obj.GetAtomicNum()

    def get_symbol(self):
        if self.is_rdkit_obj():
            return self._atm_obj.GetSymbol()
        elif self.is_openbabel_obj():
            return ob.GetSymbol(self._atm_obj.GetAtomicNum())

    def get_charge(self):
        # Both RDKit and Open Babel have the same function name.
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

    def get_neighbors(self, wrapped=True):
        atoms = []
        if self.is_rdkit_obj():
            atoms = self._atm_obj.GetNeighbors()
        elif self.is_openbabel_obj():
            atoms = ob.OBAtomAtomIter(self._atm_obj)

        if atoms and wrapped:
            return [AtomWrapper(atom) for atom in atoms]
        return atoms

    def get_degree(self):
        if self.is_rdkit_obj():
            # In RDKit, GetDegree() returns the degree of the atom.
            return self._atm_obj.GetTotalDegree()
        elif self.is_openbabel_obj():
            return self._atm_obj.GetTotalDegree()

    def get_valence(self):
        if self.is_rdkit_obj():
            return self._atm_obj.GetExplicitValence() + self._atm_obj.GetImplicitValence()
        elif self.is_openbabel_obj():
            return self._atm_obj.GetTotalValence()

            bonds = [b.GetBondOrder() for b in ob.OBAtomBondIter(self._atm_obj)]
            return sum(bonds) + self._atm_obj.ImplicitHydrogenCount()

    def get_h_count(self):
        if self.is_rdkit_obj():
            return self._atm_obj.GetTotalNumHs(includeNeighbors=True)
        elif self.is_openbabel_obj():
            return self._atm_obj.GetImplicitHCount() + self._atm_obj.ExplicitHydrogenCount()

    def get_bonds(self, wrapped=True):
        bonds = []
        if self.is_rdkit_obj():
            bonds = self._atm_obj.GetBonds()
        elif self.is_openbabel_obj():
            bonds = ob.OBAtomBondIter(self._atm_obj)

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
        # Both RDKit and Open Babel have the same function name.
        return self._atm_obj.IsInRing()

    def is_aromatic(self):
        if self.is_rdkit_obj():
            return self._atm_obj.GetIsAromatic()
        elif self.is_openbabel_obj():
            return self._atm_obj.IsAromatic()

    def has_bond_type(self, bond_type):
        if isinstance(bond_type, BondType):
            return any([bond_obj.get_bond_type() == bond_type for bond_obj in self.get_bonds()])
        else:
            raise IllegalArgumentError("The informed bond type must be an instance of '%s'." % BondType)

    def has_only_bond_type(self, bond_type):
        if isinstance(bond_type, BondType):
            return all([bond_obj.get_bond_type() == bond_type for bond_obj in self.get_bonds()])
        else:
            raise IllegalArgumentError("The informed bond type must be an instance of '%s'." % BondType)

    def set_charge(self, charge):
        # Both RDKit and Open Babel have the same function name.
        self._atm_obj.SetFormalCharge(charge)

    def set_as_aromatic(self, is_aromatic):
        if self.is_rdkit_obj():
            self._atm_obj.SetIsAromatic(is_aromatic)
        elif self.is_openbabel_obj():
            self._atm_obj.SetAromatic(is_aromatic)

    def set_in_ring(self, in_ring):
        if self.is_rdkit_obj():
            # TODO
            pass
        elif self.is_openbabel_obj():
            self._atm_obj.SetInRing(in_ring)

    def matches_smarts(self, smarts):

        if self.is_rdkit_obj():
            # TODO: Implement
            pass

        elif self.is_openbabel_obj():
            ob_smart = ob.OBSmartsPattern()
            ob_smart.Init(smarts)

            if ob_smart.Match(self.parent):
                for match in ob_smart.GetMapList():
                    if match[0] == self.get_idx():
                        return True
            return False

    def unwrap(self):
        return self._atm_obj

    def is_rdkit_obj(self):
        if isinstance(self._atm_obj, RDAtom):
            return True
        return False

    def is_openbabel_obj(self):
        if isinstance(self._atm_obj, ob.OBAtom):
            return True
        return False

    def __getattr__(self, attr):
        return getattr(self._atm_obj, attr)


class BondWrapper:
    def __init__(self, bond_obj):
        if isinstance(bond_obj, self.__class__):
            bond_obj = bond_obj.unwrap()

        if not isinstance(bond_obj, RDBond) and not isinstance(bond_obj, ob.OBBond):
            logger.exception("Objects of type '%s' are not currently accepted." % bond_obj.__class__)
            raise MoleculeObjectTypeError("Objects of type '%s' are not currently accepted." % bond_obj.__class__)

        self._bond_obj = bond_obj

    @property
    def bond_obj(self):
        return self._bond_obj

    @bond_obj.setter
    def bond_obj(self, bond_obj):
        if not isinstance(bond_obj, RDBond) and not isinstance(bond_obj, ob.OBBond):
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
        # Both RDKit and Open Babel have the same function name.
        if wrapped:
            return AtomWrapper(self._bond_obj.GetBeginAtom())
        return self._bond_obj.GetBeginAtom()

    def get_end_atom(self, wrapped=True):
        # Both RDKit and Open Babel have the same function name.
        if wrapped:
            return AtomWrapper(self._bond_obj.GetEndAtom())
        return self._bond_obj.GetEndAtom()

    def get_bond_type(self):
        if self.is_rdkit_obj():
            return BondType[self._bond_obj.GetBondType().name]
        elif self.is_openbabel_obj():
            # Map Open Babel bonds to BondType.
            return BondType[OBBondType(self._bond_obj.GetBondOrder()).name]

    def is_aromatic(self):
        if self.is_rdkit_obj():
            return self._bond_obj.GetIsAromatic()
        elif self.is_openbabel_obj():
            return self._bond_obj.IsAromatic()

    def set_bond_type(self, bond_type):

        # Convert BondType to a valid type accepted by Open Babel or RDKit.
        if isinstance(bond_type, BondType):
            try:
                if self.is_rdkit_obj():
                    bond_type = RDBondType.values[bond_type.value]
                elif self.is_openbabel_obj():
                    bond_type = OBBondType[bond_type.name].value
            except KeyError as e:
                logger.exception(e)

                tool = "Open Babel" if self.is_openbabel_obj() else "RDKit"
                raise KeyError("The bond type '%s' is not a valid %s bond type." % (bond_type.name, tool))
        else:
            raise IllegalArgumentError("The informed bond type must be an instance of '%s'." % BondType)

        if self.is_rdkit_obj():
            try:
                self._bond_obj.SetBondType(bond_type)
            except Exception as e:
                logger.exception(e)
                raise
        elif self.is_openbabel_obj():
            self._bond_obj.SetBondOrder(bond_type)

    def set_as_aromatic(self, is_aromatic):
        if self.is_rdkit_obj():
            self._bond_obj.SetIsAromatic(is_aromatic)
        elif self.is_openbabel_obj():
            self._bond_obj.SetAromatic(is_aromatic)

    def unwrap(self):
        return self._bond_obj

    def is_rdkit_obj(self):
        if isinstance(self._bond_obj, RDBond):
            return True
        return False

    def is_openbabel_obj(self):
        if isinstance(self._bond_obj, ob.OBBond):
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

        if not isinstance(mol_obj, RDMol) and not isinstance(mol_obj, ob.OBMol):
            logger.exception("Objects of type '%s' are not currently accepted." % mol_obj.__class__)
            raise MoleculeObjectTypeError("Objects of type '%s' are not currently accepted." % mol_obj.__class__)

        self._mol_obj = mol_obj

    @classmethod
    def from_mol_file(cls, mol_file, mol_id=None, mol_file_ext=None, mol_obj_type="rdkit"):
        pass

        # TODO: Test below code:

        """
            mol_file_ext = mol_file_ext or get_file_format(mol_file)

            tool = "Open Babel" if mol_obj_type == "openbabel" else "RDKit"

            try:
                if mol_obj_type == "openbabel":
                    mols = readfile(mol_file_ext, mol_file)
                    # If it is a multimol file, then we need to loop over the molecules to find the target one.
                    # Note that in this case, the ids must match.
                    if mol_id is not None:
                        for ob_mol in mols:
                            if mol_id == get_filename(ob_mol.OBMol.GetTitle()):
                                mol_obj = ob_mol
                                break
                    else:
                        mol_obj = mols.__next__()
                else:
                    if mol_file_ext == "pdb":
                        mol_obj = read_mol_from_file(mol_file, mol_format=mol_file_ext, removeHs=False)
                    else:
                        # If 'targets' is None, then the entire Mol file will be read.
                        targets = None
                        # If it is a multimol file than loop through it until the informed molecule (by its mol_id) is found.
                        if mol_id:
                            targets = [mol_id]

                        for rdk_mol, mol_id in read_multimol_file(mol_file, mol_format=mol_file_ext, targets=targets, removeHs=False):
                            # It returns None if the molecule parsing generated errors.
                            mol_obj = rdk_mol
                            break
            except Exception as e:
                logger.exception(e)
                raise MoleculeObjectError("An error occurred while parsing the molecular file '%s' with %s." % (mol_file, tool))

            if mol_obj is None:
                if mol_id:
                    raise MoleculeNotFoundError("The ligand '%s' was not found in the molecular file '%s' or some errors "
                                                "were found while parsing it with %s." % (mol_id, mol_file, tool))
                else:
                    raise MoleculeNotFoundError("Some errors were found while parsing the molecular file '%s' with %s." % (mol_file, tool))

            return cls(mol_obj)
        """

    @classmethod
    def from_mol_block(cls, block, mol_format, mol_obj_type="rdkit"):
        if mol_obj_type == "rdkit":
            return cls(new_mol_from_block(block, mol_format))
        elif mol_obj_type == "openbabel":
            return cls(readstring(mol_format, block))

    @classmethod
    def from_smiles(cls, smiles, mol_obj_type="rdkit"):
        if mol_obj_type == "rdkit":
            mol = MolFromSmiles(smiles)
        elif mol_obj_type == "openbabel":
            mol = readstring("smi", smiles)

        # Raise an error if the molecule could not be created.
        if mol is None:
            raise MoleculeObjectError("It could not create a molecule from the provided SMILES '%s'." % smiles)
        return cls(mol)

    @property
    def mol_obj(self):
        return self._mol_obj

    @mol_obj.setter
    def mol_obj(self, mol_obj):
        if isinstance(mol_obj, PybelMol):
            mol_obj = mol_obj.OBMol

        if not isinstance(mol_obj, RDMol) and not isinstance(mol_obj, ob.OBMol):
            logger.exception("Objects of type '%s' are not currently accepted." % mol_obj.__class__)
            raise MoleculeObjectTypeError("Objects of type '%s' are not currently accepted." % mol_obj.__class__)

        self._mol_obj = mol_obj

    def as_rdkit(self):
        if self.is_rdkit_obj():
            return self._mol_obj
        elif self.is_openbabel_obj():
            new_mol = MolWrapper.from_smiles(self.to_smiles(), mol_obj_type="rdkit")
            new_mol.set_name(self.get_name())
            return new_mol.unwrap()

    def as_openbabel(self):
        if self.is_openbabel_obj():
            return self._mol_obj
        elif self.is_rdkit_obj():
            new_mol = MolWrapper.from_smiles(self.to_smiles(), mol_obj_type="openbabel")
            new_mol.set_name(self.get_name())
            return new_mol.unwrap()

    def get_name(self):
        if self.is_rdkit_obj():
            if self._mol_obj.HasProp("_Name"):
                return self._mol_obj.GetProp("_Name")
            return ""
        elif self.is_openbabel_obj():
            return self._mol_obj.GetTitle()

    def get_atoms(self, wrapped=True):
        atoms = []
        if self.is_rdkit_obj():
            atoms = self._mol_obj.GetAtoms()
        elif self.is_openbabel_obj():
            atoms = ob.OBMolAtomIter(self._mol_obj)

        if atoms and wrapped:
            return [AtomWrapper(atm, self) for atm in atoms]
        return atoms

    def get_bonds(self, wrapped=True):
        bonds = []
        if self.is_rdkit_obj():
            bonds = self._mol_obj.GetBonds()
        elif self.is_openbabel_obj():
            bonds = ob.OBMolBondIter(self._mol_obj)

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

    def get_obj_type(self):
        if self.is_rdkit_obj():
            return "rdkit"
        elif self.is_openbabel_obj():
            return "openbabel"

    def has_name(self):
        if self.get_name():
            return True
        return False

    def set_name(self, name):
        if self.is_rdkit_obj():
            self._mol_obj.SetProp("_Name", name)
        elif self.is_openbabel_obj():
            self._mol_obj.SetTitle(name)

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

    def to_mol_block(self):
        if self.is_rdkit_obj():
            return MolToMolBlock(self._mol_obj)
        elif self.is_openbabel_obj():
            return PybelMol(self._mol_obj).write("mol")

    def unwrap(self):
        return self._mol_obj

    def is_rdkit_obj(self):
        if isinstance(self._mol_obj, RDMol):
            return True
        return False

    def is_openbabel_obj(self):
        if isinstance(self._mol_obj, ob.OBMol):
            return True
        elif isinstance(self._mol_obj, PybelMol):
            return True
        return False

    def is_pybel_obj(self):
        if isinstance(self._mol_obj, PybelMol):
            return True
        return False

    def __getattr__(self, attr):
        if self.mol_obj is None:
            return None
        return getattr(self._mol_obj, attr)

    def __getstate__(self):
        # Creates a copy of the class' dictionary in case we need to modify the molecular object to pickle it.
        my_dict = self.__dict__.copy()
        if self.is_openbabel_obj():
            my_dict["_mol_block"] = self.to_mol_block()
            my_dict["_mol_obj_type"] = self.get_obj_type()
            my_dict["_mol_obj"] = None
        return my_dict

    def __setstate__(self, state):
        if "_mol_obj_type" in state and "_mol_block" in state:
            state["_mol_obj"] = self.from_mol_block(state["_mol_block"], "mol", "openbabel").unwrap()
            del state["_mol_obj_type"]
            del state["_mol_block"]

        self.__dict__.update(state)
