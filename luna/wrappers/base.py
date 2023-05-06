from enum import Enum, unique
from rdkit.Chem import Atom as RDAtom
from rdkit.Chem import Mol as RDMol
from rdkit.Chem import Bond as RDBond
from rdkit.Chem import BondType as RDBondType
from rdkit.Chem import (MolFromSmiles, MolToSmiles, MolToPDBBlock,
                        MolToMolBlock, GetPeriodicTable, GetFormalCharge)
from rdkit.Chem import rdFMCS

from openbabel import openbabel as ob
from openbabel.pybel import Molecule as PybelMol
from openbabel.pybel import readstring, readfile

from luna.wrappers.rdkit import new_mol_from_block, read_mol_from_file
from luna.util.math import euclidean_distance
from luna.util.exceptions import (AtomObjectTypeError, BondObjectTypeError,
                                  IllegalArgumentError, MoleculeObjectError,
                                  MoleculeObjectTypeError)

import logging
logger = logging.getLogger()


@unique
class BondType(Enum):
    """An enumeration of bond types available at RDKit."""

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
    """An enumeration of bond types available at Open Babel."""

    # Same values and bond types available at OpenBabel.
    SINGLE = 1
    DOUBLE = 2
    TRIPLE = 3
    AROMATIC = 5


class AtomWrapper:
    """This class provides util functions to access atomic properties and
    other information from RDKit and Open Babel objects.

    Parameters
    ----------
    atm_obj : :class:`AtomWrapper`, :class:`rdkit.Chem.rdchem.Atom`, \
                    or :class:`openbabel.OBAtom`
        An atom to wrap.
    mol_obj : `MolWrapper`, :class:`rdkit.Chem.rdchem.Mol`, \
                    :class:`openbabel.pybel.Molecule`, or None
        The molecule that contains the atom ``atom``. If None, the molecule is
        recovered directly from the atom object.

    Raises
    ------
    AtomObjectTypeError
        If the atom object is not an instance
            of `MolWrapper`, :class:`rdkit.Chem.rdchem.Mol`, or
            :class:`openbabel.pybel.Molecule`.
    """

    def __init__(self, atm_obj, mol_obj=None):
        if isinstance(atm_obj, self.__class__):
            atm_obj = atm_obj.unwrap()

        if (not isinstance(atm_obj, RDAtom)
                and not isinstance(atm_obj, ob.OBAtom)):
            error_msg = ("Objects of type '%s' are not currently accepted."
                         % atm_obj.__class__)
            logger.exception(error_msg)
            raise AtomObjectTypeError(error_msg)

        self._atm_obj = atm_obj
        self._parent = None

        if mol_obj is None:
            mol_obj = self.get_parent()
        self._parent = MolWrapper(mol_obj)

    @property
    def atm_obj(self):
        """:class:`rdkit.Chem.rdchem.Atom` or :class:`openbabel.OBAtom`: \
                The wrapped atom object."""
        return self._atm_obj

    @atm_obj.setter
    def atm_obj(self, atm_obj):
        if (not isinstance(atm_obj, RDAtom)
                and not isinstance(atm_obj, ob.OBAtom)):
            error_msg = ("Objects of type '%s' are not currently accepted."
                         % atm_obj.__class__)
            logger.exception(error_msg)
            raise AtomObjectTypeError(error_msg)
        else:
            self._atm_obj = atm_obj

    @property
    def parent(self):
        """`MolWrapper`, :class:`rdkit.Chem.rdchem.Mol`, \
                or :class:`openbabel.pybel.Molecule`: \
            The molecule that contains this atom."""
        return self.get_parent()

    @parent.setter
    def parent(self, mol_obj):
        self._parent = MolWrapper(mol_obj)

    def get_idx(self):
        """Get this atom's internal index within a molecule.

        Returns
        -------
         : int
        """

        # Both RDKit and Open Babel have the same function name.
        return self._atm_obj.GetIdx()

    def get_id(self):
        """Get this atom's unique id.

        When using RDKit, :py:meth:`get_id` and :py:meth:`get_idx` returns
        the same value.

        Returns
        -------
         : int
        """
        if self.is_rdkit_obj():
            return self._atm_obj.GetIdx()
        elif self.is_openbabel_obj():
            return self._atm_obj.GetId()

    def get_parent(self, wrapped=True):
        """Get the molecule that contains this atom.

        Parameters
        ----------
        wrapped : bool
            If True, wrap the molecule with
            :class:`~luna.wrappers.base.MolWrapper`.

        Returns
        -------
         : `MolWrapper`, :class:`rdkit.Chem.rdchem.Mol`, \
                            or :class:`openbabel.pybel.Molecule`
        """
        if self._parent is None:
            if self.is_rdkit_obj():
                mol_obj = self._atm_obj.GetOwningMol()
            elif self.is_openbabel_obj():
                mol_obj = self._atm_obj.GetParent()

            self._parent = mol_obj

        if wrapped:
            return MolWrapper(self._parent)
        return self._parent

    def get_atomic_num(self):
        """Get this atom's atomic number.

        Returns
        -------
         : int
        """

        # Both RDKit and Open Babel have the same function name.
        return self._atm_obj.GetAtomicNum()

    def get_symbol(self):
        """Get the element symbol of this atom.

        Returns
        -------
         : str
        """
        if self.is_rdkit_obj():
            return self._atm_obj.GetSymbol()
        elif self.is_openbabel_obj():
            return ob.GetSymbol(self._atm_obj.GetAtomicNum())

    def get_charge(self):
        """Get this atom's formal charge.

        Returns
        -------
         : int
        """

        # Both RDKit and Open Babel have the same function name.
        return self._atm_obj.GetFormalCharge()

    def get_isotope(self):
        """Get this atom's isotope number.

        Returns
        -------
         : int
        """
        if self.is_rdkit_obj():
            return int(round(self._atm_obj.GetMass()))
        elif self.is_openbabel_obj():
            return int(round(self._atm_obj.GetExactMass()))

    def get_mass(self):
        """Get this atom's exact atomic mass, which can vary given the isotope.

        Returns
        -------
         : float
        """
        if self.is_rdkit_obj():
            return self._atm_obj.GetMass()
        elif self.is_openbabel_obj():
            return self._atm_obj.GetExactMass()

    def get_atomic_mass(self):
        """Get this atom's atomic mass given by standard IUPAC average
        molar mass.

        Returns
        -------
         : float
        """
        if self.is_rdkit_obj():
            return GetPeriodicTable().GetAtomicWeight(self.get_atomic_num())
        elif self.is_openbabel_obj():
            return self._atm_obj.GetAtomicMass()

    def get_neighbors_number(self, only_heavy_atoms=False):
        """Get the number of atoms bound to this atom.

        Parameters
        ----------
        only_heavy_atoms : bool
            If True, count only heavy atoms.

        Returns
        -------
         : int
        """

        # Subtract the number of hydrogens from the degree if only heavy atoms
        # must be considered.
        penalty = 0 if only_heavy_atoms is False else self.get_h_count()
        return self.get_degree() - penalty

    def get_neighbors(self, wrapped=True):
        """Get all atoms that are bound to this atom.

        Parameters
        ----------
        wrapped : bool
            If True, wrap each atom with
            :class:`~luna.wrappers.base.AtomWrapper`.

        Returns
        -------
         : iterable of `AtomWrapper`, :py:class:`rdkit.Chem.rdchem.Atom`, \
                or :class:`openbabel.OBAtom`
        """
        atoms = []
        if self.is_rdkit_obj():
            atoms = self._atm_obj.GetNeighbors()
        elif self.is_openbabel_obj():
            atoms = ob.OBAtomAtomIter(self._atm_obj)

        if atoms and wrapped:
            return [AtomWrapper(atom) for atom in atoms]
        return atoms

    def get_degree(self):
        """Get this atom's total degree.

        Returns
        -------
         : int
        """
        if self.is_rdkit_obj():
            return self._atm_obj.GetTotalDegree()
        elif self.is_openbabel_obj():
            return self._atm_obj.GetTotalDegree()

    def get_valence(self):
        """Get this atom's total valence (implicit and explicit)."""
        if self.is_rdkit_obj():
            valence = (self._atm_obj.GetExplicitValence()
                       + self._atm_obj.GetImplicitValence())
            return valence
        elif self.is_openbabel_obj():
            return self._atm_obj.GetTotalValence()

    def get_h_count(self):
        """Get the total number of hydrogens (implicit and explicit) bound to
        this atom.

        Returns
        -------
         : int
        """
        if self.is_rdkit_obj():
            return self._atm_obj.GetTotalNumHs(includeNeighbors=True)
        elif self.is_openbabel_obj():
            h_count = (self._atm_obj.GetImplicitHCount()
                       + self._atm_obj.ExplicitHydrogenCount())
            return h_count

    def get_bonds(self, wrapped=True):
        """Get this atom’s bonds.

        Parameters
        ----------
        wrapped : bool
            If True, wrap each bond with
            :class:`~luna.wrappers.base.BondWrapper`.

        Returns
        -------
         : iterable of `BondWrapper`, :class:`rdkit.Chem.rdchem.Bond`,
            or :class:`openbabel.OBBond`
        """
        bonds = []
        if self.is_rdkit_obj():
            bonds = self._atm_obj.GetBonds()
        elif self.is_openbabel_obj():
            bonds = ob.OBAtomBondIter(self._atm_obj)

        if bonds and wrapped:
            return [BondWrapper(bond) for bond in bonds]
        return bonds

    def get_coord(self):
        return self.parent.get_atom_coord_by_id(self.get_id())

    def get_atomic_invariants(self):
        """Get the atomic invariants of this atom.

        Atomic invariants are derived from ECFP [1]_ and E3FP [2]_ and
        consists of seven fields:

            * Number of neighboring heavy atoms;
            * Valence - Number of hydrogens;
            * Atomic number;
            * Isotope number;
            * Formal charge;
            * Number of hydrogens;
            * If the atom belongs to a ring or not.

        Returns
        -------
         : list
        """
        return [self.get_neighbors_number(only_heavy_atoms=True),
                (self.get_valence() - self.get_h_count()),
                self.get_atomic_num(),
                self.get_isotope(),
                self.get_charge(),
                self.get_h_count(),
                int(self.is_in_ring())]

    def is_in_ring(self):
        """Check if this atom is in a ring."""

        # Both RDKit and Open Babel have the same function name.
        return self._atm_obj.IsInRing()

    def is_aromatic(self):
        """Check if this atom is aromatic."""

        if self.is_rdkit_obj():
            return self._atm_obj.GetIsAromatic()
        elif self.is_openbabel_obj():
            return self._atm_obj.IsAromatic()

    def has_bond_type(self, bond_type):
        """Check if this atom has a bond of type ``bond_type``.

        Parameters
        ----------
        bond_type : `BondType`

        Raises
        ------
        IllegalArgumentError
            If the informed bond type is not an instance of `BondType`.
        """
        if isinstance(bond_type, BondType):
            return any([bond_obj.get_bond_type() == bond_type
                        for bond_obj in self.get_bonds()])
        else:
            msg = ("The informed bond type must be an instance of '%s'."
                   % BondType)
            raise IllegalArgumentError(msg)

    def has_only_bond_type(self, bond_type):
        """Check if this atom has only bonds of type ``bond_type``.

        Parameters
        ----------
        bond_type : `BondType`

        Raises
        ------
        IllegalArgumentError
            If the informed bond type is not an instance of `BondType`.
        """
        if isinstance(bond_type, BondType):
            return all([bond_obj.get_bond_type() == bond_type
                        for bond_obj in self.get_bonds()])
        else:
            msg = ("The informed bond type must be an instance of '%s'."
                   % BondType)
            raise IllegalArgumentError(msg)

    def set_charge(self, charge):
        """Set the formal charge of this atom.

        Parameters
        ----------
        charge : int
        """

        # Both RDKit and Open Babel have the same function name.
        self._atm_obj.SetFormalCharge(charge)

    def set_as_aromatic(self, is_aromatic):
        """Set whether this atom is aromatic or not.

        Parameters
        ----------
        is_aromatic : bool
        """
        if self.is_rdkit_obj():
            self._atm_obj.SetIsAromatic(is_aromatic)
        elif self.is_openbabel_obj():
            self._atm_obj.SetAromatic(is_aromatic)

    def set_in_ring(self, in_ring):
        """Set whether this atom belongs to a ring or not.

        Parameters
        ----------
        in_ring : bool
        """
        if self.is_rdkit_obj():
            # TODO
            error_msg = ("Currently, there is no function in RDKit to define "
                         "if an atom belongs to a ring or not. Please, use "
                         "Open Babel instead.")
            raise NotImplementedError(error_msg)

        elif self.is_openbabel_obj():
            self._atm_obj.SetInRing(in_ring)

    def matches_smarts(self, smarts):
        """Check if this atom matches the substructure through a SMARTS
        substructure search.

        **Note:** currently, this function only works with molecules
        read with Open Babel.

        Parameters
        ----------
        smarts : str
            A substructure defined as SMARTS.

        Returns
        -------
         : bool or None
            Whether matches occurred. Return None if the molecule was read
            with RDKit.

        Examples
        --------
        First, let's read a molecule (glutamine) using Open Babel.

        >>> from luna.wrappers.base import MolWrapper
        >>> mol_obj = MolWrapper.from_smiles("N[C@@H](CCC(N)=O)C(O)=O",
        ...                                  mol_obj_type="openbabel")

        Now, we'll loop over the list of atoms in the glutamine and check which
        atom is the amide's carbon. To do so, we can call the function
        :py:meth:`MolWrapper.get_atoms`, which, by default, returns
        `AtomWrapper` objects and then call :py:meth:`matches_smarts`.

        >>> for atm in mol_obj.get_atoms():
        >>>     print("%d\t%s\t%s" % (atm.get_idx(),
        ...                           atm.get_symbol(),
        ...                           atm.matches_smarts("C(N)(C)=O")))
        1   N   False
        2   C   False
        3   C   False
        4   C   False
        5   C   True
        6   N   False
        7   O   False
        8   C   False
        9   O   False
        10  O   False
        """
        if self.is_rdkit_obj():
            # TODO: Implement
            error_msg = ("Currently, matches_smarts() does not support RDKit "
                         "objects. Please, use Open Babel objects instead.")
            raise NotImplementedError(error_msg)

        elif self.is_openbabel_obj():
            ob_smart = ob.OBSmartsPattern()
            ob_smart.Init(smarts)

            if ob_smart.Match(self.parent):
                for match in ob_smart.GetMapList():
                    if match[0] == self.get_idx():
                        return True
            return False

    def unwrap(self):
        """Return the original atomic object.

        Returns
        -------
         : :class:`rdkit.Chem.rdchem.Atom` or :class:`openbabel.OBAtom`
        """
        return self._atm_obj

    def is_rdkit_obj(self):
        """Check if this atom is an RDKit object."""
        if isinstance(self._atm_obj, RDAtom):
            return True
        return False

    def is_openbabel_obj(self):
        """Check if this atom is an Open Babel object."""
        if isinstance(self._atm_obj, ob.OBAtom):
            return True
        return False

    def __getattr__(self, attr):
        return getattr(self._atm_obj, attr)

    def __sub__(self, other):
        if not isinstance(other, AtomWrapper):
            other = AtomWrapper(other)

        return euclidean_distance(self.get_coords(), other.get_coords())


class BondWrapper:
    """This class provides util functions to access bond properties and
    other information from RDKit and Open Babel objects.

    Parameters
    ----------
    bond_obj : `BondWrapper`, :py:class:`rdkit.Chem.rdchem.Bond`, \
            or :py:class:`openbabel.OBBond`
        A bond to wrap.

    Raises
    ------
    BondObjectTypeError
        If the bond object is not an instance
            of `BondWrapper`, :class:`rdkit.Chem.rdchem.Bond`, \
                or :class:`openbabel.pybel.OBBond`.
    """

    def __init__(self, bond_obj):
        if isinstance(bond_obj, self.__class__):
            bond_obj = bond_obj.unwrap()

        if (not isinstance(bond_obj, RDBond)
                and not isinstance(bond_obj, ob.OBBond)):
            error_msg = ("Objects of type '%s' are not currently accepted."
                         % bond_obj.__class__)
            logger.exception(error_msg)
            raise BondObjectTypeError(error_msg)

        self._bond_obj = bond_obj

    @property
    def bond_obj(self):
        """:class:`rdkit.Chem.rdchem.Bond` or :class:`openbabel.OBBond`: \
                The wrapped bond object."""
        return self._bond_obj

    @bond_obj.setter
    def bond_obj(self, bond_obj):
        if (not isinstance(bond_obj, RDBond)
                and not isinstance(bond_obj, ob.OBBond)):
            error_msg = ("Objects of type '%s' are not currently accepted."
                         % bond_obj.__class__)
            logger.exception(error_msg)
            raise AtomObjectTypeError(error_msg)
        else:
            self._bond_obj = bond_obj

    def get_partner_atom(self, atm, wrapped=True):
        """Get the partner atom that forms this bond with ``atm``.

        Parameters
        ----------
        atm : `AtomWrapper`, :py:class:`rdkit.Chem.rdchem.Atom`, \
                    or :class:`openbabel.OBAtom`
            Get the partner of this atom.
        wrapped : bool
            If True, wrap the partner atom with `AtomWrapper`.

        Returns
        -------
         : `AtomWrapper`, :py:class:`rdkit.Chem.rdchem.Atom`, \
                or :class:`openbabel.OBAtom`
        """
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
        """Return the bond’s first atom.

        Parameters
        ----------
        wrapped : bool
            If True, wrap the atom with `AtomWrapper`.

        Returns
        -------
         : `AtomWrapper`, :py:class:`rdkit.Chem.rdchem.Atom`, \
                or :class:`openbabel.OBAtom`
        """

        # Both RDKit and Open Babel have the same function name.
        if wrapped:
            return AtomWrapper(self._bond_obj.GetBeginAtom())
        return self._bond_obj.GetBeginAtom()

    def get_end_atom(self, wrapped=True):
        """Return the bond’s second atom.

        Parameters
        ----------
        wrapped : bool
            If True, wrap the atom with `AtomWrapper`.

        Returns
        -------
         : `AtomWrapper`, :py:class:`rdkit.Chem.rdchem.Atom`, \
                or :class:`openbabel.OBAtom`
        """

        # Both RDKit and Open Babel have the same function name.
        if wrapped:
            return AtomWrapper(self._bond_obj.GetEndAtom())
        return self._bond_obj.GetEndAtom()

    def get_bond_type(self):
        """Get the bond type (e.g., single bond).

        Returns
        -------
         : `BondType`
        """
        if self.is_rdkit_obj():
            return BondType[self._bond_obj.GetBondType().name]
        elif self.is_openbabel_obj():
            # Open Babel does not have a type specific for aromatic bonds.
            # Instead, they use a flag to specify if a bond is aromatic
            # or not. Therefore, if a bond is aromatic, then it returns
            # BondType.AROMATIC.
            if self._bond_obj.IsAromatic():
                return BondType["AROMATIC"]

            # Map Open Babel bonds to BondType.
            return BondType[OBBondType(self._bond_obj.GetBondOrder()).name]

    def is_aromatic(self):
        """Check if this bond is aromatic or not."""
        if self.is_rdkit_obj():
            return self._bond_obj.GetIsAromatic()
        elif self.is_openbabel_obj():
            return self._bond_obj.IsAromatic()

    def set_bond_type(self, bond_type):
        """Set the type of the bond as a ``bond_type``.

        Parameters
        ----------
        bond_type : `BondType`
        """

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
                raise KeyError("The bond type '%s' is not a valid %s "
                               "bond type." % (bond_type.name, tool))
        else:
            error_msg = ("The informed bond type must be an instance "
                         "of '%s'." % BondType)
            raise IllegalArgumentError(error_msg)

        if self.is_rdkit_obj():
            try:
                self._bond_obj.SetBondType(bond_type)
            except Exception as e:
                logger.exception(e)
                raise
        elif self.is_openbabel_obj():
            self._bond_obj.SetBondOrder(bond_type)

    def set_as_aromatic(self, is_aromatic):
        """Set if this bond is aromatic or not.

        Parameters
        ----------
        is_aromatic : bool
        """
        if self.is_rdkit_obj():
            self._bond_obj.SetIsAromatic(is_aromatic)
        elif self.is_openbabel_obj():
            self._bond_obj.SetAromatic(is_aromatic)

    def unwrap(self):
        """Return the original bond object.

        Returns
        -------
         : :class:`rdkit.Chem.rdchem.Bond` or :class:`openbabel.OBBond`
        """
        return self._bond_obj

    def is_rdkit_obj(self):
        """Check if this bond is an RDKit object."""
        if isinstance(self._bond_obj, RDBond):
            return True
        return False

    def is_openbabel_obj(self):
        """Check if this bond is an Open Babel object."""
        if isinstance(self._bond_obj, ob.OBBond):
            return True
        return False

    def __getattr__(self, attr):
        return getattr(self._bond_obj, attr)


class MolWrapper:
    """This class provides util functions to access molecule properties and
    other information from RDKit and Open Babel objects.

    Parameters
    ----------
    mol_obj : `MolWrapper`, :class:`rdkit.Chem.rdchem.Mol`, \
                    or :class:`openbabel.pybel.Molecule`
        A molecule to wrap.

    Raises
    ------
    MoleculeObjectTypeError
        If the molecular object is not an instance
            of `MolWrapper`, :class:`rdkit.Chem.rdchem.Mol`, \
                or :class:`openbabel.pybel.Molecule`.
    """

    def __init__(self, mol_obj):
        if isinstance(mol_obj, self.__class__):
            mol_obj = mol_obj.unwrap()
        elif isinstance(mol_obj, PybelMol):
            mol_obj = mol_obj.OBMol

        if (not isinstance(mol_obj, RDMol)
                and not isinstance(mol_obj, ob.OBMol)):
            msg = ("Objects of type '%s' are not currently accepted."
                   % mol_obj.__class__)
            logger.exception(msg)
            raise MoleculeObjectTypeError(msg)

        self._mol_obj = mol_obj

    @classmethod
    def from_mol_file(cls, file, mol_format, mol_obj_type="rdkit"):
        """Initialize a molecule from a molecular file.

        Parameters
        ----------
        file : str
            The molecular file.
        mol_format : str
            Define the format in which the molecule is represented
                (e.g., 'mol2' or 'mol').
        mol_obj_type : {'rdkit', 'openbabel'}
            Define which library (RDKit or Open Babel) to use to parse the
            molecular block. The default value is RDKit.

        Returns
        -------
         : `MolWrapper`
        """
        if mol_obj_type == "rdkit":
            return cls(read_mol_from_file(file, mol_format))
        elif mol_obj_type == "openbabel":
            return cls(list(readfile(mol_format, file))[0])

    @classmethod
    def from_mol_block(cls, block, mol_format, mol_obj_type="rdkit"):
        """Initialize a molecule from a string block.

        Parameters
        ----------
        block : str
            The molecular string block.
        mol_format : str
            Define the format in which the molecule is represented
            (e.g., 'mol2' or 'mol').
        mol_obj_type : {'rdkit', 'openbabel'}
            Define which library (RDKit or Open Babel) to use to parse the
            molecular block. The default value is RDKit.

        Returns
        -------
         : `MolWrapper`
        """
        if mol_obj_type == "rdkit":
            return cls(new_mol_from_block(block, mol_format))
        elif mol_obj_type == "openbabel":
            return cls(readstring(mol_format, block))

    @classmethod
    def from_smiles(cls, smiles, mol_obj_type="rdkit", name=None):
        """Initialize a molecule from a SMILES string.

        Parameters
        ----------
        smiles : str
            The SMILES string.
            Define the format in which the molecule is represented
            (e.g., 'mol2' or 'mol').
        mol_obj_type : {'rdkit', 'openbabel'}
            Define which library (RDKit or Open Babel) to use to parse the
            molecular block. The default value is RDKit.
        name : str, optional
            A name to identify the molecule.

        Returns
        -------
         : `MolWrapper`

        Raises
        ------
        MoleculeObjectError
            If it could not create a molecule from the provided SMILES.

        Examples
        --------

        >>> from luna.wrappers.base import MolWrapper
        >>> mol_obj = MolWrapper.from_smiles("N[C@@H](CCC(N)=O)C(O)=O",
        ...                                  mol_obj_type="openbabel",
        ...                                  name="Glutamine")

        """
        if mol_obj_type == "rdkit":
            mol = MolFromSmiles(smiles)
        elif mol_obj_type == "openbabel":
            mol = readstring("smi", smiles)

        # Raise an error if the molecule could not be created.
        if mol is None:
            msg = ("It could not create a molecule from the "
                   "provided SMILES '%s'." % smiles)
            raise MoleculeObjectError(msg)

        mol = cls(mol)
        mol.set_name(name or smiles)

        return mol

    @property
    def mol_obj(self):
        """:class:`rdkit.Chem.rdchem.Mol` or \
                :class:`openbabel.pybel.Molecule`: \
            The wrapped molecular object."""
        return self._mol_obj

    @mol_obj.setter
    def mol_obj(self, mol_obj):
        if isinstance(mol_obj, PybelMol):
            mol_obj = mol_obj.OBMol

        if (not isinstance(mol_obj, RDMol)
                and not isinstance(mol_obj, ob.OBMol)):
            msg = ("Objects of type '%s' are not currently accepted."
                   % mol_obj.__class__)
            logger.exception(msg)
            raise MoleculeObjectTypeError(msg)

        self._mol_obj = mol_obj

    def as_rdkit(self):
        """If the molecule is an Open Babel object, convert it to an
        RDKit object."""
        if self.is_rdkit_obj():
            return self._mol_obj
        elif self.is_openbabel_obj():
            new_mol = MolWrapper.from_smiles(self.to_smiles(),
                                             mol_obj_type="rdkit")
            new_mol.set_name(self.get_name())
            return new_mol.unwrap()

    def as_openbabel(self):
        """If the molecule is an RDKit object, convert it to an
        Open Babel object."""
        if self.is_openbabel_obj():
            return self._mol_obj
        elif self.is_rdkit_obj():
            new_mol = MolWrapper.from_smiles(self.to_smiles(),
                                             mol_obj_type="openbabel")
            new_mol.set_name(self.get_name())
            return new_mol.unwrap()

    def get_name(self):
        """Get the molecule name.

        Returns
        -------
         : str
        """
        if self.is_rdkit_obj():
            if self._mol_obj.HasProp("_Name"):
                return self._mol_obj.GetProp("_Name")
            return ""
        elif self.is_openbabel_obj():
            return self._mol_obj.GetTitle()

    def get_atoms(self, wrapped=True):
        """Get all molecule's atoms.

        Parameters
        ----------
        wrapped : bool
            If True, wrap all atoms with `AtomWrapper`.

        Returns
        -------
         : iterable of `AtomWrapper`, :py:class:`rdkit.Chem.rdchem.Atom`, \
                or :class:`openbabel.OBAtom`
        """
        atoms = []
        if self.is_rdkit_obj():
            atoms = self._mol_obj.GetAtoms()
        elif self.is_openbabel_obj():
            atoms = ob.OBMolAtomIter(self._mol_obj)

        if atoms and wrapped:
            return [AtomWrapper(atm, self) for atm in atoms]
        return atoms

    def get_atom_by_idx(self, idx, wrapped=True):
        """Return the atom whose index is ``idx``."""
        try:
            atm = [a for a in self.get_atoms() if a.get_idx() == idx][0]

            if wrapped:
                return atm
            return atm.unwrap()
        except Exception:
            return None

    def get_bonds(self, wrapped=True):
        """Get all molecule's bonds.

        Parameters
        ----------
        wrapped : bool
            If True, wrap all bonds with `BondWrapper`.

        Returns
        -------
         : iterable of `BondWrapper`, :class:`rdkit.Chem.rdchem.Bond`, \
                or :class:`openbabel.OBBond`
        """
        bonds = []
        if self.is_rdkit_obj():
            bonds = self._mol_obj.GetBonds()
        elif self.is_openbabel_obj():
            bonds = ob.OBMolBondIter(self._mol_obj)

        if bonds and wrapped:
            return [BondWrapper(bond) for bond in bonds]
        return bonds

    def get_total_charge(self):
        """Get total molecular charge."""
        if self.is_rdkit_obj():
            return GetFormalCharge(self._mol_obj)
        elif self.is_openbabel_obj():
            return self._mol_obj.GetTotalCharge()

    def get_num_heavy_atoms(self):
        """Get the number of heavy atoms in this molecule."""
        if self.is_rdkit_obj():
            return self._mol_obj.GetNumHeavyAtoms()
        elif self.is_openbabel_obj():
            return self._mol_obj.NumHvyAtoms()

    def get_atom_coord_by_id(self, atm_id):
        """Get the coordinates of an atom given by its id.

        Returns
        -------
         : array_like of float (size 3)
            Atomic coordinates (x, y, z).
        """
        if self.is_rdkit_obj():
            return list(self._mol_obj.GetConformer().GetAtomPosition(atm_id))
        elif self.is_openbabel_obj():
            atm = self._mol_obj.GetAtomById(atm_id)
            return [atm.GetX(), atm.GetY(), atm.GetZ()]

    def get_obj_type(self):
        """Get the object type ("rdkit" or "openbabel")."""
        if self.is_rdkit_obj():
            return "rdkit"
        elif self.is_openbabel_obj():
            return "openbabel"

    def find_mcs(self, other, match_valences=True,
                 ring_matches_ring_only=True, complete_rings_only=True):
        other = MolWrapper(other)
        other = other.unwrap() if other.is_rdkit_obj() else other.as_rdkit()

        mol_obj = self.unwrap() if self.is_rdkit_obj() else self.as_rdkit()

        # Maximum Common Substructure.
        return rdFMCS.FindMCS([other, mol_obj],
                              matchValences=match_valences,
                              ringMatchesRingOnly=ring_matches_ring_only,
                              completeRingsOnly=complete_rings_only)

    def has_name(self):
        """Check if this molecule has a name."""
        if self.get_name():
            return True
        return False

    def set_name(self, name):
        """Set a name to this molecule.

        Parameters
        ----------
        name : str
        """
        if self.is_rdkit_obj():
            self._mol_obj.SetProp("_Name", name)
        elif self.is_openbabel_obj():
            self._mol_obj.SetTitle(name)

    def to_smiles(self):
        """Return the canonical SMILES string for this molecule."""
        if self.is_rdkit_obj():
            return MolToSmiles(self._mol_obj)
        elif self.is_openbabel_obj():
            return PybelMol(self._mol_obj).write("smi").split("\t")[0]

    def to_pdb_block(self):
        """Return the PDB string block for this molecule."""
        if self.is_rdkit_obj():
            return MolToPDBBlock(self._mol_obj)
        elif self.is_openbabel_obj():
            return PybelMol(self._mol_obj).write("pdb")

    def to_mol_block(self):
        """Return the MOL string block for this molecule."""
        if self.is_rdkit_obj():
            return MolToMolBlock(self._mol_obj)
        elif self.is_openbabel_obj():
            return PybelMol(self._mol_obj).write("mol")

    def unwrap(self):
        """Return the original molecular object.

        Returns
        -------
         : :class:`rdkit.Chem.rdchem.Mol` or :class:`openbabel.pybel.Molecule`.
        """
        return self._mol_obj

    def is_rdkit_obj(self):
        """Check if this molecule is an RDKit object."""
        if isinstance(self._mol_obj, RDMol):
            return True
        return False

    def is_openbabel_obj(self):
        """Check if this molecule is an Open Babel object."""
        if isinstance(self._mol_obj, ob.OBMol):
            return True
        elif isinstance(self._mol_obj, PybelMol):
            return True
        return False

    def is_pybel_obj(self):
        """Check if this molecule is a Pybel object."""
        if isinstance(self._mol_obj, PybelMol):
            return True
        return False

    def __getattr__(self, attr):
        if self.mol_obj is None:
            return None
        return getattr(self._mol_obj, attr)

    def __getstate__(self):
        # Creates a copy of the class' dictionary in case we need to modify the
        # molecular object to pickle it.
        my_dict = self.__dict__.copy()
        if self.is_openbabel_obj():
            my_dict["_mol_block"] = self.to_mol_block()
            my_dict["_mol_obj_type"] = self.get_obj_type()
            my_dict["_mol_obj"] = None
        return my_dict

    def __setstate__(self, state):
        if "_mol_obj_type" in state and "_mol_block" in state:
            state["_mol_obj"] = \
                self.from_mol_block(state["_mol_block"],
                                    "mol", "openbabel").unwrap()
            del state["_mol_obj_type"]
            del state["_mol_block"]

        self.__dict__.update(state)
