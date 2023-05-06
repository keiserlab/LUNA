import numpy as np
from openbabel import openbabel as ob

import logging
logger = logging.getLogger()


class AtomData:
    """Store atomic data (atomic number, coordinates, bond type,
    and serial number).

    Parameters
    ----------
    atomic_num : int
        Atomic number.
    coord : array_like of float (size 3)
        Atomic coordinates (x, y, z).
    bond_type : int
        Bond type.
    full_id : tuple
        The atom' full id as in :class:`~luna.MyBio.PDB.Atom.Atom`.
    serial_number : int, optional
        Atom serial number.

    Attributes
    ----------
    atomic_num : int
        The atomic number.
    bond_type : int
        The bond type.
    full_id : tuple
        The atom' full id
    serial_number : int or None
        The atom serial number.
    """

    def __init__(self,
                 atomic_num,
                 coord,
                 bond_type,
                 full_id=None,
                 serial_number=None):
        self.atomic_num = atomic_num
        # Standardize all coordinate data to the same Numpy data type
        # for consistence.
        self._coord = np.array(coord, "f")
        self.bond_type = bond_type
        self.full_id = full_id
        self.serial_number = serial_number

    @property
    def x(self):
        """float: The orthogonal coordinates for ``x`` in Angstroms."""
        return self._coord[0]

    @property
    def y(self):
        """float: The orthogonal coordinates for ``y`` in Angstroms."""
        return self._coord[1]

    @property
    def z(self):
        """float: The orthogonal coordinates for ``z`` in Angstroms."""
        return self._coord[2]

    @property
    def coord(self):
        """array_like of float (size 3): The atomic coordinates (x, y, z)."""
        return self._coord

    @coord.setter
    def coord(self, xyz):
        # Standardize all coordinate data to the same Numpy data type for
        # consistence.
        self._coord = np.array(xyz, "f")

    def __repr__(self):
        full_atom_name = ""
        if self.full_id is not None:
            full_atom_name = "%s/%s/%s" % self.full_id[0:3]
            res_name = "%d%s" % (self.full_id[3][1],
                                 self.full_id[3][2].strip())
            atom_name = "%s" % self.full_id[4][0]

            if self.full_id[4][1] != " ":
                atom_name += "-%s" % self.full_id[4][1]
            full_atom_name += "/%s/%s" % (res_name, atom_name)

        return ("<ExtendedAtomData: atomic number=%d, "
                "coord=(%.3f, %.3f, %.3f), atom='%s', serial number=%s>"
                % (self.atomic_num, self.x, self.y, self.z, full_atom_name,
                   str(self.serial_number)))

    def __eq__(self, other):
        """Overrides the default implementation"""
        if isinstance(self, other.__class__):
            return (self.atomic_num == other.atomic_num
                    and np.all(self._coord == other._coord)
                    and self.full_id == other.full_id
                    and self.serial_number == other.serial_number)
        return False

    def __ne__(self, other):
        """Overrides the default implementation"""
        return not self.__eq__(other)

    def __hash__(self):
        """Overrides the default implementation"""
        return hash((self.atomic_num, tuple(self._coord),
                     self.full_id, self.bond_type,
                     self.serial_number))


class ExtendedAtom:
    """Extend :class:`~luna.MyBio.PDB.Atom.Atom` with additional properties
    and methods.

    Parameters
    ----------
    atom : :class:`~luna.MyBio.PDB.Atom.Atom`
        An atom.
    nb_info : iterable of `AtomData`, optional
        A sequence of `AtomData` containing information about atoms covalently
        bound to ``atom``.
    atm_grps : iterable of :class:`~luna.groups.AtomGroup`, optional
        A sequence of atom groups that contain ``atom``.
    invariants : list or tuple, optional
        Atomic invariants.
    """

    def __init__(self, atom, nb_info=None, atm_grps=None, invariants=None):
        self._atom = atom
        self._nb_info = nb_info or []
        self._atm_grps = atm_grps or []
        self._invariants = invariants

    @property
    def atom(self):
        """:class:`~luna.MyBio.PDB.Atom.Atom`, read-only."""
        return self._atom

    @property
    def neighbors_info(self):
        """list of `AtomData`, read-only: The list of `AtomData`
        containing information about atoms covalently bound to ``atom``.

        To add or remove neighbors information from ``neighbors_info``
        use :py:meth:`add_nb_info` or :py:meth:`remove_nb_info`,
        respectively."""
        return self._nb_info

    @neighbors_info.setter
    def neighbors_info(self, nb_info):
        self._nb_info = nb_info

    @property
    def atm_grps(self):
        """list of :class:`~luna.groups.AtomGroup`, read-only: The list of
        atom groups that contain ``atom``.

        To add or remove atom groups from ``atm_grps`` use
        :py:meth:`add_atm_grps` or :py:meth:`remove_atm_grps`, respectively."""
        return self._atm_grps

    @property
    def invariants(self):
        """list: The list of atomic invariants."""
        return self._invariants

    @invariants.setter
    def invariants(self, invariants):
        self._invariants = invariants

    @property
    def atomic_num(self):
        """int, read-only: This atom's atomic number. This information is
        obtained from Open Babel."""
        return ob.GetAtomicNum(self.element)

    @property
    def electronegativity(self):
        """float, read-only: The Pauling electronegativity for this atom.
        This information is obtained from Open Babel."""
        return ob.GetElectroNeg(ob.GetAtomicNum(self.element))

    @property
    def full_id(self):
        """tuple, read-only: The full id of an atom is the tuple (structure id,
        model id, chain id, residue id, atom name, alternate location)."""
        return self._atom.get_full_id()

    @property
    def full_atom_name(self):
        """str, read-only: The full name of an atom is composed by the
                structure id, model id, chain id, residue name, residue id,
                atom name, and alternate location if available.
                Fields are slash-separated."""
        full_atom_name = "%s/%s/%s" % self.get_full_id()[0:3]
        res_name = "%s/%d%s" % (self._atom.parent.resname,
                                self._atom.parent.id[1],
                                self._atom.parent.id[2].strip())
        atom_name = "%s" % self._atom.name
        if self.altloc != " ":
            atom_name += "-%s" % self.altloc
        full_atom_name += "/%s/%s" % (res_name, atom_name)

        return full_atom_name

    def add_nb_info(self, nb_info):
        """ Add `AtomData` objects to ``neighbors_info``."""
        self._nb_info = list(set(self._nb_info + list(nb_info)))

    def add_atm_grps(self, atm_grps):
        """ Add :class:`~luna.groups.AtomGroup` objects to ``atm_grps``."""
        self._atm_grps = list(set(self._atm_grps + list(atm_grps)))

    def remove_nb_info(self, nb_info):
        """ Remove `AtomData` objects from ``neighbors_info``."""
        self._nb_info = list(set(self._nb_info) - set(nb_info))

    def remove_atm_grps(self, atm_grps):
        """ Remove :class:`~luna.groups.AtomGroup` objects from
        ``atm_grps``."""
        self._atm_grps = list(set(self._atm_grps) - set(atm_grps))

    def get_neighbor_info(self, atom):
        """Get information from a covalently bound atom."""
        for info in self._nb_info:
            if atom.get_full_id() == info.full_id:
                return info

    def is_neighbor(self, atom):
        """Check if a given atom is covalently bound to it."""
        return atom.get_full_id() in [i.full_id for i in self._nb_info]

    def as_json(self):
        """Represent the atom as a dict containing the structure id, model id,
           chain id, residue name, residue id, and atom name.

           The dict is defined as follows:

            * ``pdb_id`` (str): structure id;
            * ``model`` (str): model id;
            * ``chain`` (str): chain id;
            * ``res_name`` (str): residue name;
            * ``res_id`` (tuple): residue id (hetflag, sequence identifier, \
                                              insertion code);
            * ``name`` (tuple): atom name (atom name, alternate location).

        """
        full_id = self.get_full_id()

        return {"pdb_id": full_id[0],
                "model": full_id[1],
                "chain": full_id[2],
                "res_name": self.parent.resname,
                "res_id": full_id[3],
                "name": full_id[4]}

    def __getattr__(self, attr):
        if hasattr(self._atom, attr):
            return getattr(self._atom, attr)
        else:
            raise AttributeError("The attribute '%s' does not exist in the "
                                 "class %s." % (attr, self.__class__.__name__))

    def __getstate__(self):
        return self.__dict__

    def __setstate__(self, state):
        self.__dict__.update(state)

    def __repr__(self):
        return "<ExtendedAtom: %s>" % self.full_atom_name

    def __sub__(self, other):
        # It calls __sub__() from Biopython.
        return self._atom - other._atom

    def __eq__(self, other):
        """Overrides the default implementation"""
        if isinstance(self, other.__class__):
            return self.full_atom_name == other.full_atom_name
        return False

    def __ne__(self, other):
        """Overrides the default implementation"""
        return not self.__eq__(other)

    def __lt__(self, a2):
        # It substitutes the residue id for its index in order to keep
        # the same order as in the PDB.
        full_id1 = self.full_id[0:2] + (self.parent.idx, ) + self.full_id[4:]
        full_id2 = a2.full_id[0:2] + (a2.parent.idx, ) + a2.full_id[4:]

        return full_id1 < full_id2

    def __hash__(self):
        """Overrides the default implementation"""
        return hash(self.full_atom_name)
