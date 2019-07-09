import numpy as np
from openbabel import etab

from graph.bellman_ford import bellman_ford


class AtomData:

    def __init__(self, atomic_num, coord, serial_number=None):
        self.atomic_num = atomic_num
        self._coord = np.array(coord)
        self.serial_number = serial_number

    @property
    def x(self):
        return self._coord[0]

    @property
    def y(self):
        return self._coord[1]

    @property
    def z(self):
        return self._coord[2]

    @property
    def vector(self):
        return self._coord

    @property
    def coord(self):
        return self._coord

    @coord.setter
    def coord(self, xyz):
        self._coord = np.array(xyz)

    def __repr__(self):
        return ("<ExtendedAtomData: atomic number=%d, coord=(%.3f, %.3f, %.3f), serial number=%s>"
                % (self.atomic_num, self.x, self.y, self.z, str(self.serial_number)))

    def __eq__(self, other):
        """Overrides the default implementation"""
        if isinstance(self, other.__class__):
            return (self.atomic_num == other.atomic_num and
                    np.all(self._coord == other._coord) and
                    self.serial_number == other.serial_number)
        return False

    def __ne__(self, other):
        """Overrides the default implementation"""
        return not self.__eq__(other)

    def __hash__(self):
        """Overrides the default implementation"""
        return hash((self.atomic_num, tuple(self._coord), self.serial_number))


class ExtendedAtom:

    def __init__(self, mybio_atom, nb_info=None, atm_grps=None, neighborhood=None):
        self._atom = mybio_atom
        self._nb_info = nb_info or []
        self._atm_grps = atm_grps or []
        self._neighborhood = neighborhood or {}

    @property
    def atom(self):
        return self._atom

    @property
    def neighbors_info(self):
        return self._nb_info

    @property
    def neighborhood(self):
        return self._neighborhood

    @property
    def atm_grps(self):
        return self._atm_grps

    @property
    def electronegativity(self):
        return etab.GetElectroNeg(etab.GetAtomicNum(self.element))

    @property
    def full_atom_name(self):
        full_atom_name = "%s/%s/%s" % self.get_full_id()[0:3]
        res_name = "%s-%d%s" % (self._atom.parent.resname, self._atom.parent.id[1], self._atom.parent.id[2].strip())
        atom_name = "%s" % self._atom.name
        if self.altloc != " ":
            atom_name += "-%s" % self.altloc
        full_atom_name += "/%s/%s" % (res_name, atom_name)

        return full_atom_name

    def set_neighborhood(self, nb_graph):
        self._neighborhood = nb_graph

    def add_nb_info(self, nb_info):
        self._nb_info = list(set(self._nb_info + list(nb_info)))

    def add_atm_grps(self, atm_grps):
        self._atm_grps = list(set(self._atm_grps + list(atm_grps)))

    def remove_nb_info(self, nb_info):
        self._nb_info = list(set(self._nb_info) - set(nb_info))

    def remove_atm_grps(self, atm_grps):
        self._atm_grps = list(set(self._atm_grps) - set(atm_grps))

    def is_neighbor(self, atom):
        return atom.serial_number in [i.serial_number for i in self._nb_info]

    def get_shortest_path_size(self, trgt_atm):
        # d stores the path size from the source to the target, and p stores the predecessors from each target.
        d, p = bellman_ford(self.neighborhood, self.serial_number)
        # Return infinite if the provided target is not in the distance object.
        return d.get(trgt_atm.serial_number, float("inf"))

    def __getattr__(self, attr):
        return getattr(self._atom, attr)

    def __repr__(self):
        return "<ExtendedAtom: %s>" % self.full_atom_name

    def __sub__(self, other):
        # It calls __sub__() from Biopython.
        return self._atom - other._atom

    def __eq__(self, other):
        """Overrides the default implementation"""
        if isinstance(self, other.__class__):
            return self._atom == other._atom
        return False

    def __ne__(self, other):
        """Overrides the default implementation"""
        return not self.__eq__(other)

    def __lt__(self, a2):
        # It calls __lt__() from Biopython
        return self._atom < a2._atom

    def __hash__(self):
        """Overrides the default implementation"""
        return hash(self._atom)
