import numpy as np


class Coordinate():

    def __init__(self, x, y, z, atomic_num=None, atom_id=None):
        self._coord = np.array([x, y, z])
        self.atomic_num = atomic_num
        self.atom_id = atom_id

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
        return self.coord

    @property
    def coord(self):
        return self._coord

    @coord.setter
    def coord(self, xyz):
        self._coord = np.array(xyz)

    def __repr__(self):
        return ("<Coordinate x=%.3f, y=%.3f, z=%.3f>" % (self.x, self.y, self.z))

    def __eq__(self, other):
        """Overrides the default implementation"""
        if isinstance(self, other.__class__):
            return np.all(self._coord == other._coord)
        return False

    def __ne__(self, other):
        """Overrides the default implementation"""
        return not self.__eq__(other)

    def __hash__(self):
        """Overrides the default implementation"""
        return hash(tuple(self._coord))
