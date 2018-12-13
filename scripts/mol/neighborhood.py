import numpy as np


class Coordinate():

    def __init__(self, x, y, z, atomic_num=None, atom_num=None):
        self._coord = (x, y, z)
        self.atomic_num = atomic_num
        self.atom_num = atom_num

    def __repr__(self):
        return ("<Coordinate x=%.3f, y=%.3f, z=%.3f>" % (self.x, self.y, self.z))

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
    def coord(self):
        return self._coord

    @coord.setter
    def coord(self, xyz):
        self._coord = np.array(xyz)


class NbCoordinates():

    def __init__(self, coords):
        self.coords = coords

    @property
    def size(self):
        return len(self.coords)


class NbAtom():

    def __init__(self, myBioPDBAtom, neighbourCoords):
        self._atom = myBioPDBAtom
        self._nbCoords = neighbourCoords

    @property
    def nbCoords(self):
        return self._nbCoords

    @property
    def numNeighbours(self):
        return self._nbCoords.size

    def __getattr__(self, attr):
        return getattr(self._atom, attr)

    def __repr__(self):
        return "<NBAtom %s>" % self._atom.__repr__()
