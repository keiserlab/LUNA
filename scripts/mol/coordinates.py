import numpy as np


class Coordinate():

    def __init__(self, x, y, z, atomicNumber=None, atomId=None):
        self._coord = (x, y, z)
        self.atomicNumber = atomicNumber
        self.atomId = atomId

    def __repr__(self):
        return ("<Coordinate x=%.3f, y=%.3f, z=%.3f>"
                % (self.x, self.y, self.z))

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


class NeighbourhoodCoordinates():

    def __init__(self, coords):
        self.coords = coords

    @property
    def size(self):
        return len(self.coords)
