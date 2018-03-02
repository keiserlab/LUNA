class VicinityAtom():

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
        return "<VicinityAtom %s>" % self._atom.__repr__()
