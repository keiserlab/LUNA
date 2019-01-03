class NbCoordinates():

    def __init__(self, coords):
        self.coords = coords

    @property
    def size(self):
        return len(self.coords)


class NbAtom():

    def __init__(self, myBioPDBAtom, nb_coods, atm_grps=None):
        self._atom = myBioPDBAtom
        self._nb_coords = nb_coods
        self._atm_grps = atm_grps or []

    @property
    def nb_coords(self):
        return self._nb_coords.coords

    @property
    def num_neighbours(self):
        return self._nb_coords.size

    @property
    def atm_grps(self):
        return self._atm_grps

    def add_atom_group(self, group):
        self._atm_grps = list(set(self.atm_grps + [group]))

    def __getattr__(self, attr):
        return getattr(self._atom, attr)

    def __repr__(self):
        return "<NBAtom %s>" % self._atom.__repr__()
