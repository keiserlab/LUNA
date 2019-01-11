class NbCoordinates():

    def __init__(self, coords):
        self.coords = coords

    @property
    def size(self):
        return len(self.coords)

    def __eq__(self, other):
        """Overrides the default implementation"""
        if isinstance(self, other.__class__):
            return self.coords == other.coords
        return False

    def __ne__(self, other):
        """Overrides the default implementation"""
        return not self.__eq__(other)

    def __hash__(self):
        """Overrides the default implementation"""
        # The list is sorted in order to avoid dependence on appending order.
        coords_tuple = tuple(sorted(self.coords, key=hash))
        return hash(coords_tuple)


class NbAtom():

    def __init__(self, mybio_atom, nb_coods, atm_grps=None):
        self._atom = mybio_atom
        self._nb_coords = nb_coods
        self._atm_grps = atm_grps or []

    @property
    def atom(self):
        return self._atom

    @property
    def nb_coords(self):
        return self._nb_coords.coords

    @property
    def num_neighbours(self):
        return self._nb_coords.size

    @property
    def atm_grps(self):
        return self._atm_grps

    def add_atm_grp(self, group):
        self._atm_grps = list(set(self.atm_grps + [group]))

    def __getattr__(self, attr):
        return getattr(self._atom, attr)

    def __repr__(self):
        return "<NBAtom %s>" % self._atom.__repr__()

    def __eq__(self, other):
        """Overrides the default implementation"""
        if isinstance(self, other.__class__):
            return self._atom == other._atom and self._nb_coords == other._nb_coords
        return False

    def __ne__(self, other):
        """Overrides the default implementation"""
        return not self.__eq__(other)

    def __hash__(self):
        """Overrides the default implementation"""
        return hash((self._atom, self._nb_coords))
