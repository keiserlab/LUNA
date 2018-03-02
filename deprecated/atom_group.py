import interaction.util as iu


class AtomGroup():

    def __init__(self, atoms, chemicalFeatures):
        self.atoms = atoms
        self.chemicalFeatures = chemicalFeatures
        self.coords = iu.get_coordinatesnp(atoms)

        self._centroid = iu.calc_centroidnp(self.coords)
        self._normal = None
        self._isNormalCalculated = False

    def get_serial_numbers(self):
        return [a.get_serial_number() for a in self.atoms]

    @property
    def centroid(self):
        return self._centroid

    @property
    def normal(self):
        if self._isNormalCalculated is False:
            self._calc_normal()

        return self._normal

    def _calc_normal(self):
        import interaction.plane_regression as plane

        if self._isNormalCalculated is False:
            self._normal = plane.calc_normal(self.coords)
            self._isNormalCalculated = True
