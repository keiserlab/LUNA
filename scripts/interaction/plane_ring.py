import interaction.util as iu
import interaction.plane_regression as plane


class RingPlane:

    def __init__(self, atoms):
        self.coordinates = iu.get_coordinatesnp(atoms)
        self.centroid = iu.calc_centroidnp(self.coordinates)
        self.normal = None

    def calc_normal(self):
        if (self.normal is None):
            self.normal = plane.calc_normal(self.coordinates)
