import interaction.util as iu
import interaction.plane_regression as plane


class RingPlane:

    coordinates = []
    centroid = None
    normal = None

    def __init__(self, atoms):
        self.coordinates = iu.get_coordinatesnp(atoms)
        self.centroid = iu.calc_centroidnp(self.coordinates)

    def calc_normal(self):
        if (self.normal is None):
            self.normal = plane.calc_normal(self.coordinates)
