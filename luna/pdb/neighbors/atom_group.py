import numpy as np
from Bio.PDB.kdtrees import KDTree


class AtomGroupNeighborhood:
    """ Class for fast neighbor atom groups searching.

    ``AtomGroupNeighborhood`` makes use of a KD Tree implemented in C,
    so it's fast.

    Parameters
    ----------
    atm_grps : iterable of `AtomGroup`, optional
        A sequence of `AtomGroup` objects, which is used in the queries.
        It can contain atom groups from different molecules.
    bucket_size : int
        Bucket size of KD tree.
        You can play around with this to optimize speed if you feel like it.
        The default value is 10.

    """

    def __init__(self, atm_grps, bucket_size=10):
        self.atm_grps = list(atm_grps)

        # get the coordinates
        coord_list = [ga.centroid for ga in self.atm_grps]

        # to Nx3 array of type float
        self.coords = np.array(coord_list).astype("d")
        assert(bucket_size > 1)
        assert(self.coords.shape[1] == 3)
        self.kdt = KDTree(self.coords, bucket_size)    

    def search(self, center, radius):
        """Return all atom groups in ``atm_grps`` that is up to a maximum of
        ``radius`` away (measured in Å) of ``center``.

        For atom groups with more than one atom, their centroid is used as a
        reference.
        """
        center = np.require(center, dtype="d", requirements="C")
        if center.shape != (3,):
            raise Exception("Expected a 3-dimensional NumPy array")
        points = self.kdt.search(center, radius)            
        n_grps_list = [self.atm_grps[point.index] for point in points]

        return n_grps_list
