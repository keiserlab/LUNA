import re
import logging

from rdkit import DataStructs
from rdkit.ML.Cluster import Butina

from luna.util.exceptions import IllegalArgumentError

logger = logging.getLogger()


def available_similarity_functions():
    """Return a list of all similarity metrics available at RDKit."""
    regex = re.compile("Bulk([a-zA-Z]+)Similarity", flags=0)
    return list(filter(regex.match, dir(DataStructs)))


def calc_distance_matrix(fps, similarity_func="BulkTanimotoSimilarity"):
    """Calculate the pairwise distance (dissimilarity) between fingerprints
     in ``fps`` using the similarity metric ``similarity_func``.

    Parameters
    ----------
    fps : iterable of RDKit \
            :class:`~rdkit.DataStructs.cDataStructs.ExplicitBitVect` or \
            :class:`~rdkit.DataStructs.cDataStructs.SparseBitVect`
        A sequence of fingerprints.
    similarity_func : str
        A similarity metric to calculate the distance between the provided
        fingerprints. The default value is 'BulkTanimotoSimilarity'.

        To check out the list of available similarity metrics, call the
        function :py:meth:`available_similarity_functions`.

    Examples
    --------

    First, let's define a set of molecules.

    >>> from luna.wrappers.base import MolWrapper
    >>> mols = [MolWrapper.from_smiles("CCCCCC").unwrap(),
    ...         MolWrapper.from_smiles("CCCCCCCC").unwrap(),
    ...         MolWrapper.from_smiles("CCCCCCCCO").unwrap()]

    Now, we generate fingerprints for those molecules.

    >>> from luna.mol.fingerprint import generate_fp_for_mols
    >>> fps = [d["fp"] for d in generate_fp_for_mols(mols, "morgan_fp")]

    Finally, calculate the distance between the molecules based on their
    fingerprints.

    >>> from luna.mol.clustering import calc_distance_matrix
    >>> print(calc_distance_matrix(fps))
    [0.125, 0.46153846153846156, 0.3846153846153846]

    Returns
    -------
    distances : list of float
        Flattened diagonal matrix.
    """
    funcs = available_similarity_functions()
    if similarity_func not in funcs:
        raise IllegalArgumentError("Similarity function not available.")

    dists = []
    for i in range(1, len(fps)):
        if (similarity_func == "BulkTverskySimilarity"):
            params = [fps[i], fps[:i], 0, 1]
        else:
            params = [fps[i], fps[:i]]

        sims = getattr(DataStructs, similarity_func)(*params)
        dists.extend([1 - x for x in sims])

    return dists


def cluster_fps(fps, cutoff=0.2, similarity_func="BulkTanimotoSimilarity"):
    """Clusterize molecules based on fingerprints using the Butina
    clustering algorithm.

    Parameters
    ----------
    fps : iterable of RDKit \
            :class:`~rdkit.DataStructs.cDataStructs.ExplicitBitVect` or
            :class:`~rdkit.DataStructs.cDataStructs.SparseBitVect`
        A sequence of fingerprints.
    cutoff : float
        Elements within this range of each other are considered
        to be neighbors.
    similarity_func : str
        A similarity metric to calculate the distance between the provided
        fingerprints.  The default value is 'BulkTanimotoSimilarity'.

        To check out the list of available similarity metrics, call the
        function :py:meth:`available_similarity_functions`.

    Examples
    --------

    First, let's define a set of molecules.

    >>> from luna.wrappers.base import MolWrapper
    >>> mols = [MolWrapper.from_smiles("CCCCCC").unwrap(),
    ...         MolWrapper.from_smiles("CCCCCCCC").unwrap(),
    ...         MolWrapper.from_smiles("CCCCCCCCO").unwrap()]

    Now, we generate fingerprints for those molecules.

    >>> from luna.mol.fingerprint import generate_fp_for_mols
    >>> fps = [d["fp"] for d in generate_fp_for_mols(mols, "morgan_fp")]

    Finally, clusterize the molecules based on their fingerprints.

    >>> from luna.mol.clustering import cluster_fps
    >>> print(cluster_fps(fps, cutoff=0.2))
    ((1, 0), (2,))

    Returns
    -------
    clusters : tuple of tuples
        Each cluster is defined as a tuple of tuples, where the first
        element for each cluster is its centroid.
    """
    logger.debug("Trying to clusterize %d molecules." % len(fps))
    logger.debug("Defined cutoff: %.2f. Defined similarity function: %s."
                 % (cutoff, similarity_func))

    # first generate the distance matrix.
    dists = calc_distance_matrix(fps, similarity_func)
    logger.debug("Distance matrix created.")

    # now cluster the data.
    cs = Butina.ClusterData(dists, len(fps), cutoff, isDistData=True)
    logger.debug("Number of cluster(s) created: %d." % len(cs))

    return cs
