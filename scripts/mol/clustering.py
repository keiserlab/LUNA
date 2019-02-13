from rdkit import DataStructs
from rdkit.ML.Cluster import Butina
from util.exceptions import IllegalArgumentError

import re
import logging

logger = logging.getLogger()


def available_similarity_functions():
    regex = re.compile("Bulk([a-zA-Z]+)Similarity", flags=0)
    funcs = list(filter(regex.match, dir(DataStructs)))

    return funcs


def calc_distance_matrix(fps, similarity_func="BulkTanimotoSimilarity"):
    funcs = available_similarity_functions()
    if (similarity_func not in funcs):
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


def cluster_fps_butina(fps, cutoff=0.2,
                       similarity_func="BulkTanimotoSimilarity"):
    logger.info("Trying to clusterize %d molecules." % len(fps))
    logger.info("Defined cutoff: %.2f. Defined similarity function: %s."
                % (cutoff, similarity_func))

    # first generate the distance matrix:
    dists = calc_distance_matrix(fps, similarity_func)
    logger.info("Distance matrix created.")

    # now cluster the data:
    cs = Butina.ClusterData(dists, len(fps), cutoff, isDistData=True)
    logger.info("Number of cluster(s) created: %d." % len(cs))

    return cs
