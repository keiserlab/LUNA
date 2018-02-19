def available_similarity_functions():
    import re
    from rdkit import DataStructs

    regex = re.compile("Bulk([a-zA-Z]+)Similarity", flags=0)
    funcs = list(filter(regex.match, dir(DataStructs)))

    return funcs


def calc_distance_matrix(fps, similarityFunc="BulkTanimotoSimilarity"):
    from rdkit import DataStructs
    from util.exceptions import IllegalArgumentError

    funcs = available_similarity_functions()
    if (similarityFunc not in funcs):
        raise IllegalArgumentError("Similarity function not available.")

    dists = []
    for i in range(1, len(fps)):
        if (similarityFunc == "BulkTverskySimilarity"):
            params = [fps[i], fps[:i], 0, 1]
        else:
            params = [fps[i], fps[:i]]

        sims = getattr(DataStructs, similarityFunc)(*params)
        dists.extend([1 - x for x in sims])

    return dists


def cluster_fps_butina(fps, cutoff=0.2,
                       similarityFunc="BulkTanimotoSimilarity"):
    from rdkit.ML.Cluster import Butina
    import logging
    logger = logging.getLogger(__name__)

    logger.info("Trying to clusterize %d molecules." % len(fps))
    logger.info("Defined cutoff: %.2f. Defined similarity function: %s."
                % (cutoff, similarityFunc))

    # first generate the distance matrix:
    dists = calc_distance_matrix(fps, similarityFunc)
    logger.info("Distance matrix created.")

    # now cluster the data:
    cs = Butina.ClusterData(dists, len(fps), cutoff, isDistData=True)
    logger.info("Number of cluster(s) created: %d." % len(cs))

    return cs
