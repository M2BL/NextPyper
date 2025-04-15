#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2025
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.
"""
Helper functions for sampling the most dissimilar paths in a graph component.
"""

__version__ = "0.1"


# =======================================================================================
#               IMPORTS
# =======================================================================================
from itertools import groupby

import kmedoids as km
import numpy as np

from gfa_graph import OrientedEdge, Assembly_graph

# =============================================================================
#                FUNCTIONS
# =============================================================================


def _jaccard_similarity(list1: list, list2: list) -> float:
    s1 = set(list1)
    s2 = set(list2)
    return float(len(s1.intersection(s2)) / len(s1.union(s2)))


def _compute_dist_matrix(paths: list[list[OrientedEdge]]):
    """
    Computes the distance matrix between paths using Jaccard similarity.
    """
    mx = np.zeros((len(paths), len(paths)))
    for i in range(len(paths)):
        for j in range(i + 1, len(paths)):
            edges_i = [x.id for x in paths[i]]
            edges_j = [x.id for x in paths[j]]
            distance = _jaccard_similarity(edges_i, edges_j)
            mx[i, j] = distance
    return mx + mx.T


def select_k_paths(
    paths: list[list[OrientedEdge]], K: int, graph: Assembly_graph
) -> list[list[OrientedEdge]]:
    """
    select K paths in order to sample as much edge diversity as possible.
    Paths are first clustered using Jaccard similarity as distance metric. Clustering algorithm is kmedoids.
    """

    def _select_longest_path(grouped_path: list[list[OrientedEdge]]):
        """
        Selects the longest path from a group of paths.
        """
        if len(grouped_path) == 1:
            return grouped_path
        length_paths = []
        for path in grouped_path:
            length = sum([len(graph.edge_dict[edge.id]) for edge in path])
            length_paths.append((length, path))
        return sorted(length_paths, key=lambda x: x[0], reverse=True)[0][1]

    if len(paths) <= K:
        return paths
    distmatrix = _compute_dist_matrix(paths)
    clusters = km.fasterpam(distmatrix, K)
    cluster_idx = sorted(zip(clusters.labels, paths), key=lambda x: x[0])
    groups = []
    for _, g in groupby(cluster_idx, lambda x: x[0]):
        groups.append([x[0] for x in g])
    best_paths = []
    for group in groups:
        best_paths.append(_select_longest_path(group))
    return best_paths


def main():
    import string
    from itertools import combinations, islice
    from random import randrange, sample

    np.set_printoptions(formatter={"all": lambda x: "{:.4g}".format(x)})
    letters = list(string.ascii_lowercase)
    edge_id = ["".join(x) for x in islice(combinations(letters, 3), 0, 20)]
    nbr_path = 10
    paths_ids = []
    while len(paths_ids) < nbr_path:
        path_length = randrange(3, 10)
        paths_ids.append(sample(edge_id, path_length))
    paths = []
    for path in paths_ids:
        paths.append([OrientedEdge(id, "+") for id in path])
    print(paths)

    distances = _compute_dist_matrix(paths)
    for row in distances:
        print(row)


if __name__ == "__main__":
    main()
