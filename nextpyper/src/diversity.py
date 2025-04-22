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
from typing import Callable, Optional, TYPE_CHECKING

import kmedoids as km
import numpy as np

if TYPE_CHECKING:
    from gfa_graph import OrientedEdge, Assembly_graph

# =============================================================================
#                FUNCTIONS
# =============================================================================


def _jaccard_similarity(list1: list, list2: list) -> float:
    s1 = set(list1)
    s2 = set(list2)
    return float(len(s1.intersection(s2)) / len(s1.union(s2)))


def _compute_dist_matrix(paths: list[list["OrientedEdge"]]):
    """
    Computes the distance matrix between paths using Jaccard similarity.
    """
    mx = np.eye(len(paths), dtype=np.float64) / 2
    for i in range(len(paths)):
        for j in range(i + 1, len(paths)):
            edges_i = [x.id for x in paths[i]]
            edges_j = [x.id for x in paths[j]]
            distance = _jaccard_similarity(edges_i, edges_j)
            mx[i, j] = distance
    return 1 - (mx + mx.T)


def select_k_paths(
    paths: list[list["OrientedEdge"]],
    K: int,
    graph: "Assembly_graph",
    key: Optional[Callable[[list["OrientedEdge"]], int]] = None,
) -> list[list["OrientedEdge"]]:
    """
    select K paths in order to sample as much edge diversity as possible.
    Paths are first clustered using Jaccard similarity as distance metric. Clustering algorithm is kmedoids.

    A custom path scoring function can be provided as key, where the highest scoring path
    will be chosen per cluster. The default function is length (longest path is chosen).
    """

    def _get_path_length(path: list["OrientedEdge"]) -> int:
        return sum(len(graph.edge_dict[edge.id]) for edge in path) - (
            graph.K * (len(path) - 1)
        )

    if len(paths) <= K:
        return paths
    distmatrix = _compute_dist_matrix(paths)
    clusters = km.fasterpam(distmatrix, K)
    cluster_idx = sorted(zip(clusters.labels, paths), key=lambda x: x[0])
    groups = []
    for _, g in groupby(cluster_idx, lambda x: x[0]):
        groups.append([x[1] for x in g])
    best_paths = []
    for group in groups:
        best_paths.append(max(group, key=_get_path_length if key is None else key))
    return best_paths


def main(): ...


if __name__ == "__main__":
    main()
