#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2024
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.
"""
Functions and classes for parsing the gfa of assembly graphs produced by SPAdes with the --custom-hmms flag.
HMM profiles for probes have been used during the assembly, so that the 'hmm_statistics.txt' file from SPAdes
contains the graph edges that have a hmm match.
#  Usage example:
    gfa_file = ".../test_data/test_clustering/assembly_graph_after_simplification.gfa"
    hmm_statistics = "../test_data/test_clustering/hmm_statistics.txt"
    # parse the SPAdes output, perform the computation:
    components = filter_components_hmm(gfa_file, hmm_statistics)
    The output consist of a list of Component objects, that correspond to the assembly subgraph that have an HMM match.
"""

__version__ = "0.1"

# =======================================================================================
#               IMPORTS
# =======================================================================================
from collections import defaultdict, namedtuple
from dataclasses import dataclass, field
from itertools import chain
import os
from pathlib import Path
import sys
from typing import Self
from union_find import UnionFind


Component_hmm = namedtuple("Component_hmm", ("component", "hmm"))


@dataclass
class BGC_candidate:
    name: str
    hmms: list[str] = field(default_factory=list)
    coordinates: list[tuple[int, int]] = field(default_factory=list)
    edges: list[str] = field(default_factory=list)
    lengths: list[int] = field(default_factory=list, init=False)
    edge_to_hmms: dict[str, str] = field(default_factory=dict, init=False)
    max_domain_length: int = field(init=False)

    def __post_init__(self):
        self._get_lengths()
        self._find_dominant_hmm()

    def _get_lengths(self) -> Self:
        for cor in self.coordinates:
            self.lengths.append(cor[1] - cor[0])
        return self

    def _find_dominant_hmm(self) -> Self:
        hmm_length = defaultdict(int)
        for idx, length in enumerate(self.lengths):
            hmm_length[self.hmms[idx]] += length
        dominant_hmm = sorted(hmm_length.items(), key=lambda x: x[1], reverse=True)[0]
        self.max_domain_length = dominant_hmm[1]
        self.edge_to_hmms = {edge: dominant_hmm[0] for edge in self.edges}
        return Self

    def get_edge_to_hmm(self) -> dict[str, str]:
        return self.edge_to_hmms

    def get_max_length(self) -> int:
        return self.max_domain_length

    def __repr__(self):
        return f"BGC_candidate(name={self.name}, hmms={set(self.hmms)}, edges={self.edges} lengths={self.lengths}) "


@dataclass
class Component:
    edges: set[str]
    hmm: str

    def get_edges(self):
        return self.edges

    def get_hmm(self):
        return self.hmm


def components_from_gfa(gfa_file: str) -> list[set[str]]:
    """
    Split nodes in a gfa file per component and for each component return a set of node ids.
    :param gfa_file:
    :return:
    """
    gfa_file = Path(gfa_file)
    all_nodes = []
    edges: list[tuple] = []
    for line in open(gfa_file):
        if line.startswith("L"):
            splt_line = line.split("\t")
            edges.append((splt_line[1], splt_line[3]))
        elif line.startswith("S"):
            split_line = line.split("\t")
            all_nodes.append(split_line[1])
    with UnionFind() as UV:
        for edge in edges:
            UV.union(edge)
    components = UV.get_components()
    components_list = [list(component) for component in components]
    single_nodes = set(all_nodes).difference(set(list(chain(*components_list))))
    components.extend([{node} for node in single_nodes])
    return components


def matched_edges_from_hmm(hmm_stat_file: str, min_domain_len=20) -> dict[str, str]:
    """
    From a 'hmm_statistics.txt' file generated from spades with the --custom-hmms flag, retrieve the edges that
    have matched a hmm profile.
    :param hmm_stat_file:
    :param min_domain_len: minimum number of matched aa by the hmm profile over the whole length of a contig.
    :return: a dict with edge id as key and hmm id as value.
    """
    hmm_stat_file = Path(hmm_stat_file)
    lines = [line for line in open(hmm_stat_file)]
    bgc_candidates = []
    hmms = []
    domain_flag = False
    edge_flag = False
    for idx, line in enumerate(lines):
        if line.startswith("BGC subgraph"):
            if hmms:
                bgc_candidates.append(BGC_candidate(name, hmms, coordinates, edges))
            name = ""
            hmms = []
            coordinates = []
            edges = []
            name += line.split()[2]
        if line.startswith("BGC candidate"):
            name += f"_candidate_{line.split()[2]}"
            hmms = lines[idx + 1].strip().split("-")
            continue
        if line.startswith("Domain coordinates:"):
            domain_flag = True
            continue
        if line.startswith("Edge order:"):
            domain_flag = False
            edge_flag = True
            continue
        if line.startswith("Path"):
            edge_flag = False
            continue
        if domain_flag:
            splt_domain = line.strip().split()
            assert len(splt_domain) == 2, f"[Error] line {line} has the wrong format"
            coordinates.append(tuple(map(int, splt_domain)))
        if edge_flag:
            edges.extend(
                line.strip()
                .replace("-", "")
                .replace("+", "")
                .replace(";", "")
                .split(",")
            )
    edge_to_hmms = {}
    Hmm_length = namedtuple("Hmm_length", ["hmm", "length"])
    for bgc in bgc_candidates:
        bgc_edge_to_hmms = bgc.get_edge_to_hmm()
        max_length = bgc.get_max_length()
        for edge in bgc_edge_to_hmms:
            hmm = bgc_edge_to_hmms[edge]
            if edge in edge_to_hmms:
                length = edge_to_hmms[edge].length
                if max_length > length:
                    edge_to_hmms[edge] = Hmm_length(hmm, length)
            else:
                edge_to_hmms[edge] = Hmm_length(hmm, max_length)
    return {
        key: value.hmm
        for key, value in edge_to_hmms.items()
        if value.length > min_domain_len
    }


def filter_components_hmm(
    gfa_file, hmm_stat_file, min_domain_len=20
) -> list[Component]:
    all_components = []
    components = components_from_gfa(gfa_file)
    matched_edges = matched_edges_from_hmm(hmm_stat_file, min_domain_len)

    for component in components:
        hmm_matches = []
        for edge in component:
            if (hmm := matched_edges.get(edge)) is not None:
                hmm_matches.append(hmm)
        if hmm_matches:
            dominant_hmm = max(set(hmm_matches), key=hmm_matches.count)
            all_components.append(Component(component, dominant_hmm))
    return all_components


def main():
    ...
    # os.chdir(
    #     "/home/yjkbertrand/Documents/projects/nextpiper/test_data/gold_standards/brassica/mapping"
    # )
    # gfa = "assembly_graph_after_simplification.gfa"
    # # print(components_from_gfa(gfa))
    # hmm_stat = "hmm_statistics.txt"
    # # print(matched_edges_from_hmm(hmm_stat))
    # print(filter_components_hmm(gfa, hmm_stat))


if __name__ == "__main__":
    main()
