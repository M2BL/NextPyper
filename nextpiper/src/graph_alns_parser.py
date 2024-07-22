#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2024
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.
""" ## ToDo: Change module docstring
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
from dataclasses import dataclass, field
from operator import itemgetter
from itertools import groupby
from pathlib import Path
from typing import Literal


@dataclass
class EdgeAln:

    name: str
    edge_id: str
    orientation: Literal["+", "-"]
    query_start: int
    query_len: int
    target_start: int
    target_len: int
    cigar: str

    @property
    def query_end(self) -> int:
        return self.query_start + self.query_len

    @property
    def tend(self) -> int:
        return self.target_start + self.target_len

    @property
    def identity(self) -> float: ...

    @classmethod
    def from_line(cls, line: str) -> "EdgeAln":
        get_fields = itemgetter(1, 6, 10, 3, 4, 7, 8)
        fname, edge, cigar, *coords = get_fields(line.split())

        return cls(fname, edge[:-1], edge[-1], *list(map(int, coords)), cigar)


@dataclass
class GraphAln:
    name: str
    fragments: list[EdgeAln]


@dataclass
class Read(GraphAln):
    mate: "Read" = field(init=False)


def parse_gpa(gpa_path: Path) -> list[GraphAln]:

    get_aln = lambda line: line.split()[2]
    with gpa_path.open() as gpa:
        return [
            GraphAln(aln, [EdgeAln.from_line(line) for line in lines])
            for aln, lines in groupby(
                (line for line in gpa if line.startswith("A")), key=get_aln
            )
        ]
