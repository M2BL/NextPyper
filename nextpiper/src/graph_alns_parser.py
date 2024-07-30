#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2024
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.
""" ## ToDo: Change module docstring
Functions and classes for parsing the gpa files.

of assembly graphs produced by SPAdes with the --custom-hmms flag.
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
from collections import namedtuple
from dataclasses import dataclass, field
from operator import itemgetter
from itertools import groupby
from pathlib import Path

OrientedEdge = namedtuple("OrientedEdge", ("id", "orientation"))


@dataclass
class EdgeAln:

    name: str
    edge: OrientedEdge
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
        edge = OrientedEdge(edge[:-1], edge[-1])

        return cls(fname, edge, *list(map(int, coords)), cigar)


@dataclass
class GraphAln:
    name: str
    fragments: list[EdgeAln]


@dataclass
class Read(GraphAln):
    mate: "Read" = field(init=False, default=None)

    def _set_mate(self, mate: "Read", pair_suffix: str | None = None) -> None:
        if self.name.rstrip(f"{pair_suffix}12") == mate.name.rstrip(f"{pair_suffix}12"):
            self.mate = mate
        else:
            raise ValueError(
                f"Read 'pairs' do not have the same name ({self.name=} and {mate.name})"
            )


def parse_gpa(
    gpa_path: Path, gpa2_path: Path | None = None, pair_suffix: str | None = None
) -> list[GraphAln | Read]:
    """Parse a given gpa and return a list of Graph alignments encoded in the given file.

    If a second gpa file is given, it is assumed that the files correspond to paired end
    reads (forward and reverse, respectively), and that they are ordered.

    Paired reads are assumed to be named identical unless pair_suffix is provided, in
    which case it will be used to trim the suffix used to distinguish the pairs.
    For instance if reads are named 'my_paired_reads/1' and 'my_paired_reads/2',
    pair_suffix is '/'.
    """

    # ToDo: Check how to properly annotate this function,
    # ToDo: It feels inappropriate because the input is a class and the output an *instance* of that class
    # A proposal (with generics):
    # def _parse_alns[C](file_path: Path, constructor: C) -> list[C]:

    # Second proposal:
    def _parse_alns[
        C: (GraphAln, Read)
    ](file_path: Path, constructor: type[C]) -> list[C]:

        # def _parse_alns(
        #         file_path: Path, constructor: GraphAln | Read
        #     ) -> list[GraphAln | Read]:

        get_aln = lambda line: line.split()[2]
        with file_path.open() as file:
            return [
                constructor(aln, [EdgeAln.from_line(line) for line in lines])
                for aln, lines in groupby(
                    (line for line in file if line.startswith("A")), key=get_aln
                )
            ]

    if gpa2_path is None:
        return _parse_alns(gpa_path, GraphAln)
    else:
        pairs_dict = {
            forward.name.removesuffix(f"{pair_suffix}1"): forward
            for forward in _parse_alns(gpa_path, Read)
        }

        for reverse in _parse_alns(gpa2_path, Read):
            if forward_mate := pairs_dict.get(
                reverse.name.removesuffix(f"{pair_suffix}2")
            ):
                forward_mate._set_mate(reverse, pair_suffix)
                reverse._set_mate(forward_mate, pair_suffix)

        return [read for read in pairs_dict.values() if read.mate is not None]
