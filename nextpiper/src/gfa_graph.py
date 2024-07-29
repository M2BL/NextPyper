#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2024
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.
"""
Functions and classes for parsing the gfa assembly graph and retrieve path sequence.
"""

__version__ = "0.1"


# =======================================================================================
#               IMPORTS
# =======================================================================================
from collections import defaultdict, namedtuple
from dataclasses import dataclass, field
from operator import itemgetter, attrgetter
from itertools import chain, groupby
from pathlib import Path
from typing import Self, NewType, List, Dict, Tuple, Optional, Literal
import sys

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from graph_alns_parser import Read

oriented_edge = NewType("oriented_edge", tuple[str, Literal["+", "-"]])


# =============================================================================
#                CLASSES
# =============================================================================
@dataclass(frozen=True)
class Path_on_graph:
    """
        Encodes the path on the graph matching either a sequence (SPAligner) or a hmm profile (PathRacer).
    Attributes
    ----------
    -start: coordinate of the starting position on the first edge.
    -end: coordinate of the ending position on the last edge.
    -edges: list of 'oriented_edge' tuples.
    """

    start: int
    end: int
    edges: list[oriented_edge]

    def get_parameters(self):
        return self.edges, self.start, self.end


@dataclass
class Edge:
    """
        Encodes the edge of the graph, i.e. the 'S' line in a gfa file.
        Keeps information about matching elements, that is sequences from SPAligner or hmm profile from PathRacer.
    Attributes
    ----------
    -id: edge number.
    -raw_seq: nucleotide sequence including k-mer overlap.
    -coverage: k-mer coverage.

    Post init
    ----------
    -seq: 'raw_seq' converted to Bio.Seq.Seq object.
    -matching_exons: dictionary of matches (hmm or sequences, needs to be worked out).
    """

    id: str
    raw_seq: str
    coverage: float
    seq: Seq = field(init=False)
    matching_exons: dict[str, int] = field(
        init=False, default_factory=list
    )  # id, length and identity of matching contigs

    def __post_init__(self):
        self.seq = Seq(self.raw_seq)

    def has_match(self) -> bool:
        return bool(self.matching_exons)

    def get_seq(self) -> Seq:
        return self.seq

    def get_length(self) -> int:
        return len(self.seq)

    def retrieve_seq(self, start: int, end: int, orientation: Literal["+", "-"]) -> Seq:
        tmp_seq = self.seq[start:end]
        if orientation == "+":
            return tmp_seq
        return tmp_seq.reverse_complement()


@dataclass
class Assembly_graph:
    """
        Encodes the gfa assembly graph. Provides methods for graph traversal and sequence retrieval.
    Attributes
    ----------
    -gfa_filename: complete path to gfa assembly file.
    Post init
    -K: k-mer used during SPAdes assembly (edge overlap value).
    -edge_dict: edge number mapping to their respective Edge objects. Used for sequence retrieval.
    -graph: encoding of the graph structure, with an edge mapping.
        Each node represented symmetrically as (1,+)->(2,+) and (2,-)->(1,-).
    -linked_edges: information about non-adjacent edges, that are connected through mate-reads.
    ----------
    """

    gfa_filename: str
    K: int = field(init=False)
    edge_dict: dict[str, Edge] = field(default_factory=dict, init=False)
    graph: dict[oriented_edge, list[oriented_edge]] = field(
        default_factory=lambda: defaultdict(list), init=False
    )
    linked_edges: dict[tuple[str], int] = field(
        default_factory=lambda: defaultdict(int), init=False
    )
    rev = {"+": "-", "-": "+"}

    def __post_init__(self):
        self._parse_graph()

    def _parse_graph(self) -> Self:
        graph_path = Path(self.gfa_filename)
        assert graph_path.exists(), f"Graph file {self.gfa_filename} not found"
        with open(self.gfa_filename, "r") as file:
            for line in file:
                match line[0]:
                    case "H":
                        pass
                    case "S":
                        edge_line = line.strip().split("\t")[1:]
                        node_id, seq, _kc = edge_line[0], edge_line[1], edge_line[-1]
                        coverage = float(_kc[len("KC:i:") :]) / len(seq)
                        self.edge_dict[node_id] = Edge(node_id, seq, coverage)

                    case "L" | "J":
                        _, node_id1, pos1, node_id2, pos2, _match = line.strip().split(
                            "\t"
                        )
                        self.K = int(_match[:-1])
                        self.graph[(node_id1, pos1)].append((node_id2, pos2))
                        self.graph[(node_id2, self.rev[pos2])].append(
                            (node_id1, self.rev[pos1])
                        )
                    case _:
                        raise NotImplementedError(
                            f"ERROR: found line of type {line[0]}"
                        )

    def link_edges(self, links: list[Read]) -> Self:
        return self

    def _retrieve_path(
        self,
        path: list[oriented_edge],
        start: int = 0,
        end: int = -1,
        extension: str = "",
        first_edge: bool = True,
    ) -> Seq:
        """

        :param path:
        :param start:
        :param end:
        :param extension:
        :param first_edge:
        :return:
        """
        id = path[0][0]
        orientation = path[0][1]
        if len(path) == 1:
            if (edge := self.edge_dict[id]) is not None:
                if first_edge:
                    extension += edge.retrieve_seq(start, end, orientation)
                else:
                    extension += edge.retrieve_seq(0, end, orientation)[self.K :]
                return extension
            else:
                sys.exit(f"failed to find edge in graph {path}")

        else:
            if (edge := self.edge_dict[id]) is not None:
                if first_edge:
                    extension += edge.retrieve_seq(start, -1, orientation)
                else:
                    extension += edge.retrieve_seq(0, -1, orientation)[self.K :]
                return self._retrieve_path(path[1:], start, end, extension, False)
            else:
                sys.exit(f"failed to find edge in graph {path}")

    def retrieve_path(self, name: str, edge_paths: Path_on_graph) -> SeqRecord:
        """
        Retrieve the sequence corresponding to a graph transversal.
        :param name: Name of the sequence, should correspond to the name of the consensus sequence or hmm profile.
        :param edge_paths: Path_on_graph object
        :return:
        """
        return SeqRecord(
            seq=self._retrieve_path(*edge_paths.get_parameters()),
            id=name,
            description="",
            name="",
        )


# =============================================================================
#                FUNCTIONS
# =============================================================================


def main():
    ...


if __name__ == "__main__":
    main()
