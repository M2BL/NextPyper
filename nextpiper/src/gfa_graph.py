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

from Bio.Seq import Seq

from graph_alns_parser import Read

oriented_edge = ("oriented_edge", tuple[str, Literal["+", "-"]])


# =============================================================================
#                CLASSES
# =============================================================================
@dataclass
class Edge:
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

    def retrieve_seq(self, start, end) -> Seq:
        return self.seq[start:end]


@dataclass
class Assembly_graph:
    gfa_filename: str
    K: int = field(init=False)
    edge_dict: dict[str, Edge] = field(default_factory=dict, init=False)
    graph: dict[oriented_edge, list[oriented_edge]] = field(
        default_factory=lambda: defaultdict(list), init=False
    )
    linked_edges: dict[tuple[str], int] = field(default_factory=lambda: defaultdict(int), init=False)
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
                        header = line
                    case "S":
                        edge_line = line.strip().split("\t")[1:]
                        node_id, seq, kc = edge_line[0], edge_line[1], edge_line[-1]
                        coverage = float(kc[len("KC:i:") :]) / len(seq)
                        self.edge_dict[node_id] = Edge(node_id, seq, coverage)

                    case "L" | "J":
                        _, node_id1, pos1, node_id2, pos2, match = line.strip().split(
                            "\t"
                        )
                        self.K = int(match[:-1])
                        self.graph[(node_id1, pos1)].append((node_id2, pos2))
                        self.graph[(node_id2, self.rev[pos2])].append(
                            (node_id1, self.rev[pos1])
                        )
                    case _:
                        raise NotImplementedError(
                            f"ERROR: found line of type {line[0]}"
                        )
    def link_edges(self, links:list[Read]) -> Self:
        return self

    def retrieve_path(self, start: int, end: int, path: list[oriented_edge]) -> Seq:
        ...


# =============================================================================
#                FUNCTIONS
# =============================================================================


def main():
    import os

    os.chdir(
        "/home/yjkbertrand/Documents/projects/nextpiper/test_data/gold_standards/brassica/spades/carinata/gold_standard_B_carinata_200_with_hmm_probe"
    )
    AG = Assembly_graph("assembly_graph_after_simplification.gfa")
    print(AG)


if __name__ == "__main__":
    main()
