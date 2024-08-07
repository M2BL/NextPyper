#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2024
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.
"""
Functions and classes for parsing the different outputs of SPAligner.
"""


__version__ = "0.1"


# =======================================================================================
#               IMPORTS
# =======================================================================================
from dataclasses import dataclass, field
from itertools import groupby
from pathlib import Path
from typing import Self

from gfa_graph import Path_on_graph
from graph_alns_parser import OrientedEdge
# =============================================================================
#                CLASSES
# =============================================================================


@dataclass
class Cigar:
    """
        Compute a distance out of a cigar string.
    Attributes
    ----------
    -cigar: raw string
    -items: list of value, operations in the string, e.g. '20H20M20S' gives [(20, 'H'), (20, 'M'), (20, 'S')].
    -distance: currently an uncorrected distance measure over aligned regions (no indels).
    """

    cigar: str
    items: list[tuple[int, str]] = field(default_factory=list, init=False)
    distance: float = field(default=1, init=False)

    def __post_init__(self):
        # query_consuming_ops = ('M', 'I', 'S', '=', 'X')
        # ref_consuming_ops = ('M', 'D', 'N', '=', 'X')
        self._itemize()
        if self.items:
            self._compute_distance()

    def _itemize(self) -> Self:
        if self.cigar != "*":
            cig_iter = groupby(self.cigar, lambda c: c.isdigit())
            self.items = [
                (int("".join(n)), "".join(next(cig_iter)[1])) for g, n in cig_iter
            ]
        return self

    def _compute_distance(self) -> Self:
        "Hamming distance, that do not take indels into account."
        ref_len = sum(l for l, op in self.items if op in ["M", "="])
        self.distance = sum(l for l, op in self.items if op == "X") / ref_len
        return self

    def get_distance(self) -> float:
        return self.distance


@dataclass
class Gpa_consensus:
    """
        Parse the output of SPAligner when matching a consensus sequence to the graph.
    Attributes
    ----------
    -gpa_filename: string indicating the path to the file.
    Post init
    ----------
    -matching_paths: list of Path_on_graph objects.
    """

    gpa_filename: str
    matching_paths: list[Path_on_graph] = field(default_factory=list, init=False)

    def __post_init__(self):
        self._parse_gpa(self.gpa_filename)

    def _parse_gpa(self, gpa_filename):
        gpa_path = Path(gpa_filename)
        assert gpa_path.exists(), f"GPA file {self.gpa_filename} not found"
        lines = gpa_path.read_text().splitlines()[1:]
        if not lines:
            return self
        target = lines[0].split("\t")
        pile = lines[1:]
        current_hit = target[2]
        start = int(target[3])
        length = int(target[4])
        end = start + length
        edges = [OrientedEdge(target[6][:-1], target[6][-1])]
        cigar = target[10]
        if not pile:
            new_path = Path_on_graph(
                current_hit, start, end, edges, length, Cigar(cigar).get_distance()
            )
            self.matching_paths.append(new_path)

        while pile:
            target = pile[0].split("\t")
            hit = target[2]
            end = int(target[3]) + int(target[4])
            if hit == current_hit:
                edges.append(OrientedEdge(target[6][:-1], target[6][-1]))
                length += int(target[4])
                cigar += target[10]
            else:
                new_path = Path_on_graph(
                    current_hit, start, end, edges, length, Cigar(cigar).get_distance()
                )
                self.matching_paths.append(new_path)
                current_hit = hit
                start = int(target[3])
                length = int(target[4])
                edges = [OrientedEdge(target[6][:-1], target[6][-1])]
                cigar = target[10]
            pile = pile[1:]
            if not pile:
                new_path = Path_on_graph(
                    current_hit, start, end, edges, length, Cigar(cigar).get_distance()
                )
                self.matching_paths.append(new_path)
        return self

    def get_matching_paths(self) -> list[Path_on_graph]:
        return self.matching_paths


# =============================================================================
#                FUNCTIONS
# =============================================================================
def main():
    c = Cigar("299=1X109=1X5=1X66=1X328=1X311=1X47=1X1=1X210=1X210=")
    gpa = Gpa_consensus(
        "/home/yjkbertrand/programs/Short-Pair/spaliner_gpa/alignment.gpa"
    )
    print(gpa.get_matching_paths())


if __name__ == "__main__":
    main()
