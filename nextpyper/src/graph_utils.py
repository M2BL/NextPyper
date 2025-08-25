#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2024
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.
"""
Utility functions to support graph operations, coloring, extensions, filtering, etc.
"""

__version__ = "0.1"


# =============================================================================
#                CLASSES
# =============================================================================


from itertools import groupby
from operator import add, attrgetter, itemgetter
import re
from typing import Literal, NamedTuple, TYPE_CHECKING
from fractions import Fraction

from more_itertools import split_when
import polars as pl
from intervaltree import IntervalTree, Interval

if TYPE_CHECKING:
    from gfa_graph import Assembly_graph


class OrientedEdge(NamedTuple):
    id: str
    orientation: Literal["+", "-"]
    jump: bool = False

    @classmethod
    def from_seg(cls, seg: str) -> "OrientedEdge":
        if match := re.match(r",|;", seg):
            return cls(seg[1:-1], seg[-1], match[0] == ";")
        else:
            return cls(seg[:-1], seg[-1])


# Maybe move it
class ColorInfo(NamedTuple):
    probe: str
    tlen: int
    cis: bool
    matches: int
    mismatches: int

    @property
    def idt(self) -> Fraction:
        return Fraction(self.matches, (self.matches + self.mismatches))

    def __add__(self, other: "ColorInfo") -> "ColorInfo":
        if self.probe != other.probe:
            raise ValueError(
                "Cannot merge intervals of different probes: {self} and {other} )"
            )

        matches = self.matches + other.matches
        mismatches = self.mismatches + other.mismatches

        return ColorInfo(self.probe, self.tlen, self.cis, matches, mismatches)


# =============================================================================
#                FUNCTIONS
# =============================================================================


def probe_cov(
    path: list[OrientedEdge],
    graph: "Assembly_graph",
    no_reduce: bool = False,
) -> float | dict[str, float]:
    """Given a path and a colored graph, compute the probe coverages the path reaches.
    The result is a dictionary of coverage for each probe touch by the path.

    The default behaviour is to return the coverage of the probe most covered (max)."""

    edges = [graph.edge_dict[edge.id] for edge in path]
    probes = {probe for edge in edges for probe in edge.get_colors()}

    if len(probes) == 0:
        return 0.0

    probe_covs = {}
    for probe in probes:
        tree = IntervalTree(
            inter for edge in edges for inter in edge.matching_exons[probe]
        )
        if tree.is_empty():
            probe_covs[probe] = 0.0
        else:
            tree.merge_overlaps(data_reducer=add)
            cov_bases = sum(interval.length() for interval in tree)
            tlen = next(iter(tree)).data.tlen
            probe_covs[probe] = cov_bases / tlen

    if no_reduce:
        return probe_covs
    else:
        return max(probe_covs.values())


def effective_cov(probe_tree: IntervalTree, colored_intervals: bool = False) -> int:
    "Compute the total number of matches on a probe tree."
    if colored_intervals:
        return round(sum(inter.length() * inter.data.idt for inter in probe_tree))
    else:
        return round(sum(inter.length() * inter.data for inter in probe_tree))


def merge_hits(hits: list[Interval]) -> Interval:
    "Given two mergeable hits, return a single Merged interval."
    return Interval(hits[0].begin, hits[-1].end, hits[0].data)


def not_mergeable(hit1: Interval, hit2: Interval) -> bool:
    "Return whether to hits are compatible for merging (adjacent and same identity)"
    return hit1.distance_to(hit2) != 0 or hit1.data != hit2.data


def build_probe_trees(
    df: pl.DataFrame, min_idt: float = 0.7
) -> dict[str, IntervalTree]:
    """Given a Dataframe of probe hits, for each probe build a probe tree with the
    coverage of all the hits, while keeping only the best hit in each region. Hits
    with similarity below min_idt are ignored."""

    TREE_COLS = ["theader", "tstart", "tend", "nident", "mismatch"]
    probe_hits = df.sort("theader").select(TREE_COLS)

    probe_trees = {}
    for probe, phits in groupby(probe_hits.iter_rows(), itemgetter(0)):
        tree = IntervalTree(
            Interval(tstart, tend, idt)
            for _, tstart, tend, n, m in phits
            if (idt := Fraction(n, n + m)) >= min_idt
        )
        tree.split_overlaps()

        best_hits = [
            max(homologous_hits, key=attrgetter("data"))
            for _, homologous_hits in groupby(sorted(tree), key=attrgetter("begin"))
        ]

        probe_trees[probe] = IntervalTree(
            map(merge_hits, split_when(best_hits, not_mergeable))
        )

    return probe_trees


def filt_probe_hits(
    df: pl.DataFrame, probe_trees: dict[str, IntervalTree]
) -> pl.DataFrame:
    """Compute the effective coverage (matches) on the probes and pick
    the best probe based on it.

    Return a filtered dataframe with only the hits of the best probe.
    """

    PROBE_COLS = ["theader", "tprobe", "tlen"]

    # Compute the "Effective coverage" on the probes
    probes_cov = (
        df.select(PROBE_COLS)
        .unique()
        .with_columns(
            glob_eff_cov=pl.col("theader").map_elements(
                lambda probe: effective_cov(probe_trees[probe]),
                return_dtype=pl.Int64,
            ),
            glob_cov=pl.col("theader").map_elements(
                lambda probe: sum(inter.length() for inter in probe_trees[probe]),
                return_dtype=pl.Int64,
            ),
        )
    )

    # Pick the best probe version
    best_probe_ver = probes_cov.group_by("tprobe").agg(
        pl.all().sort_by(["tprobe", "glob_eff_cov"], descending=True).first()
    )

    # Keep only the hits of the best probes
    return df.join(best_probe_ver, on=["theader", "tprobe", "tlen"])
