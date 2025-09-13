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
from collections import Counter, defaultdict
from fractions import Fraction
from operator import attrgetter, itemgetter, add
from dataclasses import dataclass, field
from itertools import chain, count, repeat, groupby
from pathlib import Path
from typing import Callable, Iterator, Self, Literal, Optional, NamedTuple
from functools import partial
import re
import sys

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Align import PairwiseAligner
import polars as pl
from intervaltree import IntervalTree, Interval
from more_itertools import split_at, transpose

from graph_utils import (
    OrientedEdge,
    ColorInfo,
    filt_probe_hits,
    probe_cov,
    effective_cov,
    build_probe_trees,
)

# from graph_alns_parser import Read
from union_find import UnionFind
from diversity import select_k_paths

# =============================================================================
#                CONSTANTS
# =============================================================================

NODE_ID_PAT = r"^NODE_(\d+)_"
PROBE_PAT = r"-(.*?)_EDGE"
LOG_COLS = ["orig_name", "ext_name", "comp_id", "probe", "orig_path", "ext_path"]
DMND_COLS = [
    "query",
    "evalue",
    "qstart",
    "qend",
    "qlen",
    "tstart",
    "tend",
    "tlen",
    "theader",
    "gapopen",
    "nident",
    "mismatch",
]

# =============================================================================
#                CLASSES
# =============================================================================
LinkSupport = dict[tuple[str], int]


@dataclass(slots=True)
class Path_on_graph:
    """
    Sequence of OrientedEdges that compose a scaffold
    Attributes
    ----------
    -name: scaffold's name.
    -edge:
    -start: possible start on the first OrientedEdge.
    -end: possible end on the last OrientedEdge.
    -length: length of the scaffold in nucleotides. To be implemented.
    """

    name: str
    edges: list[OrientedEdge]
    start: Optional[int] = field(default=0)
    end: Optional[int] = field(default=None)
    length: Optional[int] = field(default=None)

    def get_parameters(self):
        return self.name, self.edges, self.start, self.end, self.length


class ExtLimits(NamedTuple):
    """Parameters to control the maximum extension of a path during
    graph exploration. The length limit is calculated by scaling the
    probe length by scaling. The resulting number is constrained to
    the range defined by [floor, ceiling].
    """

    floor: int
    ceiling: int
    scaling: float


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
    raw_seq: str = field(repr=False)
    coverage: float
    seq: Seq = field(init=False)
    matching_exons: dict[str, int] = field(
        init=False, default_factory=lambda: defaultdict(list)
    )  # name (key), interval (value) and extra info of the matching probe

    def __post_init__(self):
        self.seq = Seq(self.raw_seq)

    def __len__(self) -> int:
        return len(self.seq)

    def has_match(self) -> bool:
        return self.matching_exons and any(self.matching_exons.values())

    def get_colors(self) -> list[str]:
        return (
            [probe for probe, inter in self.matching_exons.items() if inter]
            if self.has_match()
            else list()
        )

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
        Each node represented symmetrically as (1,+)->[(2,+),(3,+)] and (2,-)->[(1,-)].
    -linked_edges: information about non-adjacent edges, that are connected through mate-reads.
    ----------
    """

    gfa_filename: str
    K: int = field(init=False)
    edge_dict: dict[str, Edge] = field(default_factory=dict, init=False)
    graph: dict[OrientedEdge, list[OrientedEdge]] = field(
        default_factory=lambda: defaultdict(list), init=False
    )
    linked_edges: LinkSupport = field(
        default_factory=lambda: defaultdict(int), init=False
    )
    paths: dict[str, Path_on_graph] = field(default_factory=dict, init=False)
    components: list[set[str]] = field(default_factory=list, init=False)
    rev = {"+": "-", "-": "+"}

    def __post_init__(self):
        self._parse_graph()

    def _parse_graph(self) -> Self:
        graph_path = Path(self.gfa_filename)
        assert graph_path.exists(), f"Graph file {self.gfa_filename} not found"
        with open(self.gfa_filename, "r") as file, UnionFind() as UV:
            seg_pat = re.compile(r"[,|;]?\d+[+-]")
            for line in file:
                match line[0]:
                    case "H":
                        pass
                    case "S":
                        edge_line = line.strip().split("\t")[1:]
                        node_id, seq, _kc = edge_line[0], edge_line[1], edge_line[-1]
                        coverage = float(_kc[len("KC:i:") :]) / len(seq)
                        self.edge_dict[node_id] = Edge(node_id, seq, coverage)

                    case "L":
                        _, node_id1, pos1, node_id2, pos2, _match = line.strip().split(
                            "\t"
                        )
                        self.K = int(_match[:-1])
                        self.graph[OrientedEdge(node_id1, pos1)].append(
                            OrientedEdge(node_id2, pos2)
                        )
                        self.graph[OrientedEdge(node_id2, self.rev[pos2])].append(
                            OrientedEdge(node_id1, self.rev[pos1])
                        )
                        UV.union((node_id1, node_id2))
                    case "J":
                        _, node_id1, pos1, node_id2, pos2, _dist, *_tags = (
                            line.strip().split("\t")
                        )
                        self.graph[OrientedEdge(node_id1, pos1)].append(
                            OrientedEdge(node_id2, pos2, True)
                        )
                        self.graph[OrientedEdge(node_id2, self.rev[pos2])].append(
                            OrientedEdge(node_id1, self.rev[pos1], True)
                        )

                        ## Add to the Jump Links (J-Lines) the connections from L-Lines to prevent graph disconnection
                        self.graph[OrientedEdge(node_id2, pos2, True)].extend(
                            self.graph[OrientedEdge(node_id2, pos2)]
                        )
                        self.graph[OrientedEdge(node_id1, self.rev[pos1], True)].extend(
                            self.graph[OrientedEdge(node_id1, self.rev[pos1])]
                        )
                        UV.union((node_id1, node_id2))

                    case "P":
                        _, name, path, _, *tags = line.split()
                        edges = [
                            OrientedEdge.from_seg(seg) for seg in seg_pat.findall(path)
                        ]
                        self.paths[name] = Path_on_graph(name, edges)
                    case _:
                        raise NotImplementedError(
                            f"ERROR: found line of type {line[0]}"
                        )
            # For ends here:
            components = UV.get_components()
            comp_list = list(map(list, UV.get_components()))
            single_edges = set(self.edge_dict).difference(set(chain(*comp_list)))
            components.extend([{edge} for edge in single_edges])
            self.components = dict(enumerate(components, 1))

        return self

    # def link_edges(self, reads: list["Read"]) -> Self:
    #     """Given a list of paired reads, compute the link support that those pairs
    #     exhibit. Each pair that connects a composition of edges is added as a
    #     supported link, increasing its count in the dictionary. Compositions of
    #     a single edge are not taken into account.
    #     """

    #     get_edges = lambda frags: (frag.edge.id for frag in frags)

    #     for read in reads:
    #         edges = get_edges(read.fragments + read.mate.fragments)
    #         if len(unique_edges := set(edges)) > 1:
    #             self.linked_edges[tuple(sorted(unique_edges))] += 1

    #     return self

    def path_support(self, path: Path_on_graph) -> LinkSupport:
        "Return the links support that are congruent with the given path."

        if not self.linked_edges:
            raise ValueError("No link information in the graph to evaluate path.")

        path_edges = {edge.id for edge in path.edges}
        return {
            link: support
            for link, support in self.linked_edges.items()
            if path_edges.issuperset(set(link))
        }

    @staticmethod
    def _ovlp_extension(
        seq1: Seq, seq2: Seq, gap: int = 10, size_only: bool = False
    ) -> Seq | int:
        """Use pairwise alignment to detect if seq1 and seq2 overlap. If they do
        Use the alignment coordinates to merge them. Otherwise, scaffold them
        with a gap of the given size.

        If size_only is True, instead of scaffolding the given sequences, return
        the size of the overlap or the size of the gap they produce."""

        aligner = PairwiseAligner(scoring="megablast", mode="local")
        alns = aligner.align(seq1, seq2)

        # Proceed as no overlap between the sequences. Add a gap of 10 N
        if alns.score < 30:
            return -1 * gap if size_only else seq1 + Seq("N" * gap) + seq2

        # Overlap detected. Use alignment coordinates to merge the sequences
        else:
            coords = alns[0].coordinates
            if size_only:
                return coords[0, -1] - coords[0, 0]
            else:
                end1, start2 = coords[:, -1]
                return seq1[:end1] + seq2[start2:]

    def _retrieve_path(
        self,
        path: list[OrientedEdge],
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
        id = path[0].id
        orientation = path[0].orientation
        jump = path[0].jump

        if len(path) == 1:
            if (edge := self.edge_dict[id]) is not None:
                if first_edge:
                    extension += edge.retrieve_seq(start, end, orientation)
                elif jump:
                    newseq = edge.retrieve_seq(0, end, orientation)
                    extension = self._ovlp_extension(extension, newseq)
                else:
                    extension += edge.retrieve_seq(0, end, orientation)[self.K :]
                return extension
            else:
                sys.exit(f"failed to find edge in graph {path}")

        else:
            if (edge := self.edge_dict[id]) is not None:
                if first_edge:
                    extension += edge.retrieve_seq(start, None, orientation)
                elif jump:
                    newseq = edge.retrieve_seq(0, None, orientation)
                    extension = self._ovlp_extension(extension, newseq)
                else:
                    extension += edge.retrieve_seq(0, None, orientation)[self.K :]
                return self._retrieve_path(path[1:], start, end, extension, False)
            else:
                sys.exit(f"failed to find edge in graph {path}")

    def retrieve_path(self, path: Path_on_graph) -> SeqRecord:
        """
        Retrieve the sequence corresponding to a graph transversal.
        :param name: Name of the sequence, should correspond to the name of the consensus sequence or hmm profile.
        :param path: Path_on_graph object
        :return:
        """
        return SeqRecord(
            seq=self._retrieve_path(*path.get_parameters()[1:4]),
            id=path.name,
            description="",
            name="",
        )

    def _get_path_intervals(self, path: str) -> list[Interval]:
        """Given a path name, return a list of intervals specifying which edge that part of the
        sequence comes from."""

        def _get_ovlp_size(
            edge1: OrientedEdge, edge2: OrientedEdge, gap: int = 10
        ) -> int:
            """Determine the size of the overlp between the given edges.
            If not overlap is found, return the gap size"""

            seq1 = self.edge_dict[edge1.id].seq
            seq2 = self.edge_dict[edge2.id].seq
            return self._ovlp_extension(seq1, seq2, gap, size_only=True)

        edges: list[OrientedEdge] = self.paths[path].edges
        edge, *edges = edges
        length = len(self.edge_dict[edge.id])
        start, end = 0, length
        intervals = [Interval(start, end, edge.id)]

        while edges:
            ovlp = _get_ovlp_size(edge, edges[0]) if edges[0].jump else self.K
            edge, *edges = edges
            start += length - ovlp
            length = len(self.edge_dict[edge.id])
            end += length - ovlp
            intervals.append(Interval(start, end, edge.id))

        return intervals

    def color_edges(self, hits: pl.DataFrame, min_idt: float = 0.7) -> Self:
        """Use a table of probe hits to color the edges of the graph, populating the attribute
        'matching exons' of each edge with the interval they match of the best probe of
        the graph's component they belong.

        If the identity of a hit is below min_idt, the hit is skipped and the htting edges
        are not colored.
        """

        KEEP_COLS = (
            "query",
            "theader",
            "cis",
            "tlen",
            "qstart",
            "qend",
            "tstart",
            "tend",
            "nident",
            "mismatch",
        )

        # ToDo: Should we also degroup this and keep only smaller more granular
        # hits? How would that afect the coloring? Compare results
        df = hits.group_by(["query", "theader", "cis"]).agg(
            pl.sum("nident"),
            pl.sum("mismatch"),
            pl.sum("gapopen"),
            pl.first("tlen"),
            pl.first("comp_id"),
            pl.first("tprobe"),
            pl.col("qstart"),
            pl.col("qend"),
            pl.col("tstart"),
            pl.col("tend"),
        )

        for (
            path,
            probe,
            direc,
            tlen,
            qstarts,
            qends,
            tstarts,
            tends,
            matches,
            mismatches,
        ) in df[KEEP_COLS].iter_rows():
            # Skip hits that have a low identity
            colorinfo = ColorInfo(probe, tlen, direc, matches, mismatches)
            if colorinfo.idt < min_idt:
                continue

            # If the hits are not in the same strand, swap starts and ends accordingly
            if not direc:
                qstarts, qends = qends, qstarts

            path_intervals = self._get_path_intervals(path)
            qinter = sorted(map(Interval, qstarts, qends))
            mod_tstarts = (tstart - 1 for tstart in tstarts)
            tinter = sorted(map(Interval, mod_tstarts, tends, repeat(colorinfo)))

            for qint, tint in zip(qinter, tinter):
                for pint in path_intervals:
                    if qint.overlaps(pint):
                        edge_id = pint.data
                        self.edge_dict[edge_id].matching_exons[probe].append(tint)

        # Simplify (merge) overlapping intervals of the same color.
        for edge in self.edge_dict.values():
            if edge.has_match():
                for probe, intervals in edge.matching_exons.items():
                    if len(intervals) > 1:
                        tree = IntervalTree(intervals)
                        tree.merge_overlaps(data_reducer=add, strict=False)
                        edge.matching_exons[probe] = sorted(tree)


@dataclass(slots=True)
class OptimalExtension:
    """
    Container of partial extensions of paths discovered during graph exploration.
    The container holds the N best extensions where N is the size of the container
    at instantiation. A function key function can be given to control how the "best"
    extensions are determined.

    extensions: a list containing the paths.
    size: How many extensions to hold.
    key: Function to score the extensions, which should receive a path as only argument.
    """

    extensions: list[list[OrientedEdge]] = field(default_factory=list)
    size: int = field(default=10)
    key: Optional[Callable[[list[OrientedEdge]], int]] = field(default=None)

    def __len__(self) -> int:
        return len(self.extensions)

    def __repr__(self) -> str:
        if self.key is None:
            return f"OptimalExtension(size={self.size}, extensions={self.extensions})"
        else:
            ranks = list(map(self.key, self.extensions))
            return f"OptimalExtension(size={self.size}, ranks={ranks})"

    def __iter__(self) -> Iterator[list[OrientedEdge]]:
        for path in iter(self.extensions):
            yield path

    def _purge(self) -> None:
        while len(self) > self.size:
            self.extensions.pop()

    def add(self, path: list[OrientedEdge]) -> None:
        self.extensions.append(path)
        self.extensions.sort(reverse=True, key=self.key)
        self._purge()
        return Self

    def extend(self, paths: list[list[OrientedEdge]]) -> None:
        self.extensions.extend(paths)
        self.extensions.sort(reverse=True, key=self.key)
        self._purge()
        return Self


# =============================================================================
#                FUNCTIONS
# =============================================================================


def dfs_track_paths(
    graph: Assembly_graph,
    start: OrientedEdge,
    probe: Optional[str] = None,
    max_len: int = 5000,
    max_extensions: int = 10,
    max_intron_size: Optional[int] = 2000,
    goal: Optional[OrientedEdge] = None,
    key: Optional[Callable[[list[OrientedEdge]], int]] = None,
):
    """
    Given a starting node (from a path), compute and return compatible with the given graph.

    The number of returned extensions is limited to the "best" n, where n is controlled
    by  _max_extensions_. The definition of "Best" can be customised by providing any
    scoring function that accepts a path and returns a score. By default, length is
    used (the longest path is the best).

    Finally, the exploration of the graph is limited to max_len. Extensions longer than
    this value are stopped early.

    Attributes
    ----------
    -graph: Assembly_graph
    -start: OrientedEdge last edge on scaffold to be extended
    -max_len: in nucleotides
    -max_extensions: max number of paths
    -goal: edge where to stop the DFS path exploration
    -key: function to score the produced paths. The function has
    to receive a single argument which has to be a path: list[OrientedEdge]
    """

    def get_path_len(path: list[OrientedEdge], graph: Assembly_graph) -> int:
        return sum(len(graph.edge_dict[edge[0]]) for edge in path) - graph.K * (
            len(path) - 1
        )

    def exons_gap(path: list[OrientedEdge], probe: str, graph: Assembly_graph) -> int:
        "Determine the maximum gap (in nucleotides) between colored edges in a path"
        if len(path) == 1:
            return 0

        def local_gap(edge: Edge, probe: str) -> int:
            if edge.has_match() and any(edge.matching_exons[probe]):
                return 0
            else:
                return len(edge)

        gaps = (local_gap(graph.edge_dict[edge.id], probe) for edge in path)
        return max(map(sum, split_at(gaps, lambda gap: gap == 0)))

    def dfs_helper(
        edge: OrientedEdge,
        current_path: list[OrientedEdge],
        extensions: OptimalExtension,
    ):
        """
        Help function to apply the DFS recursion.

        edge: current edge being evaluated.
        current_path: path (ordered list of edges being explored).
        extensions: container with the best paths.
        """
        current_path.append(edge)

        # If we reach the goal, add the path
        if goal is not None and edge == goal:
            extensions.add(current_path[:])
            return

        # Maximum len size check (avoids long recursions)
        if max_len is not None and get_path_len(current_path, graph) > max_len:
            extensions.add(current_path[:])
            return

        # Colored graph exploration
        if probe:
            # No-mix color policy (abort a path if changing color)
            if (
                edge_colors := set(graph.edge_dict[edge.id].get_colors())
            ) and probe not in edge_colors:
                extensions.add(current_path[:-1])
                return

            # Gap size is bigger than the limit
            if (
                max_intron_size
                and exons_gap(current_path, probe, graph) > max_intron_size
            ):
                extensions.add(current_path[:])
                return

        # Reached a dead-end
        if not (neighbors := graph.graph[edge]):
            extensions.add(current_path[:])
            return

        # Explore neighbors
        for neighbor in neighbors:
            if neighbor.id not in [edge.id for edge in current_path]:
                dfs_helper(neighbor, current_path[:], extensions)
            # Avoid getting stuck in loops and reversing your steps.
            else:
                extensions.add(current_path[:])
                return

    if key is None:
        key = partial(get_path_len, graph=graph)

    current_path = []
    extensions = OptimalExtension(size=max_extensions, key=key)

    dfs_helper(start, current_path, extensions)

    return extensions


def extend_path(
    path: Path_on_graph,
    graph: Assembly_graph,
    max_len: int = 5000,
    max_ext: int = 10,
    max_intron_size: int = 2000,
    key: Optional[Callable[[list[OrientedEdge]], int]] = None,
    allow_gray: bool = False,
    probe: str | None = None,
) -> list[OrientedEdge]:
    """Given an assembly graph and a path, extend the path following the graph topology
    using a Depth First Search. Up to max_ext extensions of the given path are returned.

    By default, the extension is guided by the color of the path. If the path is colored by
    more than one probe, the probe with the highest effective coverage (matches) is chosen.
    If the path is gray (colorless) and allow_gray is True, colorless extension follows.
    Otherwise a ValueError is performed.

    To avoid long recursions, an extension bigger than max_len stops the exploration of that
    extension. Path extensions are limited to max_len.

    Color guided extension is further limited by max_intron_size, which is the maximum gap
    allowed between colored edges. If an extension would introduce a gap bigger than this
    value, the extension is terminated.

    Since a path can have multiple possible extensions and these are limited to max_ext,
    The extensions are scored using the key funtion and onyl the best are kept. The key
    function should receive a path and return a score. By the default, the length of the
    path is used (longer paths are better).

    Return a list with all the extended paths, represented as list of oriented edges.
    """

    edges = [graph.edge_dict[edge.id] for edge in path.edges]
    probes = {probe for edge in edges for probe in edge.get_colors()}

    # A probe was provided to use for the exploration
    if probe is not None and probe not in probes:
        raise ValueError(f"path {path.name} is not colored by {probe} but by {probes}.")

    # Probe was not provided but it is going to be discovered for colored exploration
    elif probe is None:
        match len(probes):
            case 0:
                probe = None
            case 1:
                probe = probes.pop()
            case _:
                probe_idts = {}
                for probe in probes:
                    tree = IntervalTree(
                        inter for edge in edges for inter in edge.matching_exons[probe]
                    )
                    probe_idts[probe] = effective_cov(tree, colored_intervals=True)

                probe = max(probe_idts.items(), key=itemgetter(1))[0]

    if not allow_gray and probe is None:
        raise ValueError(f"All edges in path {path.name} are colorless. Exiting.")

    # In case of colored extension, deactivate the max_len filter.
    if probe is not None:
        max_len = None

    *edges, tail = path.edges
    extensions = dfs_track_paths(
        graph, tail, probe, max_len, max_ext, max_intron_size, key=key
    )
    return [edges + ext for ext in extensions]


def get_seq_atts(path: list[OrientedEdge], graph: Assembly_graph) -> tuple[int, float]:
    """Given a list of OrientedEdges which represent a path, infer the length and coverage of
    the sequence it encodes. Returns a tuple with both values.
    """

    edges = [graph.edge_dict[edge.id] for edge in path]
    path_length = sum(len(edge) for edge in edges) - (graph.K * (len(path) - 1))
    path_cov = sum(edge.coverage * len(edge) for edge in edges) / (
        path_length - graph.K
    )

    return path_length, path_cov


def make_path_name(path: list[OrientedEdge], idx: int, graph: Assembly_graph):
    length, cov = get_seq_atts(path, graph)
    return f"EDGE_{idx}_length_{length}_cov_{cov:.3f}"


def find_best_probe_hits(
    df: pl.DataFrame, min_idt: float = 0.7, div_mod: float = 0.1
) -> pl.DataFrame:
    """From a table of probe hits to the paths of a graph, compute which probe hits
    are the best or compatible with the best.

    To determine the best hits of each probe, the probe hits are stacked in the
    probe coordinates, and overlapping hits are split, defining local subregions.
    The, the best hit (the one with the highest similarity to the probe) is selected
    for each sub region.

    All probe hits are reevaluated and kept if its identity is higher than the identity
    threshold set by the best hit of the subregion they belong to. The identity threshold
    is a function of the divergence (1 - idt) of the best hit, which makes it more relaxed
    for lower identity values and stricter for higher identities. It is computed as:

        margin = best_hit_idt * sqrt(divergence) * div_mod
        idt_threshold = best_hit_idt - margin

    Hits with lower identity than min_idt are always filtered.

    Return simplified table with only the best (or compatible) hits.
    """

    def is_main_hit(
        hit: Interval, probe_tree: IntervalTree, div_mod: float = 0.98
    ) -> bool:
        """Given a hit and a probe_tree, determine if the hit identity is within the accepted
        identity threshold of the best homologous hit in the probe_tree.
        """

        ovlps = probe_tree.overlap(hit.begin, hit.end)
        if len(ovlps) == 0:
            return False

        homolog_hit = max(ovlps, key=lambda inter: hit.overlap_size(inter))
        idt_threshold = homolog_hit.data * (1 - (1 - homolog_hit.data) ** 2 * div_mod)

        return hit.data >= idt_threshold

    # Compute the hits coverage of the probes and which are the best hits in
    # each subregion
    TREE_COLS = ["theader", "tstart", "tend", "nident", "mismatch"]
    probe_hits = df.sort("theader").select(TREE_COLS)

    probe_trees = build_probe_trees(df, min_idt)

    # Filter the hits that are not compatible (low-identity) with the best hits on each probe
    mask = [
        (idt := Fraction(n, n + m)) >= min_idt
        and is_main_hit(Interval(start - 1, end, idt), probe_trees[probe], div_mod)
        for probe, start, end, n, m in probe_hits.iter_rows()
    ]
    clean_hits = df.sort(by="theader").filter(mask)

    # This is a multiprobe set, so pick the best probe version:
    if (clean_hits["theader"] != clean_hits["tprobe"]).any():
        clean_hits = filt_probe_hits(clean_hits, probe_trees)

    return clean_hits.sort(by="theader")


def colored_paths_extension(
    graph: Assembly_graph,
    hits: pl.DataFrame,
    max_extensions: int = 10,
    max_intron_size: int = 2000,
    extension_len_limits: ExtLimits = ExtLimits(3000, 7000, 3.0),
):
    """
    Given a graph colored by probe hits, and a set of best probe hits, proceed to
    extend all the paths that keep at least one probe hit.

    Extensions are done per probe and component, using probe coverage as key function.
    For each combination, locally generated extensions are filtered first to remove
    duplicates and contained paths. Then, only the best max_extensions are kept,
    which are selected using kmedioids on the jaccard distance of the edges of the paths
    (to select a more diverse set of extensions).

    Return a dictionary with the names of the original paths as keys and the list of
    possible extensions as values.
    """

    def _get_max_len(tlen: int, extension_limits: ExtLimits) -> int:
        "Compute the maximum extension length allowed given the extension limits."
        floor_len, ceil_len, plen_scaling = extension_limits
        return min(max(floor_len, tlen * 3 * plen_scaling), ceil_len)

    def _hash_path(path) -> int:
        return hash(tuple(sorted(edge.id for edge in path)))

    def is_contained(ext: list[OrientedEdge], ext_sets: tuple[set[str]]) -> bool:
        """Given an extension (path), determine if it is contained within another one,
        by comparing with a collection of sets representing the alternative paths."""

        qext = set(map(attrgetter("id"), ext))
        return any(qext.issubset(edge_set) for edge_set in ext_sets if qext != edge_set)

    get_max_len = partial(_get_max_len, extension_limits=extension_len_limits)
    comp_pcov = partial(probe_cov, graph=graph)
    get_path_length = partial(get_seq_atts, graph=graph)

    # Drop path duplicates (e.g. due to multiple colors) and take the longest tlen for the exploration
    hits_iter = (
        hits.sort(["query", "tlen"], descending=True)
        .group_by("query")
        .agg(
            pl.first("tlen"),
            pl.first("comp_id"),
            pl.first("tprobe"),
            pl.first("theader"),
        )
        .select(["query", "tlen", "comp_id", "tprobe", "theader"])
        .sort(["comp_id", "tprobe", "query"])
        .iter_rows()
    )

    # path_extensions = {}
    path_extensions = defaultdict(dict)
    for (comp, probe), queries in groupby(hits_iter, key=itemgetter(2, 3)):
        # Compute extensions for all paths
        comp_ext = {}
        for name, tlen, _, _, probe_ver in queries:
            try:
                comp_ext[name] = extend_path(
                    graph.paths[name],
                    graph,
                    key=comp_pcov,
                    max_len=get_max_len(tlen),
                    max_ext=max_extensions,
                    max_intron_size=max_intron_size,
                    probe=probe_ver,
                )
            except RecursionError as err:
                err.add_note(
                    f"While extending: path={name}, probe={probe_ver} ({tlen=})"
                )
                raise

        # Remove duplicated paths
        signs = set()
        filt_paths = defaultdict(list)
        for name, paths in comp_ext.items():
            for path in paths:
                if (sign := _hash_path(path)) not in signs:
                    signs.add(sign)
                    filt_paths[name].append(path)

        # Some extensions would be contained by others, we need to get rid of those
        ext_sets = tuple(
            set(map(attrgetter("id"), ext))
            for ext in chain.from_iterable(filt_paths.values())
        )

        filt_comp_ext = defaultdict(list)
        for name, extensions in comp_ext.items():
            for ext in extensions:
                if not is_contained(ext, ext_sets):
                    filt_comp_ext[name].append(ext)

        # Now use Kmedioids to select K extensions
        selected_ext = select_k_paths(
            list(chain.from_iterable(filt_comp_ext.values())),
            max_extensions,
            graph,
            key=lambda path: (comp_pcov(path), get_path_length(path)),
        )

        # Finally recover the original path names for logging
        final_comp_ext = defaultdict(list)
        path2name = {
            _hash_path(path): name
            for name, paths in filt_comp_ext.items()
            for path in paths
        }
        for ext in selected_ext:
            final_comp_ext[path2name[_hash_path(ext)]].append(ext)

        path_extensions[probe].update(final_comp_ext)

    return path_extensions


def summarize_extensions(
    ext_dict: dict[str, list[Path_on_graph]],
    path2comp: dict[str, int],
    graph: Assembly_graph,
) -> pl.DataFrame:
    """Summarize the final extended paths into a dataframe for logging.
    Information includes the original and extended names and paths, the component
    in the graph and the probe used to do the extension."""

    def log_path(path: Path_on_graph) -> str:
        return ",".join(edge.id for edge in path.edges)

    def extract_info(path: str, ext: Path_on_graph) -> tuple[str, str, int, str, str]:
        probe = re.search(PROBE_PAT, ext.name)[1]
        orig_path = log_path(graph.paths[path])
        ext_path = log_path(ext)
        return path, ext.name, path2comp[path], probe, orig_path, ext_path

    log_iter = transpose(
        extract_info(path, ext) for path, exts in ext_dict.items() for ext in exts
    )
    log_df = pl.DataFrame(log_iter, orient="col", schema=LOG_COLS)

    return (
        log_df.with_columns(
            id=pl.col("orig_name").str.extract(NODE_ID_PAT).cast(pl.Int64)
        )
        .sort("id")
        .drop("id")
        .insert_column(3, (pl.col("orig_path") != pl.col("ext_path")).alias("extended"))
    )


def snakemake_call(snakemake):
    with open(snakemake.log[0], "w") as outlog:
        sys.stdout = sys.stderr = outlog
        sys.setrecursionlimit(5000)

        graph_path = Path(snakemake.input.graph)
        table_path = Path(snakemake.input.table)
        out = Path(snakemake.output[0])
        floor_len = (
            snakemake.params.floor_len
        )  # maximum extension of a path in nucleotide
        ceil_len = (
            snakemake.params.ceil_len
        )  # maximum extension of a path in nucleotide
        plen_scaling = (
            snakemake.params.plen_scaling
        )  # float that multiplies the length of the probe for an alternative extension
        # The selected extension is the max of these two thresholds.
        max_intron_size = snakemake.params.max_intron_size
        pat = snakemake.params.probe_pattern

        extension_limits = ExtLimits(floor_len, ceil_len, plen_scaling)
        max_extensions = snakemake.params.max_extensions

        # Load the matches
        df = pl.read_csv(
            table_path, separator="\t", has_header=False, new_columns=DMND_COLS
        )

        # Load the graph and a dictionary edge -> components
        graph = Assembly_graph(graph_path)
        comp_ids = {
            edge: comp_id
            for comp_id, edges in graph.components.items()
            for edge in edges
        }

        # Add component information to the table of matches
        path2comp = {
            path: comp_ids[graph.paths[path].edges[0].id]
            for path in df["query"].unique()
        }

        pre_comp = df.with_columns(
            cis=pl.col("qend") > pl.col("qstart"),
            comp_id=pl.col("query").replace_strict(path2comp),
            tprobe=pl.col("theader").str.extract(pat),
        )

        # Find which probe is best for each component and color the graph
        final_hits = find_best_probe_hits(pre_comp)
        graph.color_edges(final_hits)

        # If there are no connections in the graph, there is nothing to extend.
        newpaths = defaultdict(list)
        if not hasattr(graph, "K"):
            for name, path in graph.paths.items():
                if not (
                    colors := Counter(
                        probe
                        for edge in path.edges
                        for probe in graph.edge_dict[edge.id].get_colors()
                    )
                ):
                    continue

                probe = re.search(pat, colors.most_common(1)[0][0])[1]
                path.name = f'{out.stem}-{probe}_{path.name.replace("NODE", "EDGE")}'
                newpaths[name].append(path)

        # Explore the colored graph:
        else:
            path_extensions = colored_paths_extension(
                graph, final_hits, max_extensions, max_intron_size, extension_limits
            )

            # Make the new paths out of the extensions
            node_pat = re.compile(NODE_ID_PAT)
            counter = count(1)
            for probe, path_ext in path_extensions.items():
                for path, extensions in sorted(
                    path_ext.items(), key=lambda x: int(node_pat.match(x[0])[1])
                ):
                    for extension in extensions:
                        name = f"{out.stem}-{probe}_"
                        name += make_path_name(extension, next(counter), graph)
                        newpaths[path].append(Path_on_graph(name, extension))

        # Log the extensions: old_path -> new_path
        # for path, exts in newpaths.items():
        #     for newpath in exts:
        #         print(f"{path}\t{newpath.name}")

        log_df = summarize_extensions(newpaths, path2comp, graph)
        log_df.write_csv(sys.stdout, separator="\t")

        # Finally write the new sequences
        seqs_iter = (
            graph.retrieve_path(path) for path in chain.from_iterable(newpaths.values())
        )
        SeqIO.write(seqs_iter, out, "fasta")


def main():
    import sys

    if len(sys.argv) != 5:
        print(
            "Usage: python gfa_graph.py <graph.gfa> <hits.tsv> <output.fasta> <extensions.log>"
        )
        sys.exit(1)

    class Run:
        def __init__(self, **kwargs):
            for key, val in kwargs.items():
                setattr(self, key, val)

    # Mock the snakemake object
    snakemake = Run(
        input=Run(graph=sys.argv[1], table=sys.argv[2]),
        output=[sys.argv[3]],
        log=[sys.argv[4]],
        params=Run(
            floor_len=3000,
            ceil_len=7000,
            plen_scaling=3,
            max_extensions=5,
            max_intron_size=2000,
            probe_pattern=r"(\d+)$",
        ),
    )

    snakemake_call(snakemake)


if __name__ == "__main__":
    if "snakemake" in globals():
        snakemake_call(snakemake)
    else:
        main()
