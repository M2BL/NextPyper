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
from collections import defaultdict
from dataclasses import dataclass, field
from itertools import chain, count, starmap
from pathlib import Path
from typing import Self, Literal, Optional, NamedTuple
from functools import partial
import re
import sys

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Align import PairwiseAligner, Alignment
import pandas as pd

from graph_alns_parser import Read


# =============================================================================
#                CLASSES
# =============================================================================
LinkSupport = dict[tuple[str], int]


class OrientedEdge(NamedTuple):
    id: str
    orientation: Literal["+", "-"]
    jump: bool = False

    @classmethod
    def from_seg(cls, seg) -> "OrientedEdge":
        if match := re.match(r",|;", seg):
            return cls(seg[1:-1], seg[-1], match[0] == ";")
        else:
            return cls(seg[:-1], seg[-1])


@dataclass(slots=True)
class Path_on_graph:
    name: str
    edges: list[OrientedEdge]
    start: Optional[int] = field(default=0)
    end: Optional[int] = field(default=None)
    length: Optional[int] = field(default=None)
    score: Optional[int] = field(default=None)

    def get_parameters(self):
        return self.name, self.edges, self.start, self.end, self.length, self.score


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

    def __len__(self) -> int:
        return len(self.seq)

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
    rev = {"+": "-", "-": "+"}

    def __post_init__(self):
        self._parse_graph()

    def _parse_graph(self) -> Self:
        graph_path = Path(self.gfa_filename)
        assert graph_path.exists(), f"Graph file {self.gfa_filename} not found"
        with open(self.gfa_filename, "r") as file:
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
                    case "J":
                        _, node_id1, pos1, node_id2, pos2, _dist, *_tags = (
                            line.strip().split("\t")
                        )
                        self.graph[OrientedEdge(node_id1, pos1)].append(
                            OrientedEdge(node_id2, pos2, True)
                        )
                        self.graph[OrientedEdge(node_id2, self.rev[pos2])].append(
                            (node_id1, self.rev[pos1], True)
                        )

                        ## Add to the Jump Links (J-Lines) the connections from L-Lines to prevent graph disconnection
                        self.graph[OrientedEdge(node_id2, pos2, True)].extend(
                            self.graph[OrientedEdge(node_id2, pos2)]
                        )
                        self.graph[OrientedEdge(node_id1, self.rev[pos1], True)].extend(
                            OrientedEdge(node_id1, self.rev[pos1])
                        )

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

    def link_edges(self, reads: list[Read]) -> Self:
        """Given a list of paired reads, compute the link support that those pairs
        exhibit. Each pair that connects a composition of edges is added as a
        supported link, increasing its count in the dictionary. Compositions of
        a single edge are not taken into account.
        """

        get_edges = lambda frags: (frag.edge.id for frag in frags)

        for read in reads:
            edges = get_edges(read.fragments + read.mate.fragments)
            if len(unique_edges := set(edges)) > 1:
                self.linked_edges[tuple(sorted(unique_edges))] += 1

        return self

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
    def _ovlp_extension(seq1: Seq, seq2: Seq, gap: int = 10) -> Seq:
        """Use pairwise alignment to detect if seq1 and seq2 overlap. If they do
        Use the alignment coordinates to merge them. Otherwise, scaffold them
        with a gap of the given size."""

        aligner = PairwiseAligner(scoring="megablast", mode="local")
        alns = aligner.align(seq1, seq2)

        # Proceed as no overlap between the sequences. Add a gap of 10 N
        if alns.score < 30:
            return seq1 + Seq("N" * gap) + seq2

        # Overlap detected. Use alignment coordinates to merge the sequences
        else:
            end1, start2 = alns[0].coordinates[0, :]
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


# =============================================================================
#                FUNCTIONS
# =============================================================================


def dfs_track_paths(
    graph: Assembly_graph,
    start: OrientedEdge,
    max_len: int = 5000,
    max_extensions: int = 10,
    goal=None,
):
    def get_path_len(path: list[OrientedEdge], graph: Assembly_graph) -> int:
        return sum(len(graph.edge_dict[node[0]]) for node in path)

    def dfs_helper(node, visited, current_path, all_dead_ends):
        visited.add(node)
        current_path.append(node)
        max_len_exceeded = False

        # If we reach the goal, add the path
        if goal is not None and node == goal:
            all_dead_ends.append(current_path[:])

        # Maximum len size check (avoids long recursions)
        if get_path_len(current_path, graph) > max_len:
            max_len_exceeded = True
            all_dead_ends.append(current_path[:])

        # Get all neighbors
        neighbors = graph.graph[node]

        # If no unvisited neighbors (dead end) and no specific goal
        if goal is None and all(n in visited for n in neighbors):
            all_dead_ends.append(current_path[:])

        # Explore neighbors
        for neighbor in neighbors:
            if not max_len_exceeded and neighbor not in visited:
                dfs_helper(neighbor, visited, current_path, all_dead_ends)

        # Backtrack: remove current node from path and visited
        current_path.pop()
        visited.remove(node)

    visited = set()
    current_path = []
    all_dead_ends = []

    dfs_helper(start, visited, current_path, all_dead_ends)
    if len(all_dead_ends) > max_extensions:
        return sorted(
            all_dead_ends, key=partial(get_path_len, graph=graph), reverse=True
        )[:max_extensions]
    else:
        return all_dead_ends


def extend_path(
    path: Path_on_graph, graph: Assembly_graph, max_len: int = 5000
) -> list[OrientedEdge]:
    """Given an assembly graph and a path, extend the given path following the graph topology
    using a Depth First Search.

    Return a list with all the extended paths, represented as list of oriented edges.
    """

    return [
        path.edges[:-1] + list(starmap(OrientedEdge, ext))
        for ext in dfs_track_paths(graph, path.edges[-1], max_len=max_len)
    ]


def get_seq_atts(
    protopath: list[OrientedEdge], graph: Assembly_graph
) -> tuple[int, float]:
    """Given a list of OrientedEdges which represent a path, infer the length and coverage of
    the sequence it encodes. Returns a tuple with both values.
    """

    ids = [edge.id for edge in protopath]
    path_length = sum(graph.edge_dict[id].get_length() for id in ids) - (
        graph.K * (len(protopath) - 1)
    )
    path_cov = sum(
        graph.edge_dict[id].coverage * graph.edge_dict[id].get_length() for id in ids
    ) / (path_length - graph.K)

    return path_length, path_cov


def make_path_name(path: list[OrientedEdge], idx: int, graph):
    atts = get_seq_atts(path, graph)
    return f"EDGE_{idx}_length_{atts[0]}_cov_{atts[1]:.3f}"


def snakemake_call(snakemake):
    with open(snakemake.log[0], "w") as outlog:
        sys.stdout = sys.stderr = outlog

        graph_path = Path(snakemake.input.graph)
        table_path = Path(snakemake.input.table)
        out = Path(snakemake.output[0])
        floor_len = snakemake.params.floor_len
        plen_scaling = snakemake.params.plen_scaling

        # Try to extend only paths that match with probes
        pat = re.compile(r"(?<=NODE_)\d+")
        df = pd.read_csv(table_path, sep="\t")
        match_paths_plen = {
            query: tlen
            for _, (query, tlen) in df.loc[:, ["query", "tlen"]]
            .drop_duplicates()
            .groupby(by="query")
            .max("tlen")
            .reset_index()
            .sort_values(
                by="query",
                key=lambda names: [int(pat.search(name)[0]) for name in names],
            )
            .iterrows()
        }

        counter = count(1)
        graph = Assembly_graph(graph_path)

        path_extensions = {
            name: extend_path(
                graph.paths[name],
                graph,
                max_len=max(floor_len, plen * 3 * plen_scaling),
            )
            for name, plen in match_paths_plen.items()
        }

        newpaths = {
            path: [
                Path_on_graph(make_path_name(ext, next(counter), graph), ext)
                for ext in extensions
            ]
            for path, extensions in path_extensions.items()
        }

        for path, exts in newpaths.items():
            for newpath in exts:
                print(f"{path}\t{newpath.name}")

        seqs_iter = (
            graph.retrieve_path(path) for path in chain.from_iterable(newpaths.values())
        )
        SeqIO.write(seqs_iter, out, "fasta")


def main(): ...


if __name__ == "__main__":
    if "snakemake" in globals():
        snakemake_call(snakemake)
    else:
        main()
