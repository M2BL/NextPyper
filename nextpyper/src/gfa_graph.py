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
from itertools import chain, count, repeat, starmap
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
import pandas as pd
from intervaltree import IntervalTree, Interval

from graph_alns_parser import Read
from union_find import UnionFind


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


class ColorInfo(NamedTuple):
    probe: str
    tlen: int
    cis: bool


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
            end1, start2 = alns[0].coordinates[:, -1]
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

        edges = self.paths[path].edges
        edge, *edges = edges
        length = len(self.edge_dict[edge.id])
        start, end = 0, length
        intervals = [Interval(start, end, edge.id)]

        while edges:
            edge, *edges = edges
            start += length - self.K
            length = len(self.edge_dict[edge.id])
            end += length - self.K
            intervals.append(Interval(start, end, edge.id))

        return intervals

    def color_edges(self, hits: pl.DataFrame) -> Self:
        """Use a table of probe hits to color the edges of the graph, populating the attribute
        'matching exons' of each edge with the interval they match of the best probe of
        the graph's component they belong.
        """

        # We group by component to start to operate in the graph
        comp_df = hits.group_by(["comp_id"]).agg(
            pl.first("theader"),
            pl.first("tlen"),
            pl.col("query"),
            pl.col("cis"),
            pl.col("qstart"),
            pl.col("qend"),
            pl.col("tstart"),
            pl.col("tend"),
        )

        for (
            comp_id,
            probe,
            tlen,
            paths,
            ciss,
            mqstarts,
            mqends,
            mtstarts,
            mtends,
        ) in comp_df.iter_rows():

            for path, direc, tstarts, tends, qstarts, qends in zip(
                paths, ciss, mtstarts, mtends, mqstarts, mqends
            ):
                if direc:
                    qstarts, qends = qends, qstarts

                colorinfo = ColorInfo(probe, tlen, direc)
                path_intervals = self._get_path_intervals(path)
                qinter = sorted(map(Interval, qstarts, qends))
                mod_tstarts = (tstart - 1 for tstart in tstarts)
                tinter = sorted(map(Interval, mod_tstarts, tends, repeat(colorinfo)))

                for qint, tint in zip(qinter, tinter):
                    for pint in path_intervals:
                        if qint.overlaps(pint):
                            edge_id = pint.data
                            self.edge_dict[edge_id].matching_exons.append(tint)

        for edge in self.edge_dict.values():
            if edge.matching_exons:
                tree = IntervalTree(edge.matching_exons)
                tree.merge_overlaps(data_reducer=lambda x, y: x)
                edge.matching_exons = sorted(tree)


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
    max_len: int = 5000,
    max_extensions: int = 10,
    goal: Optional[OrientedEdge] = None,
    key: Optional[Callable[[list[OrientedEdge]], int]] = None,
):
    """

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
        return sum(len(graph.edge_dict[edge[0]]) for edge in path)

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
        if get_path_len(current_path, graph) > max_len:
            extensions.add(current_path[:])
            return

        # Reached a dead-end
        if not (neighbors := graph.graph[edge]):
            extensions.add(current_path[:])
            return

        # Explore neighbors
        for neighbor in neighbors:
            dfs_helper(neighbor, current_path[:], extensions)

    if key is None:
        key = partial(get_path_len, graph=graph)

    current_path = []
    extensions = OptimalExtension(size=max_extensions, key=key)

    dfs_helper(start, current_path, extensions)

    return extensions


def probe_cov(path: list[OrientedEdge], graph: Assembly_graph) -> float:
    tree = IntervalTree(
        chain.from_iterable(graph.edge_dict[edge.id].matching_exons for edge in path)
    )
    if tree.is_empty():
        return 0.0

    tree.merge_overlaps(data_reducer=lambda x, _: x)
    cov_bases = sum(interval.length() for interval in tree)
    tlen = next(iter(tree)).data.tlen

    return cov_bases / tlen


def extend_path(
    path: Path_on_graph,
    graph: Assembly_graph,
    max_len: int = 5000,
    max_ext: int = 10,
    key: Optional[Callable[[list[OrientedEdge]], int]] = None,
) -> list[OrientedEdge]:
    """Given an assembly graph and a path, extend the given path following the graph topology
    using a Depth First Search.

    Return a list with all the extended paths, represented as list of oriented edges.
    """

    return [
        path.edges[:-1] + list(starmap(OrientedEdge, ext))
        for ext in dfs_track_paths(
            graph, path.edges[-1], max_len=max_len, max_extensions=max_ext, key=key
        )
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


def effective_cov(starts: list[int], ends: list[int]) -> int:
    """
    Compute the effective coverage of a set of matches (intervals), by
    merging the intervals they define with their starts and ends.

    Return the total amount of bases covered [int].
    """

    gen_iter = zip(chain.from_iterable(starts), chain.from_iterable(ends))
    tree = IntervalTree.from_tuples(gen_iter)
    tree.merge_overlaps()
    return sum(inter.length() for inter in tree) + 1


def find_components_best_probe(df: pl.DataFrame) -> pl.DataFrame:
    """From a table of probe hits to the paths of a graph, compute which probe hits
    best each of the components on the graph and add this information to the table.

    Return simplified table with only the hits of the best probe per component.
    """

    # Collect all the hits of the same probe on the same path
    df2 = df.group_by(["query", "theader", "cis"]).agg(
        pl.sum("nident"),
        pl.sum("mismatch"),
        pl.sum("gapopen"),
        pl.first("tlen"),
        pl.first("comp_id"),
        pl.col("qstart"),
        pl.col("qend"),
        pl.col("tstart"),
        pl.col("tend"),
    )

    # Collect all the hits of a probe in a given component
    df3 = df2.group_by(["theader", "comp_id"]).agg(
        pl.sum("nident"),
        pl.sum("mismatch"),
        pl.first("tlen"),
        pl.concat_list("tstart"),
        pl.concat_list("tend"),
    )

    # Compute the "effective coverage" on the probe by hits on the same components
    df4 = df3.with_columns(
        eff_cov=pl.Series(
            values=starmap(effective_cov, df3[["tstart", "tend"]].iter_rows()),
            dtype=pl.Int64,
        )
    )

    # Select the nest probe for each component (the dominant probe), by "alternative coverage"
    # criterion. The alternative coverage is the Non-redundant number of matches found in the probe,
    # over the length of the probe.
    best_probes = (
        df4.with_columns(
            alt_cov=pl.col("eff_cov")
            * (pl.col("nident") / (pl.col("nident") + pl.col("mismatch")))
            / pl.col("tlen")
        )
        .sort(by=["comp_id", "alt_cov"], descending=True)
        .group_by(["comp_id"])
        .agg(pl.first("theader"))
    )

    return (
        df2.join(best_probes, on="comp_id", suffix="_best")
        .filter(pl.col("theader") == pl.col("theader_best"))
        .sort(by="comp_id")
    )


def snakemake_call2(snakemake):
    with open(snakemake.log[0], "w") as outlog:
        sys.stdout = sys.stderr = outlog

        graph_path = Path(snakemake.input.graph)
        table_path = Path(snakemake.input.table)
        out = Path(snakemake.output[0])
        floor_len = (
            snakemake.params.floor_len
        )  # maximum extension of a path in nucleotide
        plen_scaling = (
            snakemake.params.plen_scaling
        )  # float that multiplies the length of the probe for an alternative extension
        # The selected extension is the max of these two thresholds.

        max_extensions = snakemake.params.max_extensions

        # Load the matches
        df = pl.read_csv(table_path, separator="\t", has_header=True)

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
        )

        # Find which probe is best for each component and color the graph
        final_hits = find_components_best_probe(pre_comp)
        graph.color_edges(final_hits)

        # Explore the colored graph:
        # ToDo: Implement


def snakemake_call(snakemake):
    with open(snakemake.log[0], "w") as outlog:
        sys.stdout = sys.stderr = outlog

        graph_path = Path(snakemake.input.graph)
        table_path = Path(snakemake.input.table)
        out = Path(snakemake.output[0])
        floor_len = (
            snakemake.params.floor_len
        )  # maximum extension of a path in nucleotide
        plen_scaling = (
            snakemake.params.plen_scaling
        )  # float that multiplies the length of the probe for an alternative extension
        # The selected extension is the max of these two thresholds.

        max_extensions = snakemake.params.max_extensions

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
                max_ext=max_extensions,
            )
            for name, plen in match_paths_plen.items()
        }

        newpaths = {
            path: [
                Path_on_graph(
                    make_path_name(extension, next(counter), graph), extension
                )
                for extension in extensions
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
