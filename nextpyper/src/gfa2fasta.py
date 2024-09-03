#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2024
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.

"""Simple functions to get the Path sequences of an assembly graph (.gfa) in fasta format.
The functions are tested to work for standard De Brujn graphs, and assume a constant overlap
between segments. Heterogeneous overlaps are not supported. Blunted graphs (e.g. variation 
graphs) should work too, but have not been tested.

#Usage example:
    gfa_file = Path("path/to/my_graph.gfa")
    out_fasta = Path("my_path_sequences.fasta")
    paths2fasta(gfa_file, out_fasta)
    
    The headers of the sequences in the fasta files correspond to the names in the P lines.
"""

__version__ = "0.1"

# =======================================================================================
#               IMPORTS
# =======================================================================================

from pathlib import Path
from typing import Any, NamedTuple
from functools import reduce, partial
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


# add docstring
class Segment(NamedTuple):
    id: str
    seq: Seq
    tags: dict[str, Any]


# add docstring
class SeqPath(NamedTuple):
    name: str
    path: tuple[str, ...]


def parse_pline(p_line: str) -> SeqPath:
    _, name, path, *_ = p_line.split()
    return SeqPath(name, tuple(path.split(",")))


def parse_tags(tag: str) -> tuple[str, Any]:
    id, type, value = tag.split(":")
    match type:
        case "f":
            value = float(value)
        case "i":
            value = int(value)
        case "Z" | "H":
            ...
        case _:
            raise ValueError(f"Found tag of {type=}")

    return id, value


def parse_sline(s_line: str) -> Segment:
    _, name, seq, *tags = s_line.split()
    return Segment(name, Seq(seq), dict(map(parse_tags, tags)))


# add docstring
def get_seq(edge: str, segments: dict[str, Segment]) -> Seq:
    seg = segments[edge[:-1]]
    return seg.seq.reverse_complement() if edge.endswith("-") else seg.seq


# add docstring
def link_edges(seq1: Seq, seq2: Seq, K: int) -> Seq:
    return seq1 + seq2[K:]


# add docstring for the attributes
def get_path_sequence(
    path: SeqPath, segments: dict[str, Segment], K: int | None = None
) -> Seq:
    """Given a Path encoded in a p_line return the Sequence given by that path
    taken from the corresponding segments (encoded in the s_lines).

    If K is given, the graph is assumed to have overlaps of the given size
    between its segments, and hence it will be accounted for when reconstructing the
    sequence.
    """

    plink_edges = partial(link_edges, K=K)

    return (
        get_seq(path.path[0], segments)
        if len(path.path) == 1
        else reduce(plink_edges, (get_seq(piece, segments) for piece in path.path))
    )


def parse_graph_lines(
    graph_lines: list[str], K: int | None = None
) -> tuple[list[str], list[str], int]:
    s_lines = []
    p_lines = []

    for line in graph_lines:
        match line[0]:
            case "S":
                s_lines.append(line)
            case "L":
                if not K:
                    K = int(line.split()[5].rstrip("M"))
            case "P":
                p_lines.append(line)
            case _:
                continue

    return s_lines, p_lines, K


def paths_to_recs(
    graph_lines: list[str], suffix_KC=False, K: int | None = None
) -> list[SeqRecord]:
    """Given a list with gfa lines, parse and generate Sequences encoded in its
    paths (P lines) and return them as SeqRecord objects.

    If suffix_KC is true, compute and include in the sequence name the Kmer count
    for the encoded sequence according to the KC values of the edges that compose it.
    """

    def path_depth(path_segs: list[Segment], K: int) -> float:

        tot_kc = sum(seg.tags.get("KC", 0) for seg in path_segs)
        len_path = sum(len(seg.seq) for seg in path_segs)

        return tot_kc / (len_path - (len(path_segs)) * K)

    s_lines, p_lines, newK = parse_graph_lines(graph_lines)
    if K is None:
        K = newK

    segments = {seg.id: seg for line in s_lines if (seg := parse_sline(line))}
    paths = [parse_pline(line) for line in p_lines]

    if suffix_KC:
        get_depth = (
            lambda path: f"_DP{path_depth([segments[p[:-1]] for p in path.path], K):.2f}"
        )

    return [
        SeqRecord(
            seq=get_path_sequence(path, segments, K),
            id=f"{path.name}{get_depth(path) if suffix_KC else ''}",
            name="",
            description="",
        )
        for path in paths
    ]


def paths_to_fasta(gfa: Path, out_fasta: Path) -> None:
    """Given a gfa file, write sequences encoded in its P lines
    to the given outfile path in fasta format."""

    assert gfa.exists(), f"Graph file {gfa} not found"
    with gfa.open() as file:
        SeqIO.write(paths_to_recs(file), out_fasta, "fasta")
