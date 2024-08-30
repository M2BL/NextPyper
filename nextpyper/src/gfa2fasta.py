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
from typing import NamedTuple
from functools import reduce, partial
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


# add docstring
class Segment(NamedTuple):
    id: str
    seq: Seq


# add docstring
class SeqPath(NamedTuple):
    name: str
    path: tuple[str, ...]


def parse_pline(p_line: str) -> SeqPath:
    _, name, path, *_ = p_line.split()
    return SeqPath(name, tuple(path.split(",")))


def parse_sline(s_line: str) -> Segment:
    _, name, seq, *_ = s_line.split()
    return Segment(name, Seq(seq))


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


def parse_graph_lines(graph_lines: list[str]) -> tuple[list[str], list[str], int]:
    s_lines = []
    p_lines = []
    K = None
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


def paths_to_recs(graph_lines: list[str]) -> list[SeqRecord]:
    s_lines, p_lines, K = parse_graph_lines(graph_lines)

    segments = {seg.id: seg for line in s_lines if (seg := parse_sline(line))}
    paths = [parse_pline(line) for line in p_lines]

    return [
        SeqRecord(
            seq=get_path_sequence(path, segments, K),
            id=path.name,
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
