#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2024
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.

"""
Functions and classes used to parse the output of TransDecoder in 'pep' format.
For each probe, the longest path obtained by combining the possible cds is calculated.
Stop codons ('*') are removed from the sequences and the cds are combined.
AA sequences are saved in fasta format.

#  Usage example:
    transdecoder_parser("../data/longest_orfs.pep", "../data/combined_longest_orfs.fasta")
"""

__version__ = "0.1"

# =======================================================================================
#               IMPORTS
# =======================================================================================
from collections import defaultdict, namedtuple
from dataclasses import dataclass, field, fields
import re
from typing import Final, Optional, Self, TypedDict, Literal, Any

from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Parse the header of a TransDecoder sequence
pattern = re.compile(
    r"""^(?P<name>.*?)\s
                    type:(?P<type>\w+?)\s
                    len:(?P<length>\d+?)\s
                    gc:universal\s
                    [^ ]+:(?P<start>\d+)
                    -(?P<end>\d+)
                    \((?P<direction>[+-])\)$""",
    re.VERBOSE,
)

# Named tuple that keeps the nodes in a path either as their suffixes or their cds and the sum of their length.
GraphPath = namedtuple("GraphPath", ["path", "length"])


@dataclass
class Node:
    suffix: str
    cds: "Cds" = field(repr=False)
    children: list["Node"] = field(init=False, default_factory=list)
    root: bool = field(init=False, default=True)

    def get_suffix(self):
        return self.suffix

    def overlap(self, other: "Node") -> bool:
        if self.cds.overlap(other.cds):
            return True
        return False

    def add_child(self, child: "Node") -> Self:
        self.children.append(child)
        return self

    def not_root(self) -> Self:
        self.root = False
        return self

    def is_root(self) -> bool:
        return self.root

    def has_children(self) -> bool:
        if self.children:
            return True
        return False

    def get_length(self) -> int:
        return self.cds.get_length()

    def get_children(self) -> list["Node"]:
        return self.children

    def get_cds(self) -> "Cds":
        return self.cds


@dataclass
class ItervalGraph:
    cds_list: list["Cds"]
    nodes: list[Node] = field(init=False, default_factory=list)
    node_dict: dict[str, Node] = field(init=False, default_factory=dict)
    root_nodes: list[str] = field(init=False, default_factory=list)
    all_graph_path: list[GraphPath] = field(init=False, default_factory=list)

    def __post_init__(self):
        self._create_nodes()
        self._create_edges()
        self._find_roots()
        self._explore_graph()

    def _create_nodes(self) -> Self:
        for cds in self.cds_list:
            node = Node(cds.get_suffix(), cds)
            self.nodes.append(node)
        return self

    def _create_edges(self) -> Self:
        """Connect the node to all other nodes whose cds is left of the current node cds and whose sequences
        do not overlap"""
        self.node_dict = {node.get_suffix(): node for node in self.nodes}
        for idx, node in enumerate(self.nodes):
            for other_node in self.nodes[idx + 1 :]:
                if node.overlap(other_node):
                    continue
                node.add_child(other_node)
                other_node.not_root()
        return self

    def _find_roots(self) -> Self:
        self.root_nodes = [node.get_suffix() for node in self.nodes if node.is_root()]
        return self

    def _explore_graph(self) -> Self:
        """
        Find all paths between root nodes (nodes without parent) and leaf nodes (nodes without children).
        :return:
        """

        def dfs_util(
            current_node: Node, current_path: list[Node] = [], current_length: int = 0
        ) -> GraphPath:
            current_path_copy = current_path.copy()
            current_length_copy = current_length
            current_path_copy.append(current_node)
            current_length_copy += current_node.get_length()
            if not current_node.has_children():
                new_path = GraphPath(
                    [node.get_suffix() for node in current_path_copy],
                    current_length_copy,
                )
                self.all_graph_path.append(new_path)
                return
            else:
                for child in current_node.get_children():
                    dfs_util(child, current_path_copy, current_length_copy)

        for suffix in self.root_nodes:
            root = self.node_dict[suffix]
            dfs_util(root)
        return self

    def get_best_path(self, orientation: Literal["+", "-"]) -> Optional[GraphPath]:
        if not self.all_graph_path:
            return
        paths = sorted(self.all_graph_path, key=lambda x: x.length, reverse=True)
        best = paths[0]
        if orientation == "+":
            cds_path = [self.node_dict[node].get_cds() for node in best.path]
        else:
            new_path = best.path
            new_path.reverse()
            cds_path = [self.node_dict[node].get_cds() for node in new_path]
        return GraphPath(cds_path, best.length)


@dataclass
class Cds:
    cds_prefix: str
    cds_suffix: int
    type: str
    length: int
    start: int
    end: int
    direction: Literal["+", "-"]
    sequence: Seq = field(repr=False)

    def get_type(self):
        return self.type

    def get_suffix(self):
        return self.cds_suffix

    def get_seq(self):
        if self.type not in ["3prime_partial", "internal"]:
            return self.sequence[:-1]
        return self.sequence

    def invert_coordinates(self):
        """For reversed complemented sequences the start and end coordinates are reversed for better sorting"""
        left = self.start
        right = self.end
        self.start = right
        self.end = left

    def get_direction(self):
        return self.direction

    def get_length(self):
        return self.length

    def get_start(self):
        return self.start

    def get_end(self):
        return self.end

    def overlap(self, other: "Cds") -> bool:
        if self.get_end() > other.get_start():
            return True
        return False


def select_best_cds(cds_list: list[Cds]) -> GraphPath:
    """
    Create an interval graph for '+' and '-' list of cds. Explore all possible paths and determine the best path.
    :param cds_list:
    :return: a namedtuple containing a list of Cds objects and a total length.
    """
    partial_plus = []
    partial_minus = []
    for cds in cds_list:
        if cds.get_direction() == "+":
            partial_plus.append(cds)
        else:
            cds.invert_coordinates()
            partial_minus.append(cds)
    optimums = []

    if partial_plus:
        sorted_partial_plus = sorted(partial_plus, key=lambda cds: cds.start)
        iter_graph_plus = ItervalGraph(sorted_partial_plus)
        optimums.append(iter_graph_plus.get_best_path("+"))
    if partial_minus:
        sorted_partial_minus = sorted(partial_minus, key=lambda cds: cds.start)
        iter_graph_minus = ItervalGraph(sorted_partial_minus)
        optimums.append(iter_graph_minus.get_best_path("-"))

    best = sorted(optimums, key=lambda x: x.length, reverse=True)[0]
    return best


def td_parser(pep_path: Path) -> list[SeqRecord]:
    """
    Take a 'longest_orfs.pep' file generated by TransDecoder,
    parse it and test all combinations of possible cds in order to get the longest
    amino acid sequence. Return for each probe the longest sequence.
    :param pep_path:
    :return:
    """
    records = SeqIO.parse(pep_path, "fasta")
    cds_dict = defaultdict(list)
    for rec in records:
        cds_name_splt = rec.id.rsplit(".p", 1)
        cds_prefix = cds_name_splt[0]
        cds_suffix = cds_name_splt[1]
        match = pattern.match(rec.description)
        cds = Cds(
            cds_prefix,
            cds_suffix,
            match.group("type"),
            int(match.group("length")),
            int(match.group("start")),
            int(match.group("end")),
            match.group("direction"),
            rec.seq,
        )
        cds_dict[cds_prefix].append(cds)

    new_records = []
    for probe, cdss in cds_dict.items():
        best = select_best_cds(cdss)
        if not best:
            continue
        best_length = best.length
        best_path = best.path
        best_suffix = "-".join([cds.get_suffix() for cds in best_path])
        sequence = Seq("".join([str(cds.get_seq()) for cds in best_path]))
        direction = best_path[0].get_direction()
        description = (
            f"cds_suffixes={best_suffix}, length={best_length}, direction={direction}"
        )
        new_record = SeqRecord(id=probe, description=description, name="", seq=sequence)
        new_records.append(new_record)
    return new_records


def transdecoder_parser(path_to_pep: str, path_to_output: str) -> None:
    """
    Takes a 'longest_orfs.pep' file generated by TransDecoder as input, parse it and test all combinations of possible
    cds in order to get the longest amino acid sequence.
    Save it to fasta
    :param path_to_pep:
    :param path_to_output:
    :return:
    """
    pep_path = Path(path_to_pep)
    fasta_path = Path(path_to_output)
    new_records = td_parser(pep_path)
    SeqIO.write(new_records, fasta_path, "fasta")


if __name__ == "__main__":
    pep = "/home/yjkbertrand/Documents/projects/nextpiper/test_data_old/batrachium/targets.fasta.transdecoder_dir/longest_orfs.pep"
    # pep = "/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_transdecoder/target_1225.pep"
    new_records = td_parser(Path(pep))
    import os

    os.chdir("/home/yjkbertrand/Documents/projects/nextpiper/temp")
    SeqIO.write(new_records, "test.fasta", "fasta")
