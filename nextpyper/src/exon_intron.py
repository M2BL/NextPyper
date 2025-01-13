#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2025
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.

__version__ = "0.1"

# =======================================================================================
#               IMPORTS
# =======================================================================================
from collections import defaultdict, namedtuple
from dataclasses import dataclass, field, fields
import re
from typing import Final, Optional, Self, TypedDict, Literal, Any

# Named tuple that keeps the nodes in a path either as their suffixes or their cds and the sum of their length.
GraphPath = namedtuple("GraphPath", ["path", "length"])

# =======================================================================================
#               FUNCTIONS
# =======================================================================================

# =======================================================================================
#               CLASSES
# =======================================================================================


@dataclass
class Exon:
    __slots__ = ["start", "end"]
    start: int
    end: int

    def overlaps(self, other: "Exon") -> bool:
        first, second = sorted([self, other], key=lambda x: x.start)
        if first.end > second.start:
            return True
        return False


@dataclass
class Node:
    key: int | str
    interval: Exon = field(repr=True)
    children: list["Node"] = field(init=False, default_factory=list)
    root: bool = field(init=False, default=True)

    def get_key(self):
        return self.key

    def overlap(self, other: "Node") -> bool:
        if self.interval.overlaps(other.interval):
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
        return self.interval.end - self.interval.start

    def get_children(self) -> list["Node"]:
        return self.children

    def get_interval(self) -> "Interval":
        return self.interval


@dataclass
class ItervalGraph:
    """
    Data structure that create a graph of non overlapping intervals (nodes).
    Intervals are ordered by starting position. Starting with each target interval, we iterate over the rest
    of the intervals that start at greater position. If these interval don't overlap with the target
    we draw an edge, otherwise we skip the interval. Once the graph is constructed, we get the starting nodes
    as those that do not have preceding intervals (i.e. no parent node). Finding the longest path of non-overlapping
    intervals requires a DSF from each starting node.
    Attributes
    ----------
    -interval_list: a list of Cds objects.
    -nodes: a list of Node objects.
    -node_dict: a dictionary of Node objects as value and cds suffix number as key.
    -root_nodes: starting node that do not have parents.
    -all_graph_path: all possible paths from a root_node to a leaf node.
    """

    interval_list: list["Interval"]
    nodes: list[Node] = field(init=False, default_factory=list)
    node_dict: dict[int, Node] = field(init=False, default_factory=dict)
    root_nodes: list[int] = field(init=False, default_factory=list)
    all_graph_path: list[GraphPath] = field(init=False, default_factory=list)

    def __post_init__(self):
        self._create_nodes()
        self._create_edges()
        self._find_roots()
        self._explore_graph()

    def _create_nodes(self) -> Self:
        nodes = []
        for idx, interval in enumerate(self.interval_list):
            node = Node(idx, interval)
            nodes.append(node)
        self.nodes = sorted(nodes, key=lambda n: n.interval.start)
        return self

    def _create_edges(self) -> Self:
        """Connect the node to all other nodes whose cds is left of the current node cds and whose sequences
        do not overlap"""
        self.node_dict = {node.get_key(): node for node in self.nodes}
        for idx, node in enumerate(self.nodes):
            for other_node in self.nodes[idx + 1 :]:
                if node.overlap(other_node):
                    continue
                node.add_child(other_node)
                other_node.not_root()
        return self

    def _find_roots(self) -> Self:
        """
        A root note has no parent node. Non-root nodes were tagged during graph construction.
        """
        self.root_nodes = [node.get_key() for node in self.nodes if node.is_root()]
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
                    [node.get_key() for node in current_path_copy],
                    current_length_copy,
                )
                self.all_graph_path.append(new_path)
                return
            else:
                for child in current_node.get_children():
                    dfs_util(child, current_path_copy, current_length_copy)

        for key in self.root_nodes:
            root = self.node_dict[key]
            dfs_util(root)
        return self

    def get_best_path(self) -> Optional[GraphPath]:
        """Out of all the path found with the DFS, select the longest"""
        if not self.all_graph_path:
            return
        paths = sorted(self.all_graph_path, key=lambda x: x.length, reverse=True)
        best = paths[0]
        interval_path = [self.node_dict[node].get_interval() for node in best.path]
        return GraphPath(interval_path, best.length)


if __name__ == "__main__":
    inputs = [
        (143, 178),
        (107, 142),
        (50, 70),
        (60, 80),
        (71, 110),
        (111, 142),
        (40, 70),
        (150, 170),
        (190, 191),
    ]
    interval_list = [Exon(*x) for x in inputs]
    print(interval_list)
    IG = ItervalGraph(interval_list)
    print(IG.get_best_path())
