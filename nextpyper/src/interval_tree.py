#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2021
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.
#
"""
Classes used to create an interval search tree.
"""
__version__ = "0.1"

import itertools

# =======================================================================================
#               IMPORTS
# =======================================================================================
from dataclasses import dataclass, field
from itertools import combinations
from math import floor
from typing import Any, Union, Literal

# =======================================================================================
#               CLASSES
# =======================================================================================

Number = Union[float, int]


@dataclass
class Overlap:
    start: Number
    end: Number
    first_segment: Any
    second_segment: Any

    def get_length(self) -> Number:
        return self.end - self.start

    def get_start(self) -> Number:
        return self.start

    def get_end(self) -> Number:
        return self.end


class Node:
    """
    Data structure use in the interval search tree. It relies on the node data
    structure found at:
    https://algs4.cs.princeton.edu/93intersection/IntervalST.java.html
    Attributes:
        -interval: Interval
        -value: Any, value associated with a node.
        -left_node: Node, node on the left branch
        -right_node: Node, node on the right branch
    """

    def __init__(
        self,
        interval: "Interval",
        value: Any,
        left_node: "Node" = None,
        right_node: "Node" = None,
    ):
        self.interval = interval
        self.value = [value]
        self.left = left_node
        self.right = right_node
        self.N = 1  # size of subtree
        self.max = interval.getMax  # max endpoint in subtree rooted at this note

    @classmethod
    def from_tuple(cls, tuple_lo_hi: tuple[Number], value: Any) -> "Node":
        """
        Create a node from a tuple (lower bound of the interval, higher bound).
        """
        if len(tuple_lo_hi) != 2:
            raise ValueError(f"Node: tuple can contain only two values {tuple_lo_hi}")
        node = cls(Interval(tuple_lo_hi[0], tuple_lo_hi[1]), value)
        return node

    def __str__(self) -> str:
        return f"Node with interval: {self.interval} and value {self.value}"

    def __repr__(self) -> str:
        return f"Node(interval={self.interval}, value={self.value})"


class Interval:

    """
    Data structure use in the interval search tree. It relies on the interval data
    structure found at:
    https://algs4.cs.princeton.edu/93intersection/IntervalST.java.html
    Attributes:
        -lo: int, lower limit
        -hi: int, upper limit
    """

    def __init__(self, lo: Number, hi: Number):
        if lo is None or hi is None:
            raise ValueError(f"Interval: None endpoint not allowed: ({lo}, {hi})")

        if lo < 0 or hi <= 0:
            raise ValueError(
                f"Interval: Negative lower endpoint and null upper endpoint are not allowed: ({lo}, {hi})"
            )

        if lo >= hi:
            raise ValueError(
                f"Interval: Lower endpoint must be higher than upper endpoint: ({lo}, {hi})"
            )
        self.lo = lo
        self.hi = hi

    @classmethod
    def from_tuple(cls, tuple_lo_hi: tuple[Number]) -> "Interval":
        """
        Create interval from tuple (lower bound, higher bound).
        """
        if len(tuple_lo_hi) != 2:
            raise ValueError(f"Node: tuple can contain only two values {tuple_lo_hi} ")
        lo = tuple_lo_hi[0]
        hi = tuple_lo_hi[1]
        interval = cls(lo, hi)
        return interval

    def compare(self, other: "Interval") -> Literal[-1, 0, 1]:
        """
        If self and other interval are the same return 0.
        If self interval starts before other return -1, otherwise return 1.
        If self and other have the same start, break the tie using the upper limit.
        """
        if floor(self.lo) == floor(other.lo) and floor(self.hi) == floor(other.hi):
            return 0
        if floor(self.lo) < floor(other.lo):
            return -1
        if floor(self.lo) > floor(other.lo):
            return 1
        if self.hi < other.hi:
            return -1
        #  self.hi > other.hi
        return 1

    def contains(self, other: "Interval") -> bool:
        """
        True if other interval is included in the self interval.
        False otherwise.
        """
        if self.lo <= other.lo and other.hi <= self.hi:
            return True
        return False

    def intersects(self, other: "Interval") -> bool:
        """
        True if the self and other interval intersect.
        False otherwise.
        """
        if self.hi <= other.lo:
            return False
        if other.hi <= self.lo:
            return False
        return True

    def overlap(self, other: "Interval") -> float:
        starts = self.lo, other.lo
        ends = self.hi, other.hi

        if (ovlp := min(ends) - max(starts)) < 0:
            ovlp = 0

        return ovlp / (max(ends) - min(starts))

    def length(self) -> int:
        """
        Return the length of the interval.
        """
        return self.hi - floor(self.lo)

    def getMax(self) -> int:
        """
        Return the upper bound of the interval.
        """
        return self.hi

    def getMin(self) -> int:
        """
        Return the lower bound of the interval.
        """
        return self.lo

    def __repr__(self) -> str:
        return f"({self.lo}, {self.hi})"

    def __getitem__(self, item) -> Number:
        if item == 0:
            return self.lo
        if item == 1:
            return self.hi
        raise IndexError("list index out of range")


class IntervalST:
    """
    Create interval search tree using a modified version of the IntervalST data object:
        https://algs4.cs.princeton.edu/93intersection/IntervalST.java.html
    """

    def __init__(self, root=None):
        self.root = root

    def put(self, interval: Interval, value: Any) -> None:
        """
        New data object in inserted at the root of the tree (no random insertion).
        If interval is alredy present in the tree, the new value is appended to
        the value list of the interval.
        """
        self.root = self._rootInsert(self.root, interval, value)

    def _rootInsert(self, node: Node, interval: Interval, value: Any) -> Node:
        """
        Insert a new node.
        """
        if node is None:
            return Node(interval, value)
        cmp = interval.compare(node.interval)
        if cmp < 0:  #  adding left
            node.left = self._rootInsert(node.left, interval, value)
            node = self._rotR(node)
        if cmp > 0:  #  adding right
            node.right = self._rootInsert(node.right, interval, value)
            node = self._rotL(node)
        else:
            if value not in node.value:
                node.value.append(value)
        return node

    def get(self, interval: Interval) -> list[Any]:
        """Return value associated with the given key
        if no such value, return null"""
        return self._get(self.root, interval)

    def _get(self, node: Node, interval: Interval) -> list[Any]:
        """
        Find recursively the value of a node.
        """
        if node is None:
            return None
        cmp = interval.compare(node.interval)
        if cmp < 0:  #  searching left
            node.left = self._get(node.left, interval)
        elif cmp > 0:  #  searching right
            node.right = self._get(node.right, interval)
        else:
            return node.value

    def _rotR(self, node: Node) -> Node:
        """
        Rotate tree right.
        """
        x = node.left
        node.left = x.right
        x.right = node
        self._fix(node)
        self._fix(x)
        return x

    def _rotL(self, node: Node) -> Node:
        """
        Rotate tree left.
        """
        x = node.right
        node.right = x.left
        x.left = node
        self._fix(node)
        self._fix(x)
        return x

    def _fix(self, node: Node) -> None:
        """
        Fix auxilliar information (subtree count and max fields).
        """
        if node is None:
            return
        node.N = 1 + self._size(node.left) + self._size(node.right)
        node.max = max(
            node.interval.hi, self._nodeMax(node.left), self._nodeMax(node.right)
        )

    def _nodeMax(self, node: Node) -> int:
        """
        Get the maximum of a node.
        """
        if node is None:
            return 0
        return node.max

    def size(self) -> int:
        """ "
        Get the size of a tree.
        """
        return self._size(self.root)

    def searchContains(self, interval: Interval) -> list[Node]:
        """
        Main algorithm for searching the interval tree. In comparison to the classical
        algorithm it has been modified to work with included interval instead of mere
        overlapping intervals.
        """
        L = []
        self._searchContains(self.root, interval, L)
        return L

    def _searchContains(
        self, node: Node, interval: Interval, L: list[Node]
    ) -> list[Node]:
        """
        Search for intersection with the query interval,
        the current node and if not found search recursively
        the tree from the left and from the right.
        """
        found1 = False
        found2 = False
        found3 = False
        if node is None:
            return False
        if node.interval.contains(interval):
            L.append(node)
            found1 = True
        if (
            node.left is not None
            and node.left.max >= interval.hi
            and node.left.max >= interval.lo
        ):
            found2 = self._searchContains(node.left, interval, L)
        if (node.right is not None) and (
            found2
            or (node.left is None)
            or node.left.max <= interval.lo
            or (node.right.max >= interval.hi and node.right.max >= interval.lo)
        ):
            found3 = self._searchContains(node.right, interval, L)

        return found1 or found2 or found3

    def searchIntersect(self, interval: Interval) -> list[Node]:
        """
        Algorithm for searching for intersecting intervals. Returns a list
        of tuples of interval bounderies and a list of values of the intersecting nodes.
        """
        L = []
        self._searchIntersect(self.root, interval, L)
        return L

    def _searchIntersect(
        self, node: Node, interval: Interval, L: list[Node]
    ) -> list[Node]:
        while node is not None:
            if interval.intersects(node.interval):
                L.append(
                    (
                        node.interval,
                        node.value,
                    )
                )
            if node.left is None:
                node = node.right
            elif node.left.max < interval.lo:
                node = node.right
            else:
                node = node.left

    def get_all_intersections(self):
        result = []
        nodes = [self.root]
        while nodes:
            node = nodes[0]
            nodes = nodes[1:]
            if node is None:
                continue
            interval = node.interval
            start = interval.getMin()
            end = interval.getMax()
            values = node.value
            if len(values) > 1:
                value_comb = list(combinations(values, 2))
                # print(f"{value_comb =}")
                result.extend([Overlap(start, end, *value) for value in value_comb])

            if node.left is None:
                target_node = node.right
            elif node.left.max < interval.lo:
                target_node = node.right
            else:
                target_node = node.left
            matching_intervals = []
            self._searchIntersect(target_node, interval, matching_intervals)
            # print("found matches", matching_intervals)
            for match in matching_intervals:
                match_interval = match[0]
                match_start = match_interval.getMin()
                match_end = match_interval.getMax()
                match_values = match[1]
                left = max(start, match_start)
                right = min(end, match_end)
                if left == right:
                    continue
                for value in values:
                    for match_value in match_values:
                        result.append(Overlap(left, right, value, match_value))
            nodes.extend([node.left, node.right])
        return sorted(result, key=lambda x: x.start)

    def tree_traversal_bsf(self) -> list[tuple[Interval]]:
        """
        Use Breadth-first search to retrieve all nodes from an interval tree.
        For each node we get a tuple.

        Returns
        -------
        list[tuple]
            Each tuple contains the interval as an Interval object and the list
            of values associated with the node.

        """
        result = []
        nodes = [self.root]
        while nodes:
            node = nodes[0]
            nodes = nodes[1:]
            if node is None:
                continue
            interval = node.interval
            values = node.value
            result.extend(
                [
                    (
                        interval,
                        v,
                    )
                    for v in values
                ]
            )
            nodes.extend([node.left, node.right])
        return sorted(result, key=lambda x: x[1])

    def height(self) -> int:
        """
        Height of the Search Tree
        """
        return self._height(self.root)

    def _height(self, node) -> int:
        if node is None:
            return 0
        return 1 + max(self._height(node.left), self._height(node.right))

    def _size(self, node) -> int:
        if node is None:
            return 0
        return node.N


if __name__ == "__main__":
    print("starting")
    n1 = Node.from_tuple((1, 2), 2)
    n1.value.extend([3])
    print(n1)
    itr = Interval.from_tuple((1, 2))
    print(f"Interval is: {itr}")
    print(f"upper limit of interval is: {itr[1]}")
    # print(Interval(0, 10).compare(Interval(2, 10)))
    # print(Interval(0.2, 10).compare(Interval(0.99, 10)))
    # s = Interval(0, 10); print(s)
    # print(Interval(1676.162799346671, 1698).contains(Interval(1690, 1698)))
    # print(Interval(1676.540526197659, 1711).contains(Interval(1690, 1699)))
    ST = IntervalST()
    ST.put(Interval(5, 10), "a")
    ST.put(Interval(6, 11), "b")
    ST.put(Interval(6, 11), "c")
    ST.put(Interval(20, 25), "d")
    ST.put(Interval(19, 21), "e")
    # print(ST.searchIntersect(Interval(6, 7)))
    # print(ST.searchIntersect(Interval(4, 7)))
    print("all intersections", ST.get_all_intersections())
    # print(ST.tree_traversal_bsf())
    # print("search:", [x.value for x in ST.searchContains(Interval(7, 10))])

    # print("result value is ", ST.get(Interval(1, 6)))
    # print("search:", [x.value for x in ST.search(Interval(20, 21))])
