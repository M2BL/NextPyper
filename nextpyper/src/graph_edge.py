#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2024
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.
#
"""
Classes used to parse a blunted compacted graph, import mappings and explore the paths on the graph.

"""

__version__ = "0.1"

# =======================================================================================
#               IMPORTS
# =======================================================================================
from collections import defaultdict
import copy
from dataclasses import dataclass, field, make_dataclass
from itertools import chain
from pathlib import Path
import sys
from typing import Final, Optional, ClassVar, TypedDict, NewType, Literal, Self

from Bio.Seq import Seq

import hmm_search
from hmm_search import Hmmer_result, Profile_hits, Domain, Path_nodes, Node_hit

# =======================================================================================
#               CUSTOM EXCEPTIONS
# =======================================================================================


# =======================================================================================
#               CLASSES
# =======================================================================================


@dataclass
class Sequence:
    """
    A dataclass that stores sequence name and Seq.
    """

    sequence_id: str
    sequence: Seq


@dataclass
class Edge(Sequence):
    """
    A dataclass that stores all the information about an edge in a graph including
    mapping from probes (hmm, minigraph, SPAligner...).
    Attributes
    ----------
        -edge_id: Name of edge.
        -orientation: Orientation of the edge.
        -in_edges: incoming edges list.
        -out_edges: outgoing edges list.
        -hmm_probe_hit: dict of name of probe as key and Domain namedtuple as value.
    """

    edge_id: str
    orientation: Literal["+", "-"]
    in_edges: list["Edge"] = field(default_factory=list)
    out_edges: list["Edge"] = field(default_factory=list)
    hmm_probe_hit: dict[str, Domain] = field(default_factory=dict)

    def __post_init__(self):
        self._correct_orientation()

    def _correct_orientation(self) -> Self:
        """
        Reverse complement sequences with '-' orientation.
        :return:
        """
        if self.orientation == "-":
            self.sequence = self.sequence.reverse_complement()
        return self

    def __repr__(self):
        if self.hmm_probe_hit.values():
            return (
                f"Edge: {self.edge_id}, {self.sequence[:20]}, {self.orientation},"
                f" in_edges={[x.edge_id for x in self.in_edges]}, "
                f"out_edges={[x.edge_id for x in self.out_edges]},"
                f" {list(self.hmm_probe_hit.items())[0]}"
            )
        return (
            f"Edge: {self.edge_id}, {self.sequence[:20]},{self.orientation},"
            f" in_edges={[x.edge_id for x in self.in_edges]}, "
            f"out_edges={[x.edge_id for x in self.out_edges]},"
            f" No probe hit"
        )


@dataclass
class Live_path:
    """
    Dataclass that represents a path that is followed on a graph and is defined
    by hmm hits on a probe.
    Attributes
    ----------
        -probe: name of probe.
        -all_edges: list of all edge names that make up the path.
        -edges_with_mapping: list of all edge names that actually have a hmm hit.
        -last_mapped_edge: name of edge that was last mapped on path.
            Used to revert the all_edges attibute to a previous state.
        -last_domain: nametuple that represents the alignment of the last edge with hmm hit.
        -probe_distance: Counter that keep track of the distance between the current edge and
            the last edge with a hmm hit. It is reset everytime a new hmm hit if found.
    """

    probe: str
    all_edges: list[str]
    edges_with_mapping: list[str]
    last_mapped_edge: str
    last_domain: Domain
    probe_distance: int = field(default=0, init=False)

    def add_all_edges(self, new_edge: str) -> Self:
        """
        Adds new edge to path.
        :param new_edge: edge name
        :return:
        """
        self.all_edges.append(new_edge)
        return self

    def get_all_edges(self) -> Self:
        return tuple(self.all_edges)

    def add_edge_with_mapping(self, new_edge: str) -> Self:
        self.edges_with_mapping.append(new_edge)
        return self

    def set_mapped_edge(self, new_edge: str) -> Self:
        self.last_mapped_edge = new_edge
        return self

    def rollback_path(self) -> Self:
        self.all_edges = self.all_edges.split(self.last_mapped_edge)[0]
        return self

    def update_domain(self, new_domain: Domain) -> Self:
        self.last_domain = new_domain
        return self

    def get_probe_distance(self):
        return self.probe_distance

    def increment_distance(self) -> Self:
        self.probe_distance += 1
        return self

    def get_probe_end(self):
        return self.last_domain.profile_end

    def reset_distance(self) -> Self:
        self.probe_distance = 0
        return self

    def __repr__(self):
        return f"Live_path('{self.probe}'), path: {self.get_all_edges()}, mapped_path: {self.edges_with_mapping}'"


@dataclass
class Gfa:
    """
    A dataclass that stores a GFA graph, with methods to retrieve paths that have been selected
    my mapping probes.
    Attributes
    ----------
        -gfa_path: Path to the gfa file.
        -sequences: dict with node names as keys and Seq object as value.
        -graph: Intermediate object.Parsing of the 'L' line in the gfa. Keys are the left node and value are the right nodes
            e.g. '4+':['3-', '6+']
        -rev: Dict used to reverse the edge direction.
    __post_init__
        -edges: Final product of the parsing. Combines the sequences and graph object.
            Keys are the edge name (e.g. '4+') and values are Edge object.
        -external_edges: Edges without neighbour at the left.
        -terminal edges: Edges without neighbour at the right.

    """

    gfa_path: Path
    sequences: dict[str, Seq] = field(default_factory=dict)
    graph: defaultdict[list] = field(default_factory=lambda: defaultdict(list))
    paths: list[Live_path] = field(default_factory=list, init=False)
    rev = {"+": "-", "-": "+"}
    edges: dict[str, Edge] = field(default_factory=dict, init=False)
    external_edges: list[Edge] = field(default_factory=list, init=False)
    terminal_edges: list[Edge] = field(default_factory=list, init=False)

    def __post_init__(self):
        self._from_file()
        self._create_edges()

    def get_edges(self) -> dict[str, Edge]:
        return self.edges

    def _create_edges(self):
        """
        Combine the node number, with its sequence and orientation into a new Edge object.
        :return: populate self.edges.
        """
        all_edge_labels = list(self.graph.keys())
        all_edge_labels.extend(list(set(chain.from_iterable(self.graph.values()))))
        for label in all_edge_labels:
            name = label[:-1]
            orientation = label[-1]
            self.edges[label] = Edge(
                edge_id=label,
                sequence_id=name,
                sequence=self.sequences[name].sequence,
                orientation=orientation,
            )
        # Populate the in and out edges
        for label, out_edges in self.graph.items():
            target_node = self.edges[label]
            out_edges_list = [self.edges[item] for item in out_edges]
            target_node.out_edges.extend(out_edges_list)
            for node in out_edges_list:
                node.in_edges.append(target_node)

        # find external (start) and terminal (end) edges
        self._find_external_terminal()
        return self

    def _from_file(self) -> Self:
        """
        Parse the GFA file and create a dict of Sequence objects and a dict of connections between nodes.
        :return: Populate self.sequences and self.graph objects.
        """
        with open(self.gfa_path, "r") as fin:
            for line in fin.readlines():
                if line.startswith("S"):
                    lst = line.strip().split("\t")[1:]
                    node_id, seq = lst[0], lst[1]
                    self.sequences[node_id] = Sequence(node_id, Seq(seq))

                elif line.startswith("L"):
                    _, node_id1, pos1, node_id2, pos2, match = line.strip().split("\t")
                    self.graph[node_id1 + pos1].append(node_id2 + pos2)
        return self

    # Add exceptions for edges that are not found
    def load_hmm(self, hmm_node_hits: defaultdict) -> Self:
        """
        Add hmm information to each edge and orient the edges.
        If the edge has been used in the graph it has an orientation.
        If the edge has not been connected to another edge, it has no orientation, but if it has
        a hmm hit, we give it the orientation of the probe and adjust the strand of the sequence
        accordingly.
        :return: Populates the self.edges dictionary.
        """
        for node_id in hmm_node_hits:
            for probe, profile_hit in hmm_node_hits[node_id].items():
                domain = profile_hit.domain
                orientation = "+" if domain.node_start < domain.node_end else "-"
                node_name = node_id + orientation
                if orientation == "-":
                    profile_hit.invert_RC()  # Invert the node start and end coordinates
                # If edge has already an orientation, i.e. it has been used in the graph or hit by hmm
                if (edge := self.edges.get(node_name, None)) is not None:
                    ##edge.hmm_probe_hit[probe].append(profile_hit)
                    edge.hmm_probe_hit[probe] = profile_hit
                # Else if only the sequence exists
                else:
                    new_edge = Edge(
                        edge_id=node_name,
                        sequence=self.sequences[node_id].sequence,
                        orientation=orientation,
                    )
                    ##new_edge.hmm_probe_hit[probe].append(profile_hit)
                    edge.hmm_probe_hit[probe] = profile_hit
                    self.edges[node_name] = new_edge
        return self

    def _find_external_terminal(self) -> Self:
        """
        Find which edge have no left neighbour (external) and which have no right neighbour (terminal)
        :return: populates the self.external_edges and self.terminal_edges attributes.
        """
        for edge in self.edges.values():
            #  Unconnected edges
            if not edge.in_edges and not edge.out_edges:
                self.external_edges.append(edge)
                self.terminal_edges.append(edge)
                continue
            if not edge.in_edges:
                self.external_edges.append(edge)
            else:
                self.terminal_edges.append(edge)

        print("external_edges", [x.edge_id for x in self.external_edges])
        print("terminal_edges", [x.edge_id for x in self.terminal_edges])
        return self

    def hmm_DSF(self) -> Self:  # list[Live_path]:
        """
        Explore the graph starting from external edges.
        :return: Populate self.paths.
        """
        all_paths = []
        print("starting DSF")

        def dfs_util(
            active_edge: Edge,
            live_paths: dict[
                str, Live_path
            ],  # key probe name, value Live_path that is the path that are still investigated
            max_nbr_edges_not_hit=2,  # this could be changed for a maximum sequence length
        ):
            """
            Method to explore all possible paths in a graph using depth-first search.
            :param active_edge: The current edge.
            :param live_paths: dict of probe names as keys and Live_path instances as value.
                When a path is created, a key/value is added; when a path is closed, it is deleted.
            :param max_nbr_edges_not_hit: The maximum number of edges without a hmm hit
                that are allowed before closing a path
            :return:
            """
            edge_id = active_edge.edge_id
            print(f"current edge {edge_id}")
            print(f"current active paths: {list(live_paths.items())}")
            # if hmm hits on this edge
            if hmm_hits := active_edge.hmm_probe_hit:
                edge_probes = set(hmm_hits.keys())  # probes that have hmm hits
                active_path_probes = set(
                    live_paths.keys()
                )  # probes that are active on path

                #  take care of unconnected edges.
                if (
                    active_edge in self.external_edges
                    and active_edge in self.terminal_edges
                ):
                    # add path object for the unconnected edge.
                    probe = list(edge_probes)[0]
                    live_path = Live_path(
                        probe=probe,
                        all_edges=[edge_id],
                        edges_with_mapping=[edge_id],
                        last_mapped_edge=edge_id,
                        last_domain=hmm_hits[probe].domain,
                    )
                    all_paths.append(live_path)
                    print(f"unconnected edge: {edge_id}")
                    return

                # not Live_path, create one for each probe
                if not live_paths:
                    for probe in edge_probes:
                        print(hmm_hits[probe])
                        live_paths[probe] = Live_path(
                            probe=probe,
                            all_edges=[edge_id],
                            edges_with_mapping=[edge_id],
                            last_mapped_edge=edge_id,
                            last_domain=hmm_hits[probe].domain,
                        )

                # there is at least one live path
                else:
                    common_probes = (
                        edge_probes & active_path_probes
                    )  # probes both in active path and on active edge
                    probes_not_path = (
                        active_path_probes - edge_probes
                    )  # probes missing in path
                    probes_not_edge = (
                        edge_probes - active_path_probes
                    )  # probes missing in active edge

                    for probe in common_probes:
                        live_path = live_paths[probe]
                        # Check for domain compatibility
                        live_probe_end = live_path.get_probe_end()
                        edge_domain = hmm_hits[probe].domain
                        edge_probe_end = edge_domain.profile_start
                        # the two domains are incompatible
                        if live_probe_end > edge_probe_end:
                            live_path.rollback_path()
                            all_paths.append(live_path)
                            del live_paths[probe]
                            return
                            # create new path starting with the current edge?

                        # the two domains are compatible extend the path,
                        # active domain becomes active_edge domain and distance counter is reset.
                        else:
                            live_path.add_all_edges(edge_id)
                            live_path.add_edge_with_mapping(edge_id)
                            live_path.set_mapped_edge(edge_id)
                            live_path.update_domain(edge_domain)
                            live_path.reset_distance()

                    # Create new paths for new probes in active_edge
                    for probe in probes_not_path:
                        live_paths[probe] = Live_path(
                            probe=probe,
                            all_edges=[edge_id],
                            edges_with_mapping=[edge_id],
                            last_mapped_edge=edge_id,
                            last_domain=hmm_hits[probe].domain,
                        )

                    # either update path with missing domain or end path
                    for probe in probes_not_edge:
                        live_path = live_paths[probe]
                        distance = live_path.get_probe_distance()
                        if distance > max_nbr_edges_not_hit:
                            live_path.rollback_path()
                            all_paths.append(live_path)
                            del live_path
                            return
                        else:
                            live_path.increment_distance()

            # if no hmm hits on this edge
            # it depends on if there are some active paths
            # if no active path go to the next node
            # if active path exist, check if we are still in range
            if not live_paths:
                pass
            else:
                for probe in live_paths:
                    live_path = live_paths[probe]
                    distance = live_path.get_probe_distance()
                    if distance > max_nbr_edges_not_hit:
                        live_path.rollback_path()
                        all_paths.append(live_path)
                        del live_path
                        return
                    else:
                        live_path.increment_distance()

            # current edge is a terminal edge.
            if not active_edge.out_edges:
                for live_path in live_paths.values():
                    all_paths.append(live_path)
                    print(f"all_paths = {all_paths}")
                print(f"reaching a terminal edge, {all_paths=}")
                return

            # current edge is not a terminal edge
            else:
                print(
                    f"moving to the neighbours: {[x.edge_id for x in active_edge.out_edges]}"
                )
                for out_edge in active_edge.out_edges:
                    dfs_util(out_edge, copy.deepcopy(live_paths))

        print(f"{self.external_edges=}")
        for edge in self.external_edges:
            print(f"starting edge is: {edge.edge_id}")
            dfs_util(edge, {})
        self.paths = all_paths
        return self


if __name__ == "__main__":
    test_graph = "/home/yjkbertrand/Documents/projects/nextpyper/test_data/batrachium/test_data/probe_3/mapping/minigraph/H1_C8_blunted_compacted.gfa"
    gfa = Gfa(test_graph)
    import pickle

    hmm_file = "/home/yjkbertrand/Documents/projects/nextpyper/test_data/batrachium/test_data/probe_3/mapping/minigraph/hmm_H1_C8.pkl"
    hmm = pickle.load(open(hmm_file, "rb"))
    gfa.load_hmm(hmm.node_hits)
    for k, v in gfa.get_edges().items():
        print(k, v)
    print(gfa.hmm_DSF())
    # print(list(gfa.get_edges().items()))
    # print(list(gfa.graph.items()))
    # print(list(gfa.sequences.keys()))
