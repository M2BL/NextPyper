#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2024
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.
"""
Functions and classes for parsing the gfa of assembly graphs produced by SPAdes with the --custom-hmms flag.
HMM profiles for probes have been used during the assembly, so that the 'hmm_statistics.txt' file from SPAdes
contains the graph edges that have a hmm match.
#  Usage example:
    gfa_file = ".../test_data/test_clustering/assembly_graph_after_simplification.gfa"
    hmm_statistics = "../test_data/test_clustering/hmm_statistics.txt"
    # parse the SPAdes output, perform the computation:
    components = filter_components_hmm(gfa_file, hmm_statistics)
    The output consist of a list of Component objects, that correspond to the assembly subgraph that have an HMM match.
"""

__version__ = "0.1"


# =======================================================================================
#               IMPORTS
# =======================================================================================
from collections import defaultdict
from dataclasses import dataclass, field
from operator import itemgetter, attrgetter
from itertools import chain, groupby
from pathlib import Path
from typing import Self
from union_find import UnionFind
from gfa2fasta import paths_to_recs
from Bio import SeqIO


@dataclass
class BGC_candidate:
    """
    Data structure for the parsing of the 'hmm_statistics.txt' file.
    """

    name: str
    hmms: list[str] = field(default_factory=list)
    coordinates: list[tuple[int, int]] = field(default_factory=list)
    paths: list[str] = field(default_factory=list)
    edges: list[str] = field(default_factory=list, init=False)
    lengths: list[int] = field(default_factory=list, init=False)
    max_domain_length: int = field(init=False)
    dominant_hmm: str = field(init=False)

    def __post_init__(self):
        self._get_edges()
        self._get_lengths()
        self._find_dominant_hmm()

    def _get_edges(self) -> Self:
        self.edges = list(
            set(
                chain(
                    *[
                        x.replace("+", "").replace("-", "").split(",")
                        for x in self.paths
                    ]
                )
            )
        )
        return self

    def _get_lengths(self) -> Self:
        for cor in self.coordinates:
            self.lengths.append(cor[1] - cor[0])
        return self

    def _find_dominant_hmm(self) -> Self:
        hmm_length = defaultdict(int)
        for idx, length in enumerate(self.lengths):
            hmm_length[self.hmms[idx]] += length
        dominant_hmm = sorted(hmm_length.items(), key=lambda x: x[1], reverse=True)[0]
        self.max_domain_length = dominant_hmm[1]
        self.dominant_hmm = dominant_hmm[0]
        return self

    def get_max_length(self) -> int:
        return self.max_domain_length

    def get_paths(self) -> list[str]:
        return self.paths

    def get_edges(self) -> list[str]:
        return self.edges

    def get_subgraph(self) -> str:
        return self.name

    def get_dominant_hmm(self) -> str:
        return self.dominant_hmm

    def __repr__(self):
        return (
            f"BGC_candidate(name={self.name}, hmms={self.hmms}, edges={self.edges} "
            f"lengths={self.lengths} paths={self.paths}, dominant hmm={self.dominant_hmm})"
        )


@dataclass
class Component:
    """
    Data structure for the components from the assembly graph that have a match with the probes HMMs.
    Attributes
    ----------
    -edges: the name of the edges that make up the component.
    -subgraph: the name of the subgraphs from the HMM search.
    -hmm: the name of the probe that have the most HMM matches in the component.
    -paths: the list of paths matched by HMM profiles. Each path is represented as a string.
    """

    hmm: str
    edges: set[str] = field(default_factory=set)
    subgraphs: list[str] = field(default_factory=list)
    paths: list[list[str]] = field(default_factory=list)

    def get_edges(self):
        return self.edges

    def get_subgraphs(self):
        return self.subgraphs

    def get_hmm(self):
        return self.hmm

    def get_paths(self):
        return self.paths


def components_from_gfa(gfa_file: str) -> list[set[str]]:
    """
    Split nodes in a gfa file per component and for each component return a set of node ids.
    :param gfa_file:
    :return:
    """
    gfa_file = Path(gfa_file)
    all_nodes = []
    edges: list[tuple] = []
    for line in open(gfa_file):
        if line.startswith("L"):
            splt_line = line.split("\t")
            edges.append((splt_line[1], splt_line[3]))
        elif line.startswith("S"):
            split_line = line.split("\t")
            all_nodes.append(split_line[1])
    with UnionFind() as UV:
        for edge in edges:
            UV.union(edge)
    components = UV.get_components()
    components_list = [list(component) for component in components]
    single_nodes = set(all_nodes).difference(set(list(chain(*components_list))))
    components.extend([{node} for node in single_nodes])
    return components


def matched_edges_from_hmm(hmm_stat_file: str, min_domain_len=20) -> dict[str, str]:
    """
    From a 'hmm_statistics.txt' file generated from spades with the --custom-hmms flag, retrieve the edges that
    have matched a hmm profile.
    :param hmm_stat_file:
    :param min_domain_len: minimum number of matched aa by the hmm profile over the whole length of a contig.
    :return: a dict with edge id as key and hmm id as value.
    """
    with Path(hmm_stat_file).open() as file:
        lines = file.readlines()

    def filter_subgraphs(list_BGS: list[BGC_candidate]) -> list[BGC_candidate]:
        """
        Remove duplicated subgraphs from a list of BGS candidates.
        """
        target = list_BGS[0]
        paths = target.get_paths()
        filtered_subgraphs = [target]
        names = [target.get_subgraph()]
        queue = list_BGS[1:]
        while queue:
            target = queue[0]
            path = target.get_paths()
            for p in path:
                if p not in paths:
                    paths.append(p)
                    if (name := target.get_subgraph()) not in names:
                        names.append(name)
                        filtered_subgraphs.append(target)
            queue = queue[1:]
        return filtered_subgraphs

    bgc_candidates = []
    subgraphs = None
    hmms = []
    coordinates = []
    paths = []
    domain_flag = False
    edge_flag = False
    for idx, line in enumerate(lines):
        if line.startswith("BGC subgraph"):
            base_name = line.split()[2]
            if subgraphs:
                bgc_candidates.extend(filter_subgraphs(subgraphs))
            subgraphs = []
            continue
        if line.startswith("BGC candidate"):
            name = f"{base_name}_candidate_{line.split()[2]}"
            hmms = lines[idx + 1].strip().split("-")
            continue

        if line.startswith("Domain coordinates:"):
            domain_flag = True
            continue
        if line.startswith("Edge order:"):
            domain_flag = False
            edge_flag = True
            continue
        if line.startswith("Path"):
            if hmms:
                subgraphs.append(BGC_candidate(name, hmms, coordinates, paths))
            name = ""
            hmms = []
            coordinates = []
            paths = []
            edge_flag = False
            continue
        if domain_flag:
            splt_domain = line.strip().split()
            assert len(splt_domain) == 2, f"[Error] line {line} has the wrong format"
            coordinates.append(tuple(map(int, splt_domain)))
        if edge_flag:
            paths.append(line.strip().replace(";", ""))
    # process the last subgraph
    if subgraphs:
        bgc_candidates.extend(filter_subgraphs(subgraphs))

    edges_dict: dict[str, list[BGC_candidate]] = defaultdict(list)
    for bgc in bgc_candidates:
        for edge in bgc.get_edges():
            edges_dict[edge].append(bgc)

    return edges_dict


def filter_components_hmm(
    gfa_file, hmm_stat_file, min_domain_len=20
) -> list[Component]:
    """
    Filter components from a gfa file based on a hmm profile file.
    :param gfa_file:
    :param hmm_stat_file:
    :param min_domain_len: minimum matching length of the HMM profile on the scaffold in aa.
    :return: a list of components with path that match at least one probe over 'min_domain_len- amino acids.
    """
    all_components = []
    components = components_from_gfa(gfa_file)
    matched_edges = matched_edges_from_hmm(hmm_stat_file, min_domain_len)

    for component in components:
        bgc_matches = []
        for edge in component:
            if (bgc_list := matched_edges.get(edge)) is not None:
                #  filter by probe matching length
                bgc_matches.extend(
                    [bgc for bgc in bgc_list if bgc.get_max_length() > min_domain_len]
                )

        if bgc_matches:
            hmm_matches = [bgc.get_dominant_hmm() for bgc in bgc_matches]
            dominant_hmm = max(set(hmm_matches), key=hmm_matches.count)
            subgraphs = list(set([bgc.get_subgraph() for bgc in bgc_matches]))
            paths = list(set(chain(*[bgc.get_paths() for bgc in bgc_matches])))
            all_components.append(Component(dominant_hmm, component, subgraphs, paths))

    return all_components


def split_into_hmms(
    gfa_path: Path,
    components: list[Component],
    outdir: Path,
    prefix: str = "",
    write_graphs: bool = True,
    write_seqs: bool = False,
) -> None:
    """Given an assembly graph and a list of components, look for the components in the graph
    group them by hmm and write them together in the output directory specified. The structure
    of the output directory is: <outdir>/<hmm>.gfa.

    :param gfa_path: path to assembly graph to split.
    :param components: list of components to find in the graph.
    :param outdir: path where to write the split components.
    :param prefix: prefix to add to the path/sequence names.
    :param write_graphs: Whether to write a gfa for each component found.
    :param write_seqs: Whether to write a fasta per component with its corresponding sequences.
    :return:
    """

    # Helper functions
    get_hmm = attrgetter("hmm")
    linked_nodes = itemgetter(1, 3)

    # Group components by hmm and name them.
    comp_dict = {
        f"{hmm}_c{i}": comp
        for hmm, group in groupby(sorted(components, key=get_hmm), key=get_hmm)
        for i, comp in enumerate(group)
    }

    # Make the reverse dictionary edge to named component.
    node2comp = {edge: name for name, comp in comp_dict.items() for edge in comp.edges}

    # Make a dictionary that will hold the lines to write for each component
    comp_lines = defaultdict(list)
    nodes = [None]
    K = None

    with gfa_path.open() as file:
        for line in file:
            match line[0]:
                case "H":
                    header = line
                case "S":
                    nodes = line.split()[1:2]
                case "L":
                    if not K:
                        K = int(line.split()[5].rstrip("M"))
                    nodes = list(linked_nodes(line.split()))
                case "P":
                    nodes = [node.rstrip("-+") for node in line.split()[2].split(",")]
                case "J":
                    continue
                case _:
                    raise NotImplementedError(f"ERROR: found line of type {line[0]}")

            main_comp = list(
                set(comp for node in nodes if (comp := node2comp.get(node)))
            )

            # Ensure that main component is a single one. More than one should be impossible.
            if len(main_comp) > 1:
                raise ValueError(
                    f"ERROR: {line=} is related to multiple components: {main_comp}, when it should be a single one."
                )

            # Append that line to the appropiate component.
            elif len(main_comp) == 1:
                comp_lines[main_comp[0]].append(line)

            # If zero components, just continue to the next line (implicitly coded).

        # Now Populate the different components with the corresponding Paths
        for name, comp in comp_dict.items():
            for i, path in enumerate(comp.paths, 1):
                # All the edges in the path have to be in the component. They might have been filtered out.
                if any(p[:-1] not in comp.edges for p in path.split(",")):
                    continue

                path_name = f"{prefix}{name}_p{i}"
                p_line = "\t".join(["P", path_name, path, "*"]) + "\n"
                comp_lines[name].append(p_line)

        # Now write the multiple files for the found components.
        get_hmm = lambda x: x[0][: x[0].rfind("_")]
        score_line = lambda line: 0 if line[0] == "S" else 1 if line[0] == "L" else 2
        outdir.mkdir(exist_ok=True, parents=True)

        for hmm, group in groupby(sorted(comp_lines.items(), key=get_hmm), key=get_hmm):
            _, lines = zip(*group)
            subgraph = sorted(chain.from_iterable(lines), key=score_line)

            # Write graphs (.gfa)
            if write_graphs:
                with (outdir / f"{hmm}.gfa").open("w") as file:
                    file.write(header)
                    file.write("".join(subgraph))

            # Write sequences (.fasta)
            if write_seqs:
                file = outdir / f"{hmm}.fasta"
                SeqIO.write(paths_to_recs(subgraph, suffix_KC=True, K=K), file, "fasta")


def main(): ...


if __name__ == "__main__":
    main()
