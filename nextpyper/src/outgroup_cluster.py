#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2024
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.
"""
Functions and classes for clustering contigs according to their pairwise similarity.
We first align the probe's AA sequence on each contig and record the coordinates of
the alignments. Contigs are then aligned pairwise along their overlapping probe regions.
Pairwise distances are used to cluster the contigs using HDBscan (or union-find for small samples).
Individual clusters can be saved in fasta format.
#  Usage example:
    probe_fasta = ".../test_data/test_clustering/probe_3_aa.fasta"
    contig_fasta = "../test_data/test_clustering/gene_3_all.fasta"
    # load the data, perform the computation:
    hdb = HDBcluster(probe_fasta, contig_fasta)
    # save the results:
    hdb.save_clusters("../test_data/test_clustering/clusters")
    The resulting fasta files are named using the suffix of the probe file name
    and the index of the cluster as suffix, e.g. probe_3_aa_45.fasta.


"""

__version__ = "0.1"

# =======================================================================================
#               IMPORTS
# =======================================================================================

from collections import deque
from dataclasses import dataclass, field
from importlib import reload
from io import StringIO
from itertools import groupby, chain
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import textwrap
import time
from typing import Optional, Self, Literal, TypeAlias

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import numpy.typing as npt
from sklearn.cluster import HDBSCAN

from alignment import Alignment_fragment
import gff_parser

gff_parser = reload(gff_parser)
from gff_parser import Fragment, Cds
from interval_tree import IntervalST, Interval
from miniprot import MiniprotInit, run_miniprot
from union_find import UnionFind


# =======================================================================================
#               EXCEPTIOMS
# =======================================================================================



# =======================================================================================
#               FUNCTIONS
# =======================================================================================




@dataclass
class Outgrp_cluster(MiniprotInit):
    """
    Data structure for the fragment object produced by the parsing of the miniprot output.
    Each fragment is one exon. The sequences are clustered with HDBscan. Clustered sequences are processed. Low similarity
    to the probe and short sequences are sieved out. The resulting sequences can be saved as fasta.
    by simi
    Attributes
    ----------
    -probe_fasta: path to the amino acid probe fasta file.
    -contigs_fasta: path to the nucleotide contigs fasta file.
    -treads: for miniprot.
    -min_probe_contig_sim: mapping similarity for miniprot.
    -min_fragment_cov: fraction of the probe covered by a contig for miniprot.
    -min_contig_overlap: minimum overlap between two contigs to compute a distance.
    -max_clustering_dist: maximum distance between two contigs to associate them with the union-find algorithm.
    -min_cluster_size: minimum size of clusters for the HDBSCAN algorithm.
    -min_contig_length: minimum length of a contig after trimming.
    -keep_singletons: if True, keep unclustered contigs.

    Post Init
    -distance_matrix: matrix of pairwise distances indexed on the cds_dict.
    -clusters: nested list of contig names that have been clustered.
        Set to the empty list, if no sequence match the probe.
    """

    outgroups: list[str] = field(default_factory=list)
    probe_coverage:float = field(default=0.5)
    min_contig_overlap: int = field(default=300) #to compute distance
    keep_singletons: bool = field(default=True)
    cds_dict: dict[str, Cds] = field(init=False,repr=False, default_factory=dict)

    distance_matrix: npt.ArrayLike = field(init=False, repr=False, default=None)

    def __post_init__(self):
        super().__post_init__()
        assert 0 < self.probe_coverage<= 1, "[ERROR] probe_coverage must be between 0 and 1"
        assert self.outgroups, "[ERROR] outgroups must be defined"
        outgroup_scaffolds = []
        for scaffold in self.contigs_dict:
            for out in self.outgroups:
                if out in scaffold:
                    outgroup_scaffolds.append(scaffold)
        assert outgroup_scaffolds, "[ERROR] no scaffold from the outgroup was found"
        print("Running miniprot")
        with tempfile.TemporaryDirectory() as tmpdirname:
            for (contig, record) in self.contigs_dict.items():
                print(f"working on contig {contig}")
                contig_path = Path(tmpdirname) / f"{contig}.fas"
                cds = self._run_miniprot(record, self.probes_path, contig_path)
                print(cds)
                if not cds.is_empty():
                    self.cds_dict[contig] = cds

        print("building distance matrix")
        tic = time.perf_counter()
        self._create_distance_matrix()
        toc = time.perf_counter()
        print(f"Computed distances in  {toc - tic:0.4f} seconds")
        print(self.distance_matrix)
        contig_names = list(self.contigs_dict.keys())
        for x in range(len(contig_names)):
            i = list(self.cds_dict.keys()).index(contig_names[x])
            for y in range(x + 1, len(self.contigs_dict)):
                j = list(self.cds_dict.keys()).index(contig_names[y])
                distance = self.distance_matrix[i, j]
                print(f"distance between {contig_names[x]} and {contig_names[y]} is {distance}")


    def _create_distance_matrix(
        self,
    ) -> Self:
        """
        Compute the distances between all fragments from different contigs. If the fragments overlap
        a distance is calculated, otherwise it is set to None.
        :return: A distance matrix
        """
        items = list(self.cds_dict.items())
        lists = deque()
        for idx, item in enumerate(items):
            query_name_node, query_cds = item
            print(f"{query_name_node=}")
            distances = [0]
            for target in items[idx + 1 :]:
                target_name_node, target_cds = target
                print(f"{target_name_node=}")
                print(f"{query_cds.get_global_sim()}, {target_cds.get_global_sim()}")
                if any([query_cds.is_empty(), target_cds.is_empty()]):
                    distance = math.nan
                elif any(
                    [
                        query_cds.get_global_sim() < self.min_probe_contig_sim,
                        target_cds.get_global_sim() < self.min_probe_contig_sim,
                    ]
                ):
                    distance = math.nan
                else:
                    aln_frag = Alignment_fragment(
                        query_cds,
                        target_cds,
                        self.probes_dict,
                        self.min_contig_overlap,
                    )
                    if (distance := aln_frag.get_similarity()) is None:
                        distance = math.nan


                distances.append(distance)
            lists.append(distances)
        pad = len(max(lists, key=len))
        half_matrix = np.array([[math.nan] * (pad - len(i)) + i for i in lists])
        matrix = np.triu(half_matrix) + np.triu(half_matrix, 1).T
        self.distance_matrix = np.nan_to_num(matrix, nan=100)
        return self

    def _cluster_unionfind(self, contig_names: list[str]) -> list[str]:
        """
        Cluster contigs with the UnionFind algorithm. Contigs need to overlap and their pairwise distance below
        the max_clustering_dist value.
        :return: list of group labels that are the names of the contigs.
        """
        with UnionFind() as UV:
            # names = list(contig_names.keys())
            for x in range(len(contig_names)):
                i = list(self.cds_dict.keys()).index(contig_names[x])
                for y in range(x + 1, len(contig_names)):
                    j = list(self.cds_dict.keys()).index(contig_names[y])
                    distance = self.distance_matrix[i, j]
                    if distance < self.max_clustering_dist:
                        UV.union((x, y))
        components = UV.get_components()
        final_groups = []
        for component in components:
            final_groups.append(sorted([contig_names[c] for c in component]))
        if self.keep_singletons:
            singletons = set(contig_names) - set(list(chain(*final_groups)))
            final_groups.extend([singleton] for singleton in singletons)
        return sorted(final_groups)

    def _cluster_hdbscan(self, contig_names: list[str]) -> list[str]:
        """
        Using the distance matrix, cluster sequences with HDBSCAN.
        'https://scikit-learn.org/stable/modules/generated/sklearn.cluster.HDBSCAN.html'
        :return: list of group labels that are the names of the contigs.
        """
        clusters = HDBSCAN(
            min_cluster_size=self.min_cluster_size,
            min_samples=1,
            metric="precomputed",
            allow_single_cluster=True,
        ).fit(self.distance_matrix)
        labels = clusters.labels_
        record_clusters = zip(contig_names, labels)
        record_clusters_sorted = sorted(record_clusters, key=lambda x: x[1])
        groups = groupby(record_clusters_sorted, lambda x: x[1])
        record_groups = []
        for key, group in groups:
            record_groups.append(list(group))
        final_groups = []
        for group in record_groups:
            if group[0][1] == -1:
                continue
            final_groups.append(sorted([x[0] for x in group]))
        if self.keep_singletons:
            singletons = set(contig_names) - set(list(chain(*final_groups)))
            final_groups.extend([singleton] for singleton in singletons)
        return sorted(final_groups)

    def _use_UF(self) -> bool:
        """
        Determine if the UnionFind algorithm should be used.
        It is used in case the dataset contains too few sequences to use HDBscan.
        :return:
        """

        if len(self.contigs_dict) < self.min_cluster_size:
            return True
        upper_triangle = np.triu(self.distance_matrix)
        i = 0
        for row in upper_triangle:
            i += 1
            if np.count_nonzero(row[i:] < 100) > self.min_cluster_size:
                return False
        else:
            return True

    def _get_clusters(self) -> Self:
        """
        Run the UnionFind algorithm for dataset with less overlapping contigs than the min_cluster_size value.
        Otherwise, use the HDBSCAN algorithm for clustering and only then run UnionFind to separate connected components.
        :return: Populate the 'cluster' attribute
        """
        # print(self.distance_matrix)
        # Case no Cds match the probe
        if all((cds.has_empty_cds() for cds in self.cds_dict.values())):
            return []
        # Case matching below thresholds
        if set(self.distance_matrix.flatten()) == {0.0, 100.0}:
            global_sims = sorted(
                [cds.get_global_sim() for cds in self.cds_dict.values()]
            )
            fill = textwrap.fill(
                "[Warning] There are matches, but they are not reported as they register as paralogs. "
                f"There are {len(global_sims)} sequences that match the probe with similarities between {global_sims[0]} and {global_sims[1]}",
                width=90,
                subsequent_indent=" " * 11,
            )
            print(fill)
            return []
        if self._use_UF():
            print("running UF")
            self.clusters = self._cluster_unionfind(list(self.contigs_dict.keys()))
            return self

        print("running HDBSCAN")
        clusters = self._cluster_hdbscan(self.contigs_dict)
        final_clusters = []
        print("running UF after HDBSCAN")
        for cluster in clusters:
            if len(cluster) == 1:
                final_clusters.append(cluster)
                continue
            # Double check that each sequence in a given cluster found by HDBscan
            # does overlap with at least another sequence from the same cluster.
            # This check is necessary in case of clusters made out of out-layers.
            new_clusters = self._cluster_unionfind(cluster)

            final_clusters.extend(new_clusters)
        self.clusters = final_clusters
        return self

    def get_clusters(self) -> list[list[str]]:
        """
        Return the nested list of cluster labels (i.e. contig names)
        :return:
        """
        return self.clusters

    def save_distance_matrix(self, path: Path) -> None:
        """
        Save the distance matrix to disk.
        :return:
        """
        if self.distance_matrix is not None:
            print(f"saving distance matrix as {path}")
            np.save(Path(path), self.distance_matrix)

    def save_clusters(self, fasta_folder: str) -> None:
        """
        Given a list of labels generated by HDBSCAN, save the clusters to a fasta file in the specified folder.
        :param labels: list from DBSCAN with index position that corresponds to the contigs in the keys of
        the 'clean_fragments' dictionary.
        :param fasta_folder: folder where the fasta files are saved.
        :return:
        """
        if not self.clusters:
            print("No clusters found.")
            return
        print(f"Saving clusters to fasta folder {fasta_folder}")
        Path(fasta_folder).mkdir(parents=True, exist_ok=True)
        fasta_prefix = self.probe_fasta.stem
        idx = 0
        for gp in self.clusters:
            trimmed_records = self._trim_msa(gp)
            if trimmed_records:
                name = f"{fasta_prefix}_{idx}.fasta"
                name_fasta = Path(fasta_folder) / name
                SeqIO.write(trimmed_records, name_fasta, "fasta")
                print(f"saving fasta {name_fasta}")
                idx += 1

def main():

    probe_fasta = "/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_paralogy_2/best_probe_6488.fasta"
    contig_fasta = "/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_paralogy_2/saute/paralogy_6488_vsearch_con.fasta"
    matrix = "/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_paralogy_2/blosum62.csv"
    #run_miniprot_boundary(probe_fasta, contig_fasta, matrix)
    OC = Outgrp_cluster(probe_fasta, contig_fasta, outgroups=['Cichorium_intybus_ERR5033750'])
    # test_gff = "/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_paralogy_2/saute/paralogy_6488_vsearch_con_aln_gff.gfa"
    # remove_gff(test_gff)

if __name__ == "__main__":

    main()