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
from io import StringIO
from itertools import groupby, chain
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import time
from typing import Union, Optional, Self, Literal, TypeAlias

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import numpy.typing as npt
from sklearn.cluster import HDBSCAN

from alignment import Alignment_fragment
from gff_parser import Fragment, Cds
from interval_tree import IntervalST, Interval
from union_find import UnionFind


# =======================================================================================
#               CLASSES
# =======================================================================================


class WrongAlphabet(Exception):
    """Exception raised when no assembled contig overlap with the probe"""


def validate_sequence(
    seq: Seq, record_name: str, alphabet: Literal["dna", "protein"] = "dna"
) -> None:
    """
    Validate if a biopython Seq object is a dna or protein sequence
    :param seq:
    :param alphabet:
    :return: None, raises WrongAlphabet exception if failed
    """
    alphabets = {
        "dna": re.compile("^[acgtn]*$", re.I),
        "protein": re.compile("^[acdefghiklmnpqrstvwy*]*$", re.I),
    }

    if alphabets[alphabet].search(str(seq)) is not None:
        return
    else:
        raise WrongAlphabet(f"[error] The sequence {record_name} is not {alphabet}")


def run_miniprot(
    probe_path: Path,
    contig_path: Path,
    treads=2,
    min_similarity=0.85,
    min_coverage=0.01,
) -> str:
    """
    Wrapper for running miniprot.
    :return:
    """

    miniprot_cmd = f"miniprot -t{treads} --gff --aln --outn 1 --outs {min_similarity} --outc {min_coverage} {contig_path} {probe_path}".split()

    try:
        miniprot = subprocess.run(
            miniprot_cmd,
            timeout=1000,
            check=True,
            capture_output=True,
            encoding="utf-8",
        )
    except FileNotFoundError as exc:
        print(f"Process failed because the executable could not be found.\n{exc}")
        raise
    except subprocess.CalledProcessError as exc:
        print(
            f"Process failed because did not return a successful return code. "
            f"Returned {exc.returncode}\n{exc}"
        )
        raise
    except subprocess.TimeoutExpired as exc:
        print(f"Process timed out.\n{exc}")
        raise
    return StringIO(miniprot.stdout)


@dataclass
class HDBcluster:
    """
    Data structure for the fragment object produced by the parsing of the miniprot output.
    Each fragment is one exon.
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
    -keep_singletons: if true, keep unclustered contigs.

    Post Init
    -probe_record: parsed probe sequence.
    -contigs_dict: dictionary of contig name as key and SeqRecord as value.
    -cds_dict: dictionary of contig name as key and Cds as value.
    -distance_matrix
    """

    probe_fasta: str
    contigs_fasta: str
    treads: int = field(default=8)
    min_probe_contig_sim: float = field(default=0.90)
    min_fragment_cov: float = field(default=0.05)
    min_contig_overlap: int = field(default=300)
    max_clustering_dist: float = field(default=0.05)
    min_cluster_size: int = field(default=3)
    min_contig_length: int = field(default=300)
    keep_singletons: bool = field(default=True)
    probe_record: SeqRecord = field(init=False, repr=False)
    contigs_dict: dict[str, list[SeqRecord]] = field(
        init=False, repr=False, default_factory=dict
    )
    cds_dict: dict[str, Cds] = field(init=False, repr=False, default_factory=dict)
    distance_matrix: npt.ArrayLike = field(init=False, repr=False, default=None)
    clusters: list[list[str]] = field(init=False, repr=False, default_factory=list)

    def __post_init__(self):
        self.probe_fasta = Path(self.probe_fasta)
        self.contigs_fasta = Path(self.contigs_fasta)
        assert self.probe_fasta.exists(), f"{self.probe_fasta} does not exist"
        assert self.contigs_fasta.exists(), f"{self.contigs_fasta} does not exist"
        self._parse_fasta()
        print("Running miniprot")
        with tempfile.TemporaryDirectory() as tmpdirname:
            for contig, record in self.contigs_dict.items():
                fasta_name = f"{contig}.fas"
                contig_fasta = Path(tmpdirname) / fasta_name
                SeqIO.write(record, contig_fasta, "fasta")
                miniprot_out = run_miniprot(
                    self.probe_fasta,
                    contig_fasta,
                    self.treads,
                    self.min_probe_contig_sim,
                    self.min_fragment_cov,
                )
                self.cds_dict[contig] = Cds(miniprot_out)
                print(self.cds_dict[contig])
        print("building distance matrix")
        tic = time.perf_counter()
        self._create_distance_matrix()
        toc = time.perf_counter()
        print(f"Computed distances in  {toc - tic:0.4f} seconds")
        self._get_clusters()

    def _parse_fasta(self) -> None:
        """
        Load the probe fasta file and the contig fasta file.
        :return:
        """
        try:
            self.probe_record = SeqIO.read(self.probe_fasta, "fasta")
            validate_sequence(self.probe_record.seq, self.probe_fasta.name, "protein")
        except Exception as err:
            sys.exit(f"[ERROR] {err}")
        try:
            contigs_dict = SeqIO.to_dict(SeqIO.parse(self.contigs_fasta, "fasta"))
        except Exception as err:
            sys.exit(f"[ERROR] {err}")
        # cleaning the name and descriptions, required for exonerate, still necessary with miniprot?
        for key, record in contigs_dict.items():
            record.description = ""
            record.name = ""
            validate_sequence(record.seq, self.contigs_fasta.name, "dna")
            self.contigs_dict[key] = record

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
            distances = [0]
            for target in items[idx + 1 :]:
                target_name_node, target_cds = target

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
                        self.probe_record,
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


    def _cluster_unionfind(
        self, contig_names: list[str], distance_matrix: npt.ArrayLike
    ) -> list[str]:
        """
        Cluster contigs with the UnionFind algorithm. Contigs need to overlap and their pairwise distance below
        the max_clustering_dist value.
        :return: list of group labels that are the names of the contigs.
        """
        with UnionFind() as UV:
            for i in range(len(contig_names)):
                for j in range(i + 1, len(contig_names)):
                    distance = distance_matrix[i, j]
                    if distance < self.max_clustering_dist:
                        UV.union((i, j))
        components = UV.get_components()
        final_groups = []
        for component in components:
            final_groups.append(sorted([list(contig_names)[c] for c in component]))
        if self.keep_singletons:
            singletons = set(contig_names) - set(list(chain(*final_groups)))
            final_groups.extend([singleton] for singleton in singletons)
        return sorted(final_groups)

    def _cluster_hdbscan(
        self, contig_names: list[str], distance_matrix: npt.ArrayLike
    ) -> list[str]:
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
        ).fit(distance_matrix)
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
        if self._use_UF():
            print("running UF")
            return self._cluster_unionfind(self.contigs_dict, self.distance_matrix)

        print("running HDBSCAN")
        clusters = self._cluster_hdbscan(self.contigs_dict, self.distance_matrix)
        final_clusters = []
        for cluster in clusters:
            if len(cluster) == 1:
                final_clusters.append(cluster)
                continue
            # Double check that each sequence in a given cluster found by HDBscan
            # does overlap with at least another sequence from the same cluster.
            # This check is necessary in case of clusters made out of out-layers.
            new_clusters = self._cluster_unionfind(cluster, self.distance_matrix)

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

    def _trim_msa(self, cluster_names: list[str]) -> Optional[list[SeqRecord]]:
        """
        Given a list of records, trim them to the smallest region of the probe shared by at least two sequences.
        In case, the group contains a single member, the sequence is trimmed to fit the probe's boundaries.
        :param cluster_names:
        :return:
        """
        if len(cluster_names) == 1:
            contig = cluster_names[0]
            cds = self.cds_dict[contig]
            record = self.contigs_dict[contig]
            if cds.is_empty():
                return
            if cds.get_global_sim() < self.min_probe_contig_sim:
                return
            strand = cds.fragments[0].get_strand()
            if strand == 1:
                trim_start, trim_end = cds.mRNA_start, cds.mRNA_end
            else:
                trim_start, trim_end = cds.mRNA_end, cds.mRNA_start
            if abs(trim_end - trim_start) < self.min_contig_length:
                return

            if strand == -1:
                new_record = SeqRecord(
                    record.seq[trim_end:trim_start].reverse_complement(),
                    id=contig,
                    name="",
                    description="",
                )
            else:
                new_record = SeqRecord(
                    record.seq[trim_start:trim_end],
                    id=contig,
                    name="",
                    description="",
                )
            return [new_record]

        trimmed_records = []
        interval_tree = IntervalST()
        #  Find overlapping edges.
        for contig in cluster_names:
            cds = self.cds_dict[contig]
            if cds.is_empty():
                continue
            fragments = cds.get_fragments()
            start = list(fragments[0].get_correspondence().keys())[0]
            end = list(fragments[-1].get_correspondence().keys())[-1]
            interval_tree.put(Interval(start, end), contig)
        common_intervals = interval_tree.get_all_intersections()
        probe_start = min([overlap.get_start() for overlap in common_intervals])
        probe_end = max([overlap.get_end() for overlap in common_intervals])

        #  Trim the sequences
        for contig in cluster_names:
            cds = self.cds_dict[contig]
            if cds.is_empty():
                continue
            record = self.contigs_dict[contig]
            fragments = cds.get_fragments()
            strand = cds.fragments[0].get_strand()
            contig_start = list(fragments[0].get_correspondence().keys())[0]
            contig_end = list(fragments[-1].get_correspondence().keys())[-1]
            if strand == 1:
                trim_start, trim_end = cds.mRNA_start, cds.mRNA_end
            else:
                trim_start, trim_end = cds.mRNA_end, cds.mRNA_start
            # Case contig within core region
            if probe_start < contig_start < contig_end < probe_end:
                pass
            if contig_start <= probe_start:
                first, last = get_nuc_coordinates(probe_start, contig_end, fragments)
                if None in [first, last]:
                    continue
                trim_start = first if strand == 1 else last

            if probe_end <= contig_end:
                first, last = get_nuc_coordinates(
                    contig_start, probe_end - 1, fragments
                )
                if None in [first, last]:
                    continue
                trim_end = last if strand == 1 else first

            if abs(trim_end - trim_start) < self.min_contig_length:
                continue

            if strand == -1:
                new_record = SeqRecord(
                    record.seq[trim_end:trim_start].reverse_complement(),
                    id=contig,
                    name="",
                    description="",
                )
            else:
                new_record = SeqRecord(
                    record.seq[trim_start:trim_end],
                    id=contig,
                    name="",
                    description="",
                )
            trimmed_records.append(new_record)
        return trimmed_records


def get_nuc_coordinates(
    probe_start: int, probe_end: int, fragments: list[Fragment]
) -> Optional[list[int, int]]:
    """
    Given a start and end position on a probe in AA space, return the nucleotide coordinates of start and end.
    :param probe_start:
    :param probe_end:
    :param fragments: list of fragments
    :return: Return [None, None] if for some reason there are no AA-nt correspondence, otherwise [start, end]
    """
    correspondences = {}
    for fragment in fragments:
        correspondences.update(fragment.get_correspondence())
    while probe_start < probe_end:
        if (first := correspondences.get(probe_start)) is None:
            probe_start += 1
        if (last := correspondences.get(probe_end)) is None:
            probe_end -= 1
        if None not in [first, last]:
            return sorted(
                [
                    first,
                    last,
                ]
            )
    return [None, None]


def main():
    ...


if __name__ == "__main__":
    main()
