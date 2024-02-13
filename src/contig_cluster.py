#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2024
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.
"""
Classes used to cluster contigs according to their pairwise distance.
The pairwise alignments contigs are obtained by first aligning them to a probe.
The probe is a single amino acid sequence in fasta format, while the contigs
are typically the assembly of SPAdes.
The common regions of alignment to the probe are used to compute their pairwise
distances. Using these distances, contigs are clustered using the DBSCAN algorithm.
For each cluster, the sequence composing it can be written in fasta format.

Example usage:

probe_fasta = Path("probe_3_aa.fasta")
contig_fasta = Path("gene_3_all.fasta")

cluster = DBcluster(probe_fasta, contig_fasta, False, True)
print("getting the matrix")
print(cluster.get_matrix())
knee = cluster.get_knee()
if isinstance(knee, float):
    print(f"labels from the knee as {type(knee)=}")
    labels = cluster.cluster(knee)
else:
    print("[ERROR] No value was found for the knee.")
    # print("labels from fixed value")
    # labels = cluster.cluster(0.051773049645390146)

cluster.save_cluster(
    labels,
    "/home/yjkbertrand/Documents/projects/nextpiper/test_data/batrachium/exonerate/clusters_2",
)
"""


from pathlib import Path
import subprocess
import os
from io import StringIO
from collections import defaultdict
from itertools import groupby
import sys

from dataclasses import dataclass, field
from typing import Final, Optional, Self, TypedDict, Literal, Any, Self
from icecream import ic

import numpy as np
import numpy.typing as npt
from collections import deque
import math
from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors
import matplotlib.pyplot as plt

from Bio import Align
from Bio import SearchIO
from Bio.SearchIO._model.hsp import HSPFragment
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import seq1
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator

from kneed import KneeLocator


@dataclass
class Probe_limits:
    """
    A dataclass used to store the coordinates of the overlap between two contigs.
    Because of indels in the probe sequence, we cannot directly deduce the contigs overlaps
    from their alignment coordinates on the probe.
    Attributes
    ----------
        -current_left: index of 5' overlap with the probe.
        -current_right: index of 3' overlap with the probe.
    Post init
        -final_left: index of 5' overlap with the probe at the end of the search.
        -final_right: index of 3' overlap with the probe at the end of the search.
    """

    current_left: int
    current_right: int
    final_left: int = field(init=False, default=None)
    final_right: int = field(init=False, default=None)

    def is_complete(self) -> bool:
        """
        Check whether an index has been found both for the 'final_left' and 'final_right' attributes.
        :return:
        """
        if None in [self.final_left, self.final_right]:
            return False
        return True

    def seach_fail(self) -> Self:
        """
        When no good value has been found for the 'final_left' and 'final_right' attributes,
        they are set to the same index as 'current_left' and 'current_right'.
        :return:
        """
        self.final_left = self.current_left
        self.final_right = self.current_right
        return self

    def get_final_left(self) -> int:
        return self.final_left

    def get_final_right(self) -> int:
        return self.final_right


# for each hsp we create an object holding the alignment and information about it
@dataclass
class Alignment_fragment:
    """
    Data structure for the fragment object produced by the parsing of the Exonerate output
    Attributes
    ----------
        -hsp_fragment: fragment obtained from Exonerate.
        -probe_record: aa record of the probe.
        -contig_record: nt record of the contig.
    Post Init
        -protein_model: model used to compute aa sequence similarity.
        -aln: protein-protein alignment.
        -query_id: header of the probe.
        -hit_id: header of the contig.
        -query_start: index in nt of the start of the alignment block on the probe.
        -query_end: index in nt of the end on the probe.
        -hit_start: index in nt of the start on the contig.
        -hit_end: index in nt of the end on the contig.
        -hit_orientation: direct or reversed complemented.
        -aa_similarity: similarity on the aligned block.
        -correspondence: dictionary in nt (!) that for each position on the probe assign the position on the contig.
    """

    hsp_fragment: HSPFragment = field(repr=False)
    probe_record: SeqRecord = field(repr=False)
    contig_record: SeqRecord = field(repr=False)
    protein_model: Literal[
        "blastp",
        "blosum45",
        "blosum50",
        "blosum62",
        "blosum80",
        "blosum90",
        "pam250",
        "pam30",
        "pam70",
    ] = field(
        default="blosum62", repr=False
    )  # for similarity computation
    aln: MultipleSeqAlignment = field(init=False, repr=False)
    query_id: str = field(init=False)
    hit_id: str = field(init=False)
    query_start: int = field(init=False)  # indexed on protein!
    query_end: int = field(init=False)  # indexed on nucleotide!
    hit_start: int = field(init=False)
    hit_end: int = field(init=False)
    hit_orientation: Literal[-1, 1] = field(init=False)
    aa_similarity: float = field(init=False)
    correspondence: dict[int, int] = field(init=False, repr=False)

    def __post_init__(self):
        # self.aln = self.hsp_fragment.aln
        self.query_id = self.hsp_fragment.query_id
        self.hit_id = self.hsp_fragment.hit_id
        self._get_query_start_end()
        self.hit_start = self.hsp_fragment.hit_start
        self.hit_end = self.hsp_fragment.hit_end
        self.hit_orientation = self.hsp_fragment.hit_strand
        self._get_alignment()
        self._get_similarity()
        self._get_correspondence()

    def _get_query_start_end(self) -> Self:
        """
        For some reason the starting index on the probe sequence is not always accurate and
        might not be divisible by 3. We take the nearest higher integer that is divisible by 3 in such
        case. Redo this with using the hit_frame attribute.
        There is also a problem with the start coordinate for the probe as it is given in protein space whereas
        the end is in nucleotide space.
        :return: populates self.query_start and self.query_end.
        """
        if self.hsp_fragment.query_start == 0:
            query_start = self.hsp_fragment.query_start
            query_end = self.hsp_fragment.query_end
        else:
            query_start = self.hsp_fragment.query_start - 1
            query_end = self.hsp_fragment.query_end - 1
        remainder_start = query_start % 3
        remainder_end = query_end % 3
        self.query_start = (
            query_start * 3
            if not remainder_start
            else (query_start - remainder_start + 3) * 3
        )
        self.query_end = (
            query_end if not remainder_end else query_end - remainder_end + 3
        )
        return self

    def get_hit_length(self) -> int:
        """
        Because reversed in complemented hits the start and the end of the hit are in
        decreasing order, we need to take the absolute value when computing the length.
        :return: the length of the overlapping region of the contig with the probe.
        """
        return abs(self.hit_end - self.hit_start)

    def _get_alignment(self) -> Self:
        """
        Convert aln from three letters code to one letter
        :return: populate the self.aln attribute.
        """
        new_records = []
        for record in self.hsp_fragment.aln:
            record.seq = seq1(record.seq)
            new_records.append(record)
        self.aln = MultipleSeqAlignment(new_records)
        return self

    def _get_similarity(self) -> Self:
        """
        Compute nucleotide pairwise similarity (i.e. 1 - distance)
        :return: populate the self.aa_similarity attribute.
        """
        assert len(self.aln) == 2, "[Error] The alignment must contain two sequences"
        calculator = DistanceCalculator(self.protein_model)
        dm = calculator.get_distance(self.aln)
        distance = dm.matrix[1][0]
        self.aa_similarity = 1 - distance
        return self

    def _get_correspondence(self) -> Self:
        """
        Create a dictionary where keys are probe indexes and values are the corresponding position on the contig.
        Indels are accounted for.
        :return: populates the self.aa_correspondence attribute.
        """
        corr = {}
        probe_idx = self.query_start
        contig_idx = self.hit_start if self.hit_orientation == 1 else self.hit_end

        for idx in range(self.aln.get_alignment_length()):
            column = self.aln[:, idx]
            if column[0] == "X":
                contig_idx += 3 * self.hit_orientation
            elif column[1] == "X":
                probe_idx += 3
            else:
                corr[probe_idx] = contig_idx
                probe_idx += 3
                contig_idx += 3 * self.hit_orientation
        self.correspondence = corr
        return self

    def get_probe_length(self) -> int:
        return len(self.probe_record.seq)

    def get_contig_record(self) -> SeqRecord:
        """
        Extract the contig record, reorient it in case of reversed complement.
        :return:
        """
        if self.hit_orientation == 1:
            return self.contig_record
        else:
            record = self.contig_record.reverse_complement()
            record.id = self.contig_record.id
            record.description = "[Reverse Complement]"
            record.name = ""
            return record


@dataclass
class Pairwise_aln:
    """
    Data structure for performing contig to contig pairwise alignments.
    Attributes
    ----------
        -fragment_0: Alignment_fragment object for the first contig.
        -fragment_1: Alignment_fragment object for the second contig.
        -min_overlap_nt: minimum overlap between the two contigs for computing their distance.
        -dna_model: model for computing the distance between contigs.
        see https://biopython.org/docs/1.75/api/Bio.Phylo.TreeConstruction.html for a list of possible models.
    Post Init
        -fragments_aln: pairwise alignment.
        -distance: distance between the two aligned contigs.
    """

    # contig to contig alignment
    fragment_0: Alignment_fragment
    fragment_1: Alignment_fragment
    min_overlap_nt: int = 300
    dna_model: DistanceCalculator.dna_models = "blastn"
    fragments_aln: Optional[MultipleSeqAlignment] = field(
        init=False, repr=False, default=None
    )
    distance: Optional[float] = field(init=False, repr=False, default=None)

    def __post_init__(self):
        assert (
            self.dna_model in DistanceCalculator.dna_models
        ), f"[ERROR] {self.dna_nodel} is not in {DistanceCalculator.dna_models}"
        self._create_nt_aln()

    def _create_nt_aln(self) -> Self:
        """
        Align the two sequences over their overlapping fragment (transitively obtained from the alignment to the probe).
        Compute their pairwise distance.
        :return: populates the self.distance attribute.
        """
        start_candidate = max(self.fragment_0.query_start, self.fragment_1.query_start)
        end_candidate = min(self.fragment_0.query_end, self.fragment_1.query_end) - 3
        limits = self._find_correspondence(start_candidate, end_candidate)
        start = limits.get_final_left()
        end = limits.get_final_right()
        overlap = end - start
        if overlap > self.min_overlap_nt:
            seq_0 = self._extract_nt_region(self.fragment_0, start, end)
            seq_1 = self._extract_nt_region(self.fragment_1, start, end)
            assert all(
                [seq_0, seq_1]
            ), f"[Error] One or both sequences are of zero length {seq_0=} , {seq_1=}"
            aligner = Align.PairwiseAligner()
            # Increase the cost of gap
            aligner.match_score = 1
            aligner.mismatch_score = -0.1
            aligner.gap_score = -1
            aligner.internal_open_gap_score = -6
            aligner.internal_extend_gap_score = -3
            alignment = aligner.align(seq_0, seq_1)[0]
            record_0 = SeqRecord(Seq(alignment[0]), id=self.fragment_0.hit_id)
            record_1 = SeqRecord(Seq(alignment[1]), id=self.fragment_1.hit_id)
            self.fragments_aln = MultipleSeqAlignment([record_0, record_1])
            # Write the alignments for Debugging
            # from Bio import AlignIO
            # import shutil
            # name = self.fragment_0.hit_id[:13] + "_" + self.fragment_1.hit_id[:13] + ".fasta"
            # AlignIO.write(self.fragments_aln, name, "fasta")
            # shutil.move(name,
            #             "/home/yjkbertrand/Documents/projects/nextpiper/test_data/batrachium/exonerate/alignments")
            calculator = DistanceCalculator(self.dna_model)
            dm = calculator.get_distance(self.fragments_aln)
            # normalize distance by alignment length fraction of the probe length
            raw_distance = dm.matrix[1][0]
            self.distance = (
                raw_distance * self.fragments_aln.get_alignment_length()
            ) / self.fragment_0.get_probe_length()
            # print(f"{self.distance=}")

    def _extract_nt_region(
        self, frag: Alignment_fragment, probe_start: int, probe_end: int
    ) -> str:
        """
        Extract matching region
        :param frag: Alignment_fragment object.
        :param probe_start: alignment start position on the probe.
        :param probe_end: alignment end position on the probe.
        :return: sequence of contig matching the probe.
        """
        record = frag.contig_record
        assert (
            frag.hit_id == record.id
        ), f"[ERROR] Mismatch between fragment {frag.hit_id} and sequence {record.id}"
        first = frag.correspondence[probe_start]
        end = frag.correspondence[probe_end]
        if frag.hit_orientation == -1:
            seq_region = str(record.seq[end:first].reverse_complement())
        else:
            seq_region = str(record.seq)[first:end]
        return seq_region

    def _find_correspondence(self, left: int, right: int) -> Probe_limits:
        """
        Because of deletion in the contig or insertions in the probe sequence,
        not all probe indices have a counterpart contig index. We search for the first
        counterpart by shrinking the overlapping region.
        :param left: start of the overlap between the contigs
        :param right: end of the overlap between the contigs
        :return: Probe_limits object
        """
        i = 0
        limits = Probe_limits(left, right)
        while limits.current_left < limits.current_right and i < 5:
            i += 1
            if limits.is_complete():
                return limits
            if limits.final_left is None:
                contig_left_0 = self.fragment_0.correspondence.get(limits.current_left)
                contig_left_1 = self.fragment_1.correspondence.get(limits.current_left)
                if None not in (contig_left_0, contig_left_1):
                    limits.final_left = limits.current_left
                else:
                    limits.current_left += 3
            if limits.final_right is None:
                contig_right_0 = self.fragment_0.correspondence.get(
                    limits.current_right
                )
                contig_right_1 = self.fragment_1.correspondence.get(
                    limits.current_right
                )
                if None not in (contig_right_0, contig_right_1):
                    limits.final_right = limits.current_right
                else:
                    limits.current_right -= 3
        limits.seach_fail()
        return limits

    def get_distance(self) -> float:
        return self.distance


@dataclass
class DBcluster:
    probe_fasta: Path
    contigs_fasta: Path
    compute_exonerate: bool = field(default=True)
    load_matrix: bool = field(default=False)
    probe_record: SeqRecord = field(init=False, repr=False)
    contig_dict: dict[int, list[SeqRecord]] = field(init=False, repr=False)
    clean_fragments: defaultdict = field(init=False, repr=False)
    distance_matrix: npt.ArrayLike = field(init=False, repr=False, default=None)

    def __post_init__(self):
        self._parse_fasta()
        if self.compute_exonerate:
            try:
                print("computing exonerate")
                exonerate_output = self._run_exonerate()
            except Exception as err:
                sys.exit(f"[ERROR] {err}")
        else:
            exonerate_file = str(self.contigs_fasta.stem) + ".exonerate"
            print("loading exonerate")
            exonerate_output = list(SearchIO.parse(exonerate_file, "exonerate-text"))[0]
        assert (
            exonerate_output
        ), f"[ERROR] Exonerate produced no output with probe file {self.probe_fasta} and contig file {self.contigs_fasta}"

        contig_exonerate = defaultdict(list)
        print("parsing exonerate")
        for hit in exonerate_output:
            hsps = hit.hsps
            for hsp in hsps:
                for fragment in hsp:
                    contig_exonerate[fragment.hit_id].append(
                        Alignment_fragment(
                            fragment,
                            self.probe_record,
                            self.contig_dict[fragment.hit_id],
                        )
                    )
        print("cleaning fragments")
        self.clean_fragments: defaultdict = self._filter_fragments(contig_exonerate)
        print("building distance matrix")
        matrix_file = self.contigs_fasta.stem + "_distance_matrix.npy"
        self.matrix_path = self.contigs_fasta.parent / matrix_file
        if self.load_matrix and self.matrix_path.exists():
            print(f"loading distance matrix from {self.matrix_path}")
            self.distance_matrix = np.load(self.matrix_path)
        else:
            self.distance_matrix = self._create_distance_matrix(self.clean_fragments)
            self.save_distance_matrix()

    def save_distance_matrix(self) -> None:
        """
        Save the distance matrix to disk.
        :return:
        """
        if self.distance_matrix is not None:
            print(f"saving distance matrix as {self.matrix_path}")
            np.save(self.matrix_path, self.distance_matrix)

    def _parse_fasta(self) -> None:
        """
        Load the probe fasta file and the contig fasta file.
        :return:
        """
        try:
            self.probe_record = SeqIO.read(self.probe_fasta, "fasta")
        except Exception as err:
            sys.exit(f"[ERROR] {err}")
        try:
            contig_dict = SeqIO.to_dict(SeqIO.parse(self.contigs_fasta, "fasta"))
        except Exception as err:
            sys.exit(f"[ERROR] {err}")
        self.contig_dict = {}
        for key, record in contig_dict.items():
            record.description = ""
            record.name = ""
            self.contig_dict[key] = record

    def _filter_fragments(
        self, contig_exonerate
    ) -> defaultdict[str, Alignment_fragment]:  ### add parameters
        """
        Filter Alignment_fragment according to similarity to the probe and to the length of the
        overlap.
        :param contig_exonerate:
        :return: a dictionary with name of contig as key and filtered Alignment_fragment as value.
        """
        clean_contig_exonerate = defaultdict(list)
        for contig, fragments in contig_exonerate.items():
            for fragment in fragments:
                if filter_fragment(fragment, self.probe_record):
                    clean_contig_exonerate[contig].append(fragment)
        return clean_contig_exonerate

    def _run_exonerate(self) -> list:
        """
        Wrapper for running exonereate. The output is parsed.
        :return: list of exonerate objects
        """
        exonerate_cmd = f"exonerate -m protein2genome --showvulgar no -V 0 --refine full {self.probe_fasta} {self.contigs_fasta}".split()
        # print(f"running exonerete command {exonerate_cmd}")
        try:
            exonerate = subprocess.run(
                exonerate_cmd,
                timeout=100,
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
        fasta_io = StringIO(exonerate.stdout)
        analysis_outcome = list(next(SearchIO.parse(fasta_io, "exonerate-text")))
        fasta_io.close()
        return analysis_outcome

    def _create_distance_matrix(
        self, contig_exonerate_dict: [str, Alignment_fragment]
    ) -> np.array:
        """
        Compute the distances between all fragments from different contigs. If the fragments overlap
        a distance is calculated, otherwise it is set to None.
        :param contig_exonerate_dict: Dictionary with name of the sequence as key, Alignment_fragment
        as value.
        :return: A distance matrix.
        """

        def select_fragment(list_frag: list[Alignment_fragment]) -> Alignment_fragment:
            # select the longest (or the most similar of the fragments - not implemented-)
            if len(list_frag) == 1:
                return list_frag[0]
            else:
                return sorted(list_frag, key=lambda x: x.get_hit_length())[-1]

        items = list(contig_exonerate_dict.items())
        lists = deque()
        for idx, item in enumerate(items):
            query_name_node, query_fragments = item
            query_fragment = select_fragment(query_fragments)
            distances = [0]
            for target in items[idx + 1 :]:
                target_name_node, target_fragments = target
                target_fragment = select_fragment(target_fragments)
                aln = Pairwise_aln(query_fragment, target_fragment)
                if (distance := aln.get_distance()) is None:
                    distance = math.nan
                distances.append(distance)
            lists.append(distances)

        pad = len(max(lists, key=len))
        half_matrix = np.array([[math.nan] * (pad - len(i)) + i for i in lists])
        half_matrix = np.nan_to_num(half_matrix)
        matrix = half_matrix + half_matrix.T - np.diag(np.diag(half_matrix))
        return matrix

    def get_matrix(self):
        return self.distance_matrix

    def get_knee(self) -> Any:
        """
        Find the inflexion point using the kneedle algorithm which produced the epsilon parameter for
        the DBSCAN algorithm.
        :return: Either a numpy array when unsuccessful or a float when successful.
        """
        if self.distance_matrix is None:
            print("[ERROR] no matrix was computed")
            return
        nbrs = NearestNeighbors(n_neighbors=3, metric="precomputed").fit(
            self.distance_matrix
        )
        neigh_dist, neigh_ind = nbrs.kneighbors(self.distance_matrix)
        sort_neigh_dist = np.sort(neigh_dist, axis=0)
        i = np.arange(len(sort_neigh_dist))
        distances = sort_neigh_dist[:, 2]
        knee = KneeLocator(
            i,
            distances,
            S=1,
            curve="convex",
            direction="increasing",
            interp_method="polynomial",
            polynomial_degree=2,
        )
        k_dist = sort_neigh_dist[:, 2]
        plt.plot(k_dist)
        knee.plot_knee()
        plt.axhline(y=0.05, linewidth=1, linestyle="dashed", color="k")
        plt.ylabel("k-NN distance")
        plt.xlabel("Sorted observations (2th NN)")
        # print(f"{distances[knee.knee]=}")
        plt.show()
        return float(distances[knee.knee])

    def cluster(self, epsilon: float) -> list[int]:
        """
        Using the distance matrice and the inflection distance, cluster sequences with DBSCAN.
        :param epsilon: the knee value
        :return: list of group labels
        """
        clusters = DBSCAN(eps=epsilon, min_samples=1, metric="precomputed").fit(
            self.distance_matrix
        )
        labels = clusters.labels_
        print(clusters.labels_)
        print(f"number of clusters is {len(set(clusters.labels_))}")
        return labels

    def save_cluster(self, labels: list[int], fasta_folder: str) -> None:
        """
        Given a list of labels generated by DBSCAN, save the clusters to a fasta file in the specified folder.
        :param labels: list from DBSCAN with index position that corresponds to the contigs in the keys of
        the 'clean_fragments' dictionary.
        :param fasta_folder: folder where the fasta files are saved.
        :return:
        """

        print(f"Saving clusters to fasta folder {fasta_folder}")
        Path(fasta_folder).mkdir(parents=True, exist_ok=True)
        assert len(labels) == len(
            self.clean_fragments.keys()
        ), f"[ERROR] Number of labels {len(labels)} does not match the number of fragments {len(self.clean_fragments)}"
        record_cluster = zip(self.clean_fragments.keys(), labels)
        record_cluster_sorted = sorted(record_cluster, key=lambda x: x[1])
        groups = groupby(record_cluster_sorted, lambda x: x[1])
        record_groups = []
        for key, group in groups:
            record_groups.append(list(group))
        fasta_prefix = self.probe_fasta.stem
        for idx, gp in enumerate(record_groups):
            names = [x[0] for x in gp]
            name = f"{fasta_prefix}_{idx}.fasta"
            name_fasta = Path(fasta_folder) / name
            print(f"saving fasta {name_fasta}")
            fragments = [self.clean_fragments[key][0] for key in names]
            records = [fragment.get_contig_record() for fragment in fragments]
            SeqIO.write(records, name_fasta, "fasta")


def filter_fragment(
    alignment_fragment: Alignment_fragment,
    probes_aa: SeqRecord,
    min_probe_similarity: float = 0.85,
    min_fraction_length: float = 0.5,
    max_fraction_length: float = 3.0,
) -> bool:
    """Filter out a sequence if it is smaller than min_fraction_length*len(seq) or
    larger than max_fraction_length*len(seq) or less similar than a given threshold,
    or not a key in the exonerate dictionary."""
    assert (
        alignment_fragment.query_id == probes_aa.id
    ), f"[ERROR] Seq_0 has not been mapped against probe {probes_aa.id} but against probe {alignment_fragment.query_id}"
    # print(alignment_fragment.hit_id)
    if alignment_fragment.aa_similarity < min_probe_similarity:
        # print(f"node similarity to probe is {alignment_fragment.similarity} which is below the threshold of {min_probe_similarity}")
        return False
    probe_sequence = [x for x in alignment_fragment.aln[0].seq if x != "-"]
    probe_overlap_fraction = len(probe_sequence) / len(probes_aa)
    if not min_fraction_length < probe_overlap_fraction < max_fraction_length:
        return False
    return True


def main():
    probe_fasta = Path("probe_3_aa.fasta")
    # contig_fasta = Path("gene_3_test_clean.fasta")
    contig_fasta = Path("gene_3_all.fasta")
    # contig_fasta = Path("A11_C8.fas")

    cluster = DBcluster(probe_fasta, contig_fasta, False, True)
    print("getting the matrix")
    print(cluster.get_matrix())
    knee = cluster.get_knee()
    if isinstance(knee, float):
        print(f"labels from the knee as {type(knee)=}")
        labels = cluster.cluster(knee)
    else:
        print("[ERROR] No value was found for the knee.")
        # print("labels from fixed value")
        # labels = cluster.cluster(0.051773049645390146)

    cluster.save_cluster(
        labels,
        "/home/yjkbertrand/Documents/projects/nextpiper/test_data/batrachium/exonerate/clusters_2",
    )


if __name__ == "__main__":
    os.chdir(
        "/home/yjkbertrand/Documents/projects/nextpiper/test_data/batrachium/exonerate"
    )
    main()
