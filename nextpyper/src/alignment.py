#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys

from collections import defaultdict
from dataclasses import dataclass, field
import os
from pathlib import Path
import pickle
from typing import Final, Optional, Self, TypedDict, Literal, Any

from Bio import Align
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
from gff_parser import Fragment, Cds
from interval_tree import IntervalST, Interval, Overlap


@dataclass
class Alignment_fragment:
    """
    Data structure for aligning two contigs over their exonic regions.
    Attributes
    ----------
        -cds_0: Cds object of the first contig.
        -cds_1: Cds object of the second contig.
        -probe: Seqrecord of the probe in amino acids.
        -min_overlap_nt: Minimum overlap length to compute a distance.
        -dna_model: DNA model used for the distance calculation.
    Post Init
        -fragments_0: fragment objects of the first contig.
        -fragments_1: fragment objects of the second contig.
        -contig_0: name of the first contig.
        -contig_1: name of the second contig.
        -common_intervals: List of Overlap objects in the form of
            Overlap(start=1434, end=1469, first_segment=('NODE_2_length_3201_cov_9.143787', 0),
            second_segment=('NODE_1_length_3428_cov_6.322630', 0)) where in the segment tuple,
            we have the name and the index in list of fragments of the matching fragment.
        -aln: Concatenated alignment of all fragments.
        -similarity: similarity computed on the pairwise alignment (max 1, min 0).
    """

    cds_0: Cds
    cds_1: Cds
    probe: SeqRecord
    min_overlap_nt: int
    dna_model: DistanceCalculator.dna_models = "identity"
    fragments_0: list[Fragment] = field(init=False, repr=False)
    fragments_1: list[Fragment] = field(init=False, repr=False)
    contig_0: str = field(init=False, repr=False)
    contig_1: str = field(init=False, repr=False)
    common_intervals: list[Overlap] = field(init=False, repr=False)
    fraction_probe_covered: int = field(init=False, repr=False, default=0)
    aln: MultipleSeqAlignment = field(init=False, repr=False, default=None)
    similarity: Optional[float] = field(init=False, repr=False, default=None)

    def __post_init__(self):
        self._find_common_intervals()
        self._produce_alignment()
        if self.aln is not None:
            self._compute_similarity()

    def _find_common_intervals(self) -> Self:
        """
        Create an interval tree from the probe fragments that match the contigs.
        Each fragment limits are defined as the first and last index of the coordinate mapping.
        Use a bfs search to find all overlapping intervals.
        :return:
        """
        interval_tree = IntervalST()

        def add_to_tree(cds_name, fragments: list[Fragment]):
            for idx, frg in enumerate(fragments):
                correspondence = frg.get_correspondence()
                start = list(correspondence.keys())[0]
                end = list(correspondence.keys())[-1]
                interval_tree.put(Interval(start, end), (cds_name, idx))

        self.fragments_0 = self.cds_0.get_fragments()
        self.fragments_1 = self.cds_1.get_fragments()
        self.contig_0 = self.fragments_0[0].get_contig_name()
        self.contig_1 = self.fragments_1[0].get_contig_name()
        add_to_tree(self.contig_0, self.fragments_0)
        add_to_tree(self.contig_1, self.fragments_1)
        self.common_intervals = interval_tree.get_all_intersections()
        for interval in self.common_intervals:
            self.fraction_probe_covered += interval.get_length()
        return self

    def _extract_nt_region(
        self,
        nt_start: int,
        nt_end: int,
        strand: Literal[-1, 1],
        contig_name: str,
    ) -> SeqRecord:
        """
        Extract the nucleotide sequence of the matching region for a given fragment.
        :param nt_start: alignment start position on the probe.
        :param nt_end: alignment end position on the probe.
        :param strand: matching strand.
        :param contig_name: name of contig.
        :return: sequence of contig matching the probe.
        """
        if contig_name == self.contig_0:
            cds = self.cds_0
        else:
            cds = self.cds_1
        contig_record = SeqRecord(
            Seq(cds.target_nucleotides.upper()), id=contig_name, name="", description=""
        )

        if strand == -1:
            nt_start_minus = cds.mRNA_end - nt_end
            nt_end_minus = cds.mRNA_end - nt_start
            return SeqRecord(
                contig_record.seq[nt_start_minus:nt_end_minus],
                id=contig_name,
                name="",
                description="",
            )
        nt_start = nt_start - cds.mRNA_start
        nt_end = nt_end - cds.mRNA_start
        return SeqRecord(
            contig_record.seq[nt_start:nt_end],
            id=contig_name,
            name="",
            description="",
        )

    def _get_fragment(self, contig_name: str, fragment_number: int) -> Fragment:
        """
        Given a contig name and fragment number, return the corresponding fragment object.
        :param contig_name:
        :param fragment_number: idx of the fragment in the fragment list.
        :return:
        """
        if contig_name == self.contig_0:
            return self.fragments_0[fragment_number]
        return self.fragments_1[fragment_number]

    def _get_nuc_coordinates(
        self, probe_start: int, probe_end: int, fragment: Fragment
    ) -> Optional[list[int, int]]:
        """
        Given a start and end position on a probe in AA space, return the nucleotide coordinates.
        :param probe_start:
        :param probe_end:
        :param fragment:
        :return:
        """
        # print(f"{probe_start=}, {probe_end=}, {fragment.correspondence=}")
        while probe_start < probe_end:
            if (first := fragment.correspondence.get(probe_start)) is None:
                probe_start += 1
            if (last := fragment.correspondence.get(probe_end)) is None:
                probe_end -= 1
            if None not in [first, last]:
                return sorted(
                    [
                        first,
                        last,
                    ]
                )
        return [None, None]

    def _produce_alignment(self) -> Self:
        """
        Populate the aln object. Extract all sequences from overlapping exons that match the probe sequence.
        Align them separately and concatenate the aligned sequences.
        :return:
        """
        if not self.common_intervals:
            self.aln = None
            return self
        aligned_intervals = defaultdict(list)
        for overlap in self.common_intervals:
            probe_start = overlap.start
            probe_end = overlap.end
            first_contig_name = overlap.first_segment[0]
            first_fragment_nbr = overlap.first_segment[1]
            second_contig_name = overlap.second_segment[0]
            second_fragment_nbr = overlap.second_segment[1]
            if first_contig_name == self.contig_0:
                frg_0 = self._get_fragment(first_contig_name, first_fragment_nbr)
                frg_1 = self._get_fragment(second_contig_name, second_fragment_nbr)
            else:
                frg_0 = self._get_fragment(second_contig_name, second_fragment_nbr)
                frg_1 = self._get_fragment(first_contig_name, first_fragment_nbr)
            start_0, end_0 = self._get_nuc_coordinates(probe_start, probe_end, frg_0)
            start_1, end_1 = self._get_nuc_coordinates(probe_start, probe_end, frg_1)
            if None in [start_0, end_0, start_1, end_1]:
                self.aln = None
                return self

            record_0 = self._extract_nt_region(
                start_0,
                end_0,
                frg_0.strand,
                self.contig_0,
            )
            record_1 = self._extract_nt_region(
                start_1,
                end_1,
                frg_1.strand,
                self.contig_1,
            )

            for contig_name, sequence in pairwise_aligner(record_0, record_1).items():
                aligned_intervals[contig_name].append(sequence)
        concatenated_sequences = {
            name: "".join(sequence) for name, sequence in aligned_intervals.items()
        }
        # Remove indels
        indel_idxs = []
        for idx in range(len(list(concatenated_sequences.values())[0])):
            for seq in concatenated_sequences.values():
                if seq[idx] == "-":
                    indel_idxs.append(idx)

        clean_sequences = {
            name: "".join([i for j, i in enumerate(seq) if j not in indel_idxs])
            for name, seq in concatenated_sequences.items()
        }

        total_alignment = MultipleSeqAlignment(
            SeqRecord(Seq(sequence), id=name)
            for (name, sequence) in clean_sequences.items()
        )
        if total_alignment.get_alignment_length() > self.min_overlap_nt:
            self.aln = total_alignment
            return self
        self.aln = None
        return self

    def _compute_similarity(self):
        """
        Compute pairwise distance, normalized by fraction of probe that is covered
        :return:
        """
        calculator = DistanceCalculator(self.dna_model)
        dm = calculator.get_distance(self.aln)
        raw_distance = dm.matrix[1][0]
        assert (
            self.fraction_probe_covered != 0.0
        ), f"[ERROR] the two targeted sequence have no overlap with the probe."
        # print(f"{raw_distance=}")
        # self.similarity = 1 - ((1 - raw_distance) * self.fraction_probe_covered) / len(
        #     self.probe
        # )
        # self.similarity = 1 - raw_distance
        self.similarity = raw_distance

    def get_similarity(self):
        if self.aln is None:
            return None
        return self.similarity

    def get_alignment(self):
        if self.aln is None:
            return None
        return self.aln

    def write_alignment(self, alignment_file: Path):
        if self.aln is None:
            pass
        else:
            AlignIO.write([self.aln], alignment_file, "fasta")


def pairwise_aligner(
    record_0: SeqRecord,
    record_1: SeqRecord,
    match_score=1,
    mismatch_score=-0.1,
    gap_score=-1,
    internal_open_gap_score=-6,
    internal_extend_gap_score=-3,
) -> dict[str, str]:
    aligner = Align.PairwiseAligner()
    # Increase the cost of gap
    aligner.match_score = match_score
    aligner.mismatch_score = mismatch_score
    aligner.gap_score = gap_score
    aligner.internal_open_gap_score = internal_open_gap_score
    aligner.internal_extend_gap_score = internal_extend_gap_score
    alignment = aligner.align(record_0, record_1)[0]
    # print({record_0.id: alignment[0], record_1.id: alignment[1]})
    return {record_0.id: alignment[0], record_1.id: alignment[1]}


def main():
    pickle_file = "/home/yjkbertrand/Documents/projects/nextpyper/test_data/test_clustering/gene_8631_6contigs.pkl"
    # pickle_file = "/home/yjkbertrand/Documents/projects/nextpyper/test_data/test_clustering/gene_8631_all.pkl"

    cds_dict = pickle.load(
        open(
            pickle_file,
            "rb",
        )
    )
    probe = "/home/yjkbertrand/Documents/projects/nextpyper/test_data/test_clustering/probe_8631_aa.fasta"
    probe_record = list(SeqIO.parse(probe, "fasta"))[0]
    print(cds_dict)

    cdss = list(cds_dict.items())
    from itertools import combinations

    idx = 0
    for comb in list(combinations(cdss, 2)):
        print(comb)
        node_0 = comb[0][1]
        node_1 = comb[1][1]
        if any([node_0.is_empty(), node_1.is_empty()]):
            continue
        aln_frag = Alignment_fragment(
            node_0,
            node_1,
            probe_record,  # dna_model="identity"
        )
        # print(aln_frag.get_similarity())

        if aln_frag.get_similarity():
            fasta = f"test_aln_{idx}.fasta"

            aln_frag.write_alignment(
                f"/home/yjkbertrand/Documents/projects/nextpyper/test_data/test_clustering/{fasta}"
            )
            print(f"saved {fasta}")
            idx += 1

    # node_0 = cdss[0]
    # node_1 = cdss[1]
    # aln_frag = Alignment_fragment(node_0, node_1, probe_record, dna_model="identity")
    # print(aln_frag.get_distance())
    # print("writing alignment")
    # aln_frag.write_alignment(
    #     "/home/yjkbertrand/Documents/projects/nextpyper/test_data/test_clustering/test_aln.fasta"
    # )


if __name__ == "__main__":
    main()
