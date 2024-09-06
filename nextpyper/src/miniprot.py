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
from operator import attrgetter
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
from union_find import UnionFind


# =======================================================================================
#               EXCEPTIOMS
# =======================================================================================


class WrongAlphabet(Exception):
    """Exception raised when no assembled contig overlap with the probe"""


# =======================================================================================
#               FUNCTIONS
# =======================================================================================


def get_nuc_coordinates(
    probe_start: int, probe_end: int, fragments: list[Fragment]
) -> list[Optional[int]]:
    """
    Given a start and end position on a probe in AA space, return the nucleotide coordinates of start and end.
    :param probe_start:
    :param probe_end:
    :param fragments: list of fragments
    :return: [None, None] if for some reason there are no AA-nt correspondence, otherwise [start, end]
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

def merge_intervals(arr: list['OverlappingSeqs']):
    """
    Combine overlapping intervals of sequences and cluster their names.
    """
    # Sorting based on the increasing order of the start intervals
    arr.sort(key=attrgetter("probe_start"))
    index = 0
    index_to_del = []
    for i in range(1, len(arr)):
        if arr[index].probe_end >= arr[i].probe_start:
            arr[index].probe_end = max(arr[index].probe_end, arr[i].probe_end)
            arr[index].combine(arr[i].seq_names)
            index_to_del.append(i)
        else:
            index = index + 1
            arr[index] = arr[i]
    return [arr[index] for index in range(len(arr)) if index not in index_to_del]



def run_miniprot(
    probe_path: Path,
    contig_path: Path,
    treads=2,
    min_similarity=0.85,
    min_coverage=0.01,
) -> StringIO:
    """
    Wrapper for running miniprot.
    :return:
    """

    miniprot_cmd = f"miniprot -t {treads} --gff --aln --outn 1 --outs {min_similarity} --outc {min_coverage} {contig_path} {probe_path}".split()

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


def validate_sequence(
    seq: Seq, record_name: str, alphabet: Literal["dna", "protein"] = "dna"
) -> None:
    """
    Validate if a biopython Seq object is a dna or protein sequence
    :param seq:
    :param record_name: name of the fasta file
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


# =======================================================================================
#               CLASSES
# =======================================================================================


@dataclass
class ProbeCds:
    """
    Data structure for the fragment object produced by the parsing of the miniprot output.
    Each fragment is one exon.
    Attributes
    ----------
    -probe_fasta: path to the amino acid probes fasta file.
    -contigs_fasta: path to the nucleotide contigs fasta file.
    -treads: for miniprot.
    -min_probe_contig_sim: mapping similarity for miniprot.
    -min_fragment_cov: fraction of the probe covered by a contig for miniprot.
    -min_contig_length: minimum length of a contig after trimming, if the sequences are saved.

    Post Init
    -probe_record: parsed probe sequence.
    -contigs_dict: dictionary of contig name as key and SeqRecord as value.
    -cds_dict: dictionary of contig name as key and Cds object as value.
    """

    probes_fasta: str
    contigs_fasta: str
    treads: int = field(default=8)
    min_probe_contig_sim: float = field(default=0.85)
    min_fragment_cov: float = field(default=0.05)
    min_contig_length: int = field(default=300)
    probes_path: Path = field(init=False, repr=False)
    contigs_path:Path = field(init=False, repr=False)
    probes_dict: dict[str, list[SeqRecord]] = field(
        init=False, repr=False, default_factory=dict)
    contigs_dict: dict[str, list[SeqRecord]] = field(
        init=False, repr=False, default_factory=dict)
    cds_dict: dict[str, Cds] = field(init=False, repr=False, default_factory=dict) # defaultdict(dict)?

    def __post_init__(self):
        self.probes_path = Path(self.probes_fasta)
        self.contigs_path = Path(self.contigs_fasta)
        assert self.probes_path.exists(), f"{self.probes_fasta} does not exist"
        assert self.contigs_path.exists(), f"{self.contigs_fasta} does not exist"
        self._parse_fasta()
        print("Running miniprot")

        # run miniprot on each scaffold/probe combination
        with tempfile.TemporaryDirectory() as tmpdirname:
            # write all the probe and scaffold sequences

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

    def _parse_fasta(self) -> None:
        """
        Load the probe fasta file and the contig fasta file.
        :return:
        """
        try:
            self.probes_dict = SeqIO.to_dict(SeqIO.parse(self.probes_path, "fasta"))
        except Exception as err:
            sys.exit(f"[ERROR] {err}")
        for record in self.probes_dict.values():
            # record.description = ""
            # record.name = ""
            validate_sequence(record.seq, self.contigs_path.name, "protein")
            # self.probes_dict[key] = record
        try:
            contigs_dict = SeqIO.to_dict(SeqIO.parse(self.contigs_fasta, "fasta"))
        except Exception as err:
            sys.exit(f"[ERROR] {err}")
        # cleaning the name and descriptions, required for exonerate, still necessary with miniprot?
        for record in self.probes_dict.values():
            # record.description = ""
            # record.name = ""
            validate_sequence(record.seq, self.contigs_path.name, "dna")
            # self.contigs_dict[key] = record

    def _trim_msa(self, cluster_names: list[str]) -> Optional[list[SeqRecord]]:
        """
        Given a list of records, trim them to the smallest region of the probe shared by at least two sequences.
        In case, the group contains a single member, the sequence is trimmed to fit the probe's boundaries.
        :param cluster_names:
        :return: None is there is no match with a probe on this cluster, otherwise returns the trim sequence(s).
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


@dataclass
class OverlappingSeqs:
    """
    Container for groups of sequences that overlap a common region of the probe
    Attributes
    ----------
    -probe_start: the smallest index on the probe (AA space) that has a matching sequence.
    -probe_end: the largest index on the probe (AA space) that has a matching sequence.
    -seq_names: a list of sequence names that overlap this region.
    """

    probe_start: int
    probe_end: int
    seq_names: list[str]

    def combine(self, other: list[str]) -> Self:
        """
        Add items to seq_names
        """
        self.seq_names.extend(other)
        return self

    def get_seq_names(self) -> list[str]:
        return self.seq_names

    def get_probe_interval(self) -> str:
        """
        Get the boundaries of the probe covered by the group of sequences.
        """
        return f"start_{self.probe_start}_end_{self.probe_end}"


@dataclass
class OverlapDetect(ProbeCds):
    """
    Use miniprot to sort overlapping sequences
    """
    probe_fasta: str
    contigs_fasta: str
    treads: int = field(default=8)
    non_overlapping: list[OverlappingSeqs] = field(init=False, default_factory=list)

    def __post_init__(self):
        super().__post_init__()
        self._sorting_overlapping()
        print(self.non_overlapping)

    def _sorting_overlapping(self) -> Self:
        """
        Separate the Cds by overlap with the probe sequence.
        :return:
        """
        intervals = []
        for contig, cds in self.cds_dict.items():
            if cds.is_empty():
                continue
            start = cds.probe_start
            end = cds.probe_end
            intervals.append(OverlappingSeqs(start, end, [contig]))
        self.non_overlapping = merge_intervals(intervals)
        return self

    def save_clusters(self, fasta_folder: str) -> None:
        """
        Given a list of labels generated by HDBSCAN, save the clusters to a fasta file in the specified folder.
        :param fasta_folder: folder where the fasta files are saved.
        :return:
        """
        assert self.non_overlapping, f"[ERROR] No overlap detected in {self.contigs_fasta}"

        print(f"Saving clusters to fasta folder {fasta_folder}")
        Path(fasta_folder).mkdir(parents=True, exist_ok=True)
        fasta_prefix = self.probe_fasta.stem
        idx = 0
        for merged in self.non_overlapping:
            probe_interval = merged.get_probe_interval()
            trimmed_records = self._trim_msa(merged.get_seq_names() )
            if trimmed_records:
                name = f"{fasta_prefix}_{probe_interval}_{idx}.fasta"
                name_fasta = Path(fasta_folder) / name
                SeqIO.write(trimmed_records, name_fasta, "fasta")
                print(f"saving fasta {name_fasta}")
                idx += 1



