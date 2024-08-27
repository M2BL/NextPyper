#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2024
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.
"""

"""

__version__ = "0.1"

# =======================================================================================
#               IMPORTS
# =======================================================================================

from collections import defaultdict, Counter
from dataclasses import dataclass, field
from importlib import reload
from io import StringIO
from itertools import groupby
from operator import attrgetter
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import textwrap
import time
from typing import Optional, Self, Literal

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment

from contig_cluster import ProbeCds
import gff_parser

gff_parser = reload(gff_parser)
from gff_parser import Fragment, Cds

# =======================================================================================
#               EXCEPTIOMS
# =======================================================================================

class EmptyConsensus(Exception):
    """Exception raised when the Vsearch consensus sequence is only made up of indels"""

# =======================================================================================
#               FUNCTIONS
# =======================================================================================

def get_vsearch_representative_consensus(vsearch_file: str) -> list[SeqRecord]:
    """
    Retrieve the longest sequence for each set of sequences forming the consensus.
    Parse the fasta file obtained through the '--msaout' flag.
    :param vsearch_file:
    :return:
    """
    records = SeqIO.parse(vsearch_file, "fasta")
    new_records = []
    idx = 0
    for record in records:
        if record.id.startswith("*"):
            record.id = record.id.replace("*", "").rsplit(" ", 1)[0] + f"_{idx}"
            record.name = record.name.replace("*", "").rsplit(" ", 1)[0] + f"_{idx}"
            record.description = ""
            new_records.append(record)
            idx += 1
    return new_records


def get_vsearch_regular_consensus(vsearch_file: str) -> list[SeqRecord]:
    """
    Retrieve the majority rule consensus obtained from vsearch.
    From the vsearch manual: 'the consensus sequence is
        constructed by taking the majority symbol (nucleotide or gap) from each column of the
        alignment. Columns containing a majority of gaps are skipped, except for terminal
        gaps'
    Parse the fasta file obtained through the '--msaout' flag.
    :param vsearch_file:
    :return:
    """
    records = list(SeqIO.parse(vsearch_file, "fasta"))
    prefix = records[0].id.rsplit("_", 1)[0].replace("*", "")
    new_records = []
    idx = 0
    for record in records:
        if record.id.startswith("consensus"):
            record.id = f"{prefix}_{idx}"
            record.name = f"{prefix}_{idx}"
            record.description = ""
            new_records.append(record)
            idx += 1

    return new_records


def _generate_kmer_consensus(
    msa: MultipleSeqAlignment, cluster_nbr: int
) -> Optional[SeqRecord]:
    """
    From a given cluster generate a kmer weighted majority rule consensus.
    :param msa: The multi-species alignment produced by vsearch.
    :param cluster_nbr: Index of the cluster, starting at 0.
    :return: None if there is a majority of gaps at all positions, raise EmptyConsensus exception.
    """
    prefix = msa[0].id.rsplit("_", 1)[0].replace("*", "")
    #  Case a single sequence with the consensus.
    if len(msa) == 2:
        record = msa[0]
        record.id = f"{prefix}_{cluster_nbr}"
        record.name = f"{prefix}_{cluster_nbr}"
        return record
    #  Case at least two sequences with the consensus.
    consensus = []
    kmers = [
        int(record.id.rsplit("-")[-1]) / len(str(record.seq).replace("-", ""))
        for record in msa[:-1]
    ]
    vsearch_consensus = [record for record in msa if record.id.startswith("consensus")][
        0
    ]
    for idx, nucl in enumerate(list(vsearch_consensus.seq)):
        if nucl != "-":
            vsearch_consensus_start = idx
            break
    else:
        raise EmptyConsensus(
            f"[ERROR] Cluster {prefix} produced an empty consensus sequence"
        )
        return

    for idx, nucl in enumerate(reversed(list(vsearch_consensus.seq))):
        if nucl != "-":
            vsearch_consensus_end = len(vsearch_consensus.seq) - idx
            break

    for idx in range(msa.get_alignment_length()):
        column = list(msa[:, idx].upper())
        seq_nucl = column[:-1]
        pairs = zip(seq_nucl, kmers)
        weights = defaultdict(int)
        for pair in pairs:
            weights[pair[0]] += pair[1]
        if idx < vsearch_consensus_start or vsearch_consensus_end <= idx:
            weights.pop("-")

        consensus.append(Counter(weights).most_common(1)[0][0])

    return SeqRecord(
        Seq("".join(consensus)),
        id=f"{prefix}_{cluster_nbr}",
        name=f"{prefix}_{cluster_nbr}",
        description="",
    )


def get_vsearch_kmer_consensus(vsearch_file: str) -> Optional[list[SeqRecord]]:
    """
    Generate a consensus sequence from a Vsearch cluster by weighting each position by the Kmer count
    and creating a column wise majority consensus. Columns containing a majority of gaps are skipped, except for
    initial and terminal gaps.
    :param vsearch_file:
    :return:
    """
    clusters_with_consensus = []
    records = list(SeqIO.parse(vsearch_file, "fasta"))
    pile = records[1:]
    cluster = [records[0]]
    while pile:
        target = pile[0]
        pile = pile[1:]
        if not pile:
            break
        if target.id.startswith("consensus"):
            cluster.append(target)
            clusters_with_consensus.append(cluster)
            cluster = []
        else:
            cluster.append(target)

    final_consensus = []
    for idx, cluster in enumerate(clusters_with_consensus):
        msa = MultipleSeqAlignment(cluster)
        final_consensus.append(_generate_kmer_consensus(msa, idx))
    return final_consensus

def merge_intervals(arr: list['OverlappingSeqs']):
    """
    Combine overlapping intervals of sequences and cluster their names.
    """
    # Sorting based on the increasing order of the start intervals
    arr.sort(key=attrgetter("probe_start"))
    index = 0
    index_to_del = []
    for i in range(1, len(arr)):
        #print(arr)
        if arr[index].probe_end >= arr[i].probe_start:
            arr[index].probe_end = max(arr[index].probe_end, arr[i].probe_end)
            arr[index].combine(arr[i].seq_names)
            index_to_del.append(i)
        else:
            index = index + 1
            arr[index] = arr[i]
    return [arr[index] for index in range(len(arr)) if index not in index_to_del]


# =======================================================================================
#               CLASSES
# =======================================================================================

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






def main():
    # tempfile.tempdir = "/temp"
    probe_fasta = "/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_paralogy/probe_105_aa.fasta"
    contig_fasta = "/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_paralogy/substample_0_no_gap.fasta"
    # probe_fasta = "/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_paralogy/probe_10052_aa.fasta"
    # contig_fasta = "/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_paralogy/final_alns/probe_10052.fasta"
    contig_fasta = "/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_paralogy/batarchium/H1_G3/H1_G3_probe_10248.fasta"
    probe_fasta = "/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_paralogy/probe_10248_aa.fasta"
    import pickle

    db = OverlapDetect(probe_fasta, contig_fasta)
    db.save_clusters("/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_paralogy/tmp")
    # with open(
    #     "/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_paralogy/probe_10052_aa.pkl",
    #     "wb",
    # ) as f:
    #     pickle.dump(db, f)
    # db = pickle.load(
    #     open(
    #         "/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_paralogy/probe_10052_aa.pkl",
    #         "rb",
    #     )
    # )
    # sorting_overlapping(db)

    saute_selected = "/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_paralogy/target_assembly/H1_A1/target_vars.fasta"
    saute_all = "/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_paralogy/target_assembly/H1_A1/all_vars.fasta"
    saute_graph = "/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_paralogy/target_assembly/H1_A1/graph.gfa"
    vsearch_file = "/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_paralogy/vsearch_aln.fasta"

    # records = get_vsearch_representative(vsearch_file)
    # records = get_vsearch_dumb_consensus(vsearch_file)
    # records = get_vsearch_kmer_consensus(vsearch_file)
    # SeqIO.write(
    #     records,
    #     "/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_paralogy/vsearch_kmer_consensus.fasta",
    #     "fasta",
    # )


if __name__ == "__main__":
    os.chdir("/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_clustering")
    # fragments = CdsParser(
    #     run_miniprot("probe_8631_aa.fasta", "gene_8631_node3.fasta")
    # ).get_fragments()
    # print(fragments)
    main()