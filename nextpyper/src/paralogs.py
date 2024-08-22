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



# =======================================================================================
#               CLASSES
# =======================================================================================
