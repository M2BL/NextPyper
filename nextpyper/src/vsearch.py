#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2024
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.
"""
Functions used to handle vsearch consensus.
Three consensus methods are available:
    get_vsearch_representative_consensus(): gets the centroid sequence of the cluster, i.e. sequence prefixed with '*'.
    get_vsearch_regular_consensus(): gets the vsearch native consensus that has problems with flanking missing sequences.
    get_vsearch_kmer_consensus(): gets a kmer count weighted consensus without using flanking missing sequences

"""

__version__ = "0.1"

# =======================================================================================
#               IMPORTS
# =======================================================================================

from collections import defaultdict, Counter, namedtuple
import os
from pathlib import Path
from typing import Optional, Literal

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


def get_vsearch_representative_consensus(vsearch_file: Path) -> list[SeqRecord]:
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


def get_vsearch_regular_consensus(vsearch_file: Path) -> list[SeqRecord]:
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
    msa: MultipleSeqAlignment,
    cluster_nbr: int,
    assembly_method: Literal["SAUTE", "SPAdes"],
) -> Optional[SeqRecord]:
    """
    From a given cluster generate a kmer weighted majority rule consensus.
    :param msa: The multi-species alignment produced by vsearch.
    :param cluster_nbr: Index of the cluster, starting at 0.
    :param assembly_method: For SAUTE the record.id ends with the kmer count e.g. probe_10010_0:1:1:25095,
            whereas currently the count in SPAdes is Hedypnois_rhagadioloides_ERR7618432-FUPX_6128_c6_p816
    :return: None if there is a majority of gaps at all positions, raise EmptyConsensus exception.
    """
    if assembly_method == "SAUTE":
        prefix = msa[0].id.rsplit("_", 1)[0].replace("*", "")
    else:
        # To be modified
        prefix = msa[0].id.rsplit("_", 1)[0].replace("*", "")
    print(f"{prefix=}")
    #  Case a single sequence with the consensus.
    if len(msa) == 2:
        record = msa[0]
        record.seq = Seq(str(record.seq).upper())
        record.id = f"{prefix}_{cluster_nbr}"
        record.name = ""
        record.description = ""
        return record

    #  Case at least two sequences with the consensus.
    consensus = []
    if assembly_method == "SAUTE":
        # Extract for each sequence the k-mer count, normalize by length
        kmers = [
            int(record.id.rsplit(":")[-1]) / len(str(record.seq).replace("-", ""))
            for record in msa[:-1]
        ]
    else:
        # Extract for each sequence the k-mer count, do not normalize by length as for SPAdes the
        # normalization has been done during the scaffold creation from the contigs.
        kmers = [
            float(record.id.split("_")[-1].removeprefix("DP")) for record in msa[:-1]
        ]

    SequenceBoundaries = namedtuple("SequenceBoundaries", ["start", "end"])

    # for each sequence find the start and end of the sequence, remove its weight if the idx is in the flanking region
    def find_start_end(seq: Seq) -> SequenceBoundaries:
        sequence = list(str(seq))
        if set(sequence) == {"-"}:
            return SequenceBoundaries(None, None)
        length = len(sequence)
        revseq = sequence[::-1]
        return SequenceBoundaries(
            next(
                i
                for i, nucleotide in enumerate(sequence)
                if nucleotide in ["A", "T", "C", "G"]
            ),
            length
            - next(
                i
                for i, nucleotide in enumerate(revseq)
                if nucleotide in ["A", "T", "C", "G"]
            ),
        )

    #  For each sequence find where in the alignment, the nucleotides start and end.
    #  Store in dict with keys the index of the sequence in the msa.
    sequences_boundaries = {}
    for idx, record in enumerate(msa[:-1]):
        sequences_boundaries[idx] = find_start_end(record.seq)

    for idx in range(msa.get_alignment_length()):
        column = list(msa[:, idx].upper())
        seq_nucl = column[:-1]
        pairs = zip(seq_nucl, kmers)
        weights = defaultdict(int)
        for rec_index, pair in enumerate(pairs):
            boundaries = sequences_boundaries[rec_index]
            #  Ignore indels from flanking regions
            if boundaries.start is None or boundaries.end is None:
                continue
            if idx < boundaries.start or idx >= boundaries.end:
                continue
            weights[pair[0]] += pair[1]
        if weights:
            consensus.append(Counter(weights).most_common(1)[0][0])

    if set(consensus) == {"-"} or not consensus:
        raise EmptyConsensus

    return SeqRecord(
        Seq("".join(consensus).replace("-", "")),
        id=f"{prefix}_{cluster_nbr}",
        name="",
        description="",
    )


def get_vsearch_kmer_consensus(
    vsearch_file: Path, assembly_method: Literal["SAUTE", "SPAdes"]
) -> list[SeqRecord]:
    """
    Generate a consensus sequence from a Vsearch cluster by weighting each position by the Kmer count
    and creating a column wise majority consensus. Columns containing a majority of gaps are skipped, except for
    initial and terminal gaps. Returns an empty list if no consensus was generated.
    :param vsearch_file:
    :param assembly_method: For SAUTE the record.id ends with the kmer count e.g. probe_10010_0:1:1:25095,
            whereas currently the count in SPAdes is Hedypnois_rhagadioloides_ERR7618432-FUPX_6128_c6_p816
    :return:
    """
    final_consensuses = []
    records = list(SeqIO.parse(vsearch_file, "fasta"))
    pile = records[1:]
    cluster = [records[0]]
    idx = 0
    while pile:
        target = pile[0]
        pile = pile[1:]
        if not pile:
            break
        if target.id.startswith("consensus"):
            cluster.append(target)
            msa = MultipleSeqAlignment(cluster)
            try:
                final_consensuses.append(
                    _generate_kmer_consensus(msa, idx, assembly_method)
                )
                cluster = []
            except EmptyConsensus:
                pass
            finally:
                idx += 1
        else:
            cluster.append(target)

    return final_consensuses


if __name__ == "__main__":
    os.chdir(
        "/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_clustering_final"
    )
    # fragments = CdsParser(
    #     run_miniprot("probe_8631_aa.fasta", "gene_8631_node3.fasta")
    # ).get_fragments()
    # print(fragments)
    vsearch_result = Path(
        "/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_clustering_final/vsearch_6487_aln.fasta"
    )
    out = "/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_clustering_final/vsearch_6487_con_test.fasta"
    vsearch_result = Path(
        "/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_clustering_final/empty.fasta"
    )
    out = "/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_clustering_final/empty_con.fasta"
    records_con = get_vsearch_kmer_consensus(Path(vsearch_result), "SPAdes")
    SeqIO.write(records_con, out, "fasta")
