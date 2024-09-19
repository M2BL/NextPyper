#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2024
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.
"""
Functions used to handle create protein consensuses out of amino acid probe sequences
that have been clustered with mmseqs2. Clusters are parse and their sequences aligned with
abPOA (https://github.com/yangao07/abPOA/tree/main/python) and their consensus is called
with the 'heaviest bunlding' algorithm. The input of generate_consensuses() is
the 'all_seqs.fasta' from mmseqs2.
#  Usage example:
    bio_seqrecords = generate_consensuses("all_seqs.fasta")
"""

__version__ = "0.1"

# =======================================================================================
#               IMPORTS
# =======================================================================================

from pathlib import Path
from typing import Optional

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq

import pyabpoa as pa


# =======================================================================================
#               EXCEPTIOMS
# =======================================================================================


# =======================================================================================
#               FUNCTIONS
# =======================================================================================
def _parse_all_seqs(fasta: str) -> list[list[SeqRecord]]:
    clusters = []
    records = list(SeqIO.parse(Path(fasta), "fasta"))
    cluster = [records[1]]
    pile = records[2:]
    while pile:
        target = pile[0]
        pile = pile[1:]
        if str(target.seq) == "":
            clusters.append(cluster)
            cluster = []
        else:
            cluster.append(target)
        if not pile:
            clusters.append(cluster)
            break
    return clusters


def _align_abpoa(cluster: list[SeqRecord], rec_id: str) -> SeqRecord:
    a = pa.msa_aligner(is_aa=True)
    seqs = [str(rec.seq) for rec in cluster]
    res = a.msa(seqs, out_cons=True, out_msa=False)
    if (sequence := res.cons_seq[0]) != "":
        return SeqRecord(Seq(sequence), id=rec_id, name="", description="")


def generate_consensuses(fasta: str) -> Optional[list[SeqRecord]]:
    mmseqs_clusters = _parse_all_seqs(fasta)
    consensuses = []
    if not mmseqs_clusters:
        return
    for idx, cluster in enumerate(mmseqs_clusters):
        rec_id = f"{cluster[0].id}_{idx}"
        if len(cluster) == 1:
            new_rec = cluster[0]
            new_rec.id = rec_id
            new_rec.name = ""
            new_rec.description = ""
            consensuses.append(new_rec)
            continue
        consensuses.append(_align_abpoa(cluster, rec_id))
    return consensuses


if __name__ == "__main__":
    mmseqs_fasta = "/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_clustering_final/all_proteins_cluster_all_seqs.fasta"
    consensuses = generate_consensuses(mmseqs_fasta)
    out = "/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_clustering_final/tmp/consensuses.fasta"
    SeqIO.write(consensuses, out, "fasta")
