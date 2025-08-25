#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2024
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.


"""Parse saute assembly and do multiple operations on them, including collapsing the alleles,
optionally limiting the number of variants per component too, and Splitting into well behaved
and "explosive" components."""

__version__ = "0.1"

# =======================================================================================
#               IMPORTS
# =======================================================================================

from pathlib import Path
import sys

import polars as pl
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from more_itertools import partition


# =======================================================================================
#                CONSTANTS
# =======================================================================================

FIELDS = {
    "seed": str,
    "probe": int,
    "seed_id": int,
    "seed_len": int,
    "cov": float,
    "comp": int,
    "ctg1": int,
    "ctg2": int,
    "kmers": int,
}

# ALL_PAT = r"^(?P<sample>.*?)-(?P<probe>.*?)_EDGE_(?P<seed_id>\d+)_length_(?P<len>\d+)_cov_(?P<cov>[\w.]+):[^ ]+:(?P<kmers>\d+)$"
TARGET_PAT = r"^Contig_(?P<seed>.*?)-(?P<probe>.*?)_EDGE_(?P<seed_id>\d+)_length_(?P<seed_len>\d+)_cov_(?P<cov>[\w.]+):(?P<comp>\d+):(?P<ctg1>\d+):(?P<ctg2>\d+):(?P<kmers>\d+)$"


# =======================================================================================
#                FUNCTIONS
# =======================================================================================


def query2df(recs: list[SeqRecord], pat: str, schema: dict[str, str]) -> pl.DataFrame:
    recs_ids, lens = zip(*map(lambda rec: (rec.id, len(rec)), recs))
    return (
        pl.DataFrame(data={"query": recs_ids, "len": lens})
        .with_columns(pl.col("query").str.extract_groups(pat).struct.unnest())
        .cast(schema)
    )


def collapse_alleles(recs: list[SeqRecord], df: pl.DataFrame) -> list[SeqRecord]:
    """Group the sequences from the same component with the same length and pick the one
    with the highest kmer count. It is assumed that these sequences are different alleles.

    Return a reduced set of records with alleles collapsed."""

    target_vars = df.group_by("probe", "seed", "seed_id", "comp", "len").agg(
        pl.col("query").gather(pl.col("kmers").arg_max()).first(),
        pl.col("kmers").max(),
        count=pl.len(),
    )

    allele_collapsed_set = set(target_vars["query"])
    return [rec for rec in recs if rec.id in allele_collapsed_set]


def snakemake_call(snakemake):
    records_path = Path(snakemake.input[0])
    out_path = Path(snakemake.output.normal)
    # expl_out_path = Path(snakemake.output.expl)

    pattern = snakemake.params.get("pattern", TARGET_PAT)

    recs = list(SeqIO.parse(records_path, "fasta"))
    df = query2df(recs, pattern, FIELDS)
    new_recs = collapse_alleles(recs, df)

    SeqIO.write(new_recs, out_path, "fasta")


def main():
    if len(sys.argv) != 3:
        print("Usage: python var_asm_parser.py <saute_asm.fasta> <output.fasta>")
        sys.exit(1)

    records_path = Path(sys.argv[1])
    out_path = Path(sys.argv[2])

    recs = list(SeqIO.parse(records_path, "fasta"))
    df = query2df(recs, TARGET_PAT, FIELDS)
    new_recs = collapse_alleles(recs, df)
    SeqIO.write(new_recs, out_path, "fasta")


if __name__ == "__main__":
    if "snakemake" in globals():
        snakemake_call(snakemake)
    else:
        main()
