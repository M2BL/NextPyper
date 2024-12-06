#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2024
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.

"""
Filter sequences based on mmseqs2 alignments against probes.
"""

__version__ = "0.1"

# =======================================================================================
#               IMPORTS
# =======================================================================================

from pathlib import Path
import polars as pl
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys


# =======================================================================================
#               CONSTANTS
# =======================================================================================

REF_PAT = r".*probe([0-9]+)$"
SAUTE_PAT = r"^Contig_(?P<sample>.*?)-(?P<probe>.*?)_(?P<cluster>\d+?)_(?P<seed>\d+?):(?P<component>\d+?):[^ ]+$"


# =============================================================================
#                FUNCTIONS
# =============================================================================


def orient_scf(rec: SeqRecord, trans: bool) -> SeqRecord:
    if trans:
        rec.seq = rec.seq.reverse_complement()
    return rec


def compute_hits(df: pl.DataFrame, min_cov: float, min_idt: float) -> pl.DataFrame:
    """Parse the alignments and determine the sequences that satisfy the min_cov
    and min_idt thresholds."""

    final_df = (
        df.with_columns(
            qprobe=pl.col("query").str.extract(SAUTE_PAT, 2).cast(pl.Int64),
            tprobe=pl.col("theader").str.extract(REF_PAT, 1).cast(pl.Int64),
            cis=pl.col("qend") > pl.col("qstart"),
        )
        .filter(pl.col("qprobe") == pl.col("tprobe"))
        .group_by(["query", "theader", "cis"])
        .agg(
            pl.sum("nident"),
            pl.sum("mismatch"),
            pl.sum("gapopen"),
            pl.first("qlen"),
        )
        .with_columns(
            cov=(pl.col("nident") + pl.col("mismatch")) * 3 / pl.col("qlen"),
            idt=pl.col("nident")
            / (pl.col("nident") + pl.col("mismatch") + pl.col("gapopen")),
        )
        .filter((pl.col("cov") > min_cov) & (pl.col("idt") > min_idt))
        .group_by("query")
        .agg(pl.all().sort_by("idt").last())
        .select(["query", "cis"])
        .with_columns(~pl.col("cis"))
    )

    return final_df


def filt_records(recs: list[SeqRecord], filt_ids: dict[str, bool]) -> list[SeqRecord]:
    return (
        orient_scf(rec, trans)
        for rec in recs
        if (trans := filt_ids.get(rec.id)) is not None
    )


def match_mmseqs_recs(
    rec_path: Path, table_path: Path, out_path: Path, min_cov: float, min_idt: float
) -> None:

    recs = list(SeqIO.parse(rec_path, "fasta"))
    df = pl.read_csv(table_path, separator="\t", has_header=True)
    filt_ids = dict(compute_hits(df, min_cov, min_idt).iter_rows())
    filt_scfs = filt_records(recs, filt_ids)
    SeqIO.write(filt_scfs, out_path, "fasta")


def snakemake_call(snakemake):
    with open(snakemake.log[0], "w") as outlog:
        sys.stdout = sys.stderr = outlog

        min_idt = snakemake.params.min_idt
        min_cov = snakemake.params.min_cov

        recs_path = Path(snakemake.input.scfs)
        table_path = Path(snakemake.input.table)
        out_recs = Path(snakemake.output[0])

        match_mmseqs_recs(recs_path, table_path, out_recs, min_idt, min_cov)


def main(): ...


if __name__ == "__main__":
    if "snakemake" in globals():
        snakemake_call(snakemake)
    else:
        main()
