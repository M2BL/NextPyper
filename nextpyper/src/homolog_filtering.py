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
SAUTE_PAT = r"^Contig_(?P<sample>.*?)-(?P<probe>.*?)_(?P<cluster>.+?)_(?P<seed>\d+?):(?P<component>\d+?):[^ ]+$"

# =============================================================================
#                FUNCTIONS
# =============================================================================


def orient_scf(rec: SeqRecord, trans: bool) -> SeqRecord:
    if trans:
        rec.seq = rec.seq.reverse_complement()
    return rec


def tag_probe(rec: SeqRecord, probe) -> SeqRecord:
    idx = rec.id.index("-")
    rec.id = f"{rec.id[:idx+1]}{probe}_{rec.id[idx+1:]}"
    return rec


def compute_hits(
    df: pl.DataFrame,
    min_cov: float,
    min_idt: float,
    qpat: str = SAUTE_PAT,
    tpat: str = REF_PAT,
) -> pl.DataFrame:
    """Parse the alignments and determine the sequences that satisfy the min_cov
    and min_idt thresholds.

    If qpat is None the sequences are not filtered by matching query and target probes.
    """

    df.replace_column(8, df["theader"].str.split(" ").list.first())
    if qpat:
        pre_df = df.with_columns(
            qprobe=pl.col("query").str.extract(qpat, 2),
            tprobe=pl.col("theader").str.extract(tpat, 1),
            cis=pl.col("qend") > pl.col("qstart"),
        ).filter(pl.col("qprobe") == pl.col("tprobe"))
    else:
        pre_df = df.with_columns(
            tprobe=pl.col("theader").str.extract(tpat, 1),
            cis=pl.col("qend") > pl.col("qstart"),
        )

    final_df = (
        pre_df.group_by(["query", "theader", "cis"])
        .agg(
            pl.sum("nident"),
            pl.sum("mismatch"),
            pl.sum("gapopen"),
            pl.first("tlen"),
            pl.first("tprobe"),
        )
        .with_columns(
            cov=(pl.col("nident") + pl.col("mismatch")) / pl.col("tlen"),
            idt=pl.col("nident")
            / (pl.col("nident") + pl.col("mismatch") + pl.col("gapopen")),
        )
        .filter((pl.col("cov") > min_cov) & (pl.col("idt") > min_idt))
        .group_by("query")
        .agg(pl.all().sort_by("idt").last())
    )

    return final_df


def filt_records(recs: list[SeqRecord], filt_ids: dict[str, bool]) -> list[SeqRecord]:
    return (
        orient_scf(rec, trans)
        for rec in recs
        if (trans := filt_ids.get(rec.id)) is not None
    )


def match_mmseqs_recs(
    rec_path: Path,
    table_path: Path,
    out_path: Path,
    min_cov: float,
    min_idt: float,
    qpat: str,
    tpat: str,
    sep_probes: bool = False,
    log_results: bool = True,
) -> None:
    """Given a set of file with sequences, and a table with mmseqs2 matches against a set
    of probes, filter the sequences to those that match the probes with at least a minimum
    coverage and a minimum identity.

    A target pattern (tpat) is required to extract the probe id from the probes db hits. Similarly a
    query pattern (qpat) can be provided, although is optional. If provided, the alignments will be
    filtered to those that match query_probe and target_probe.

    Finally, is sep_probes is true the sequences will be separated per probe in an output directory.
    Otherwise, all sequences will be written together in a single file.
    """

    recs = list(SeqIO.parse(rec_path, "fasta"))
    df = pl.read_csv(table_path, separator="\t", has_header=True)

    filt_df = compute_hits(df, min_cov, min_idt, qpat, tpat)

    # Save the computation results in the log
    if log_results:
        filt_df.write_csv(sys.stdout, separator="\t")

    # Separate the surviving sequences by probe
    if sep_probes:
        out_path.mkdir(exist_ok=True)
        recs_dict = {rec.id: rec for rec in recs}
        iter_df = filt_df.group_by("tprobe").all().select(["tprobe", "query", "cis"])

        for probe, ids, orient_list in iter_df.iter_rows():
            probe_recs = [
                tag_probe(orient_scf(recs_dict[rec_id], not cis), probe)
                for rec_id, cis in zip(ids, orient_list)
            ]
            SeqIO.write(probe_recs, out_path / f"{probe}.fasta", "fasta")

    # Output all sequences in a single file
    else:
        filt_ids = dict(
            filt_df.select(["query", "cis"]).with_columns(~pl.col("cis")).iter_rows()
        )
        filt_scfs = filt_records(recs, filt_ids)
        SeqIO.write(filt_scfs, out_path, "fasta")


def snakemake_call(snakemake):
    with open(snakemake.log[0], "w") as outlog:
        sys.stdout = sys.stderr = outlog

        min_idt = snakemake.params.min_idt
        min_cov = snakemake.params.min_cov
        sep_probes = snakemake.params.separate_probes

        recs = Path(snakemake.input.scfs)
        table = Path(snakemake.input.table)
        out = Path(snakemake.output[0])

        qpat = snakemake.params.get("qpat", SAUTE_PAT)
        tpat = snakemake.params.get("tpat", REF_PAT)

        match_mmseqs_recs(recs, table, out, min_cov, min_idt, qpat, tpat, sep_probes)


def main(): ...


# min_idt = 0.1
# min_cov = 0.1
# sep_probes = False

# recs = Path(sys.argv[1])
# table = Path(sys.argv[2])
# out = Path(sys.argv[3])

# qpat = SAUTE_PAT
# tpat = r"(\d+)$"

# match_mmseqs_recs(recs, table, out, min_cov, min_idt, qpat, tpat, sep_probes)


if __name__ == "__main__":
    if "snakemake" in globals():
        snakemake_call(snakemake)
    else:
        main()
