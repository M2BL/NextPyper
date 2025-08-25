#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2024
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.

"""
This module filter sequences based on:

- GC content.
- Depth of coverage.
- MMseqs2 alignment hits against probes.

Returning the surviving sequences, and a table of metrics of the best constructed alignment
of each sequence with its corresponding probe. The table includes the identity and coverage
of the sequence to its probe.

Each filter has its own parameters described below.

1. GC filtering: removes sequences whose GC content falls outside the the range [min_gc, max_gc].
2. Coverage filtering: Remove sequences that do not meet depth of coverage criteria as
   estimated by mapping. The filtering can be performed in three modes (absolute, relative,
   and dynamic) and it is controlled by two parameters (cov_threshold and cov_dynamic_filt)

    - In absolute mode (cov_dynamic_filt=False, cov_threshold > 1): All sequences with a depth
      below the given threshold (e.g. 10) are filtered.
    - In relative mode (cov_dynamic_filt=False, cov_threshold in [0 - 1]): The threshold indicates
      the fraction of sequences to be removed. For instance, if cov_threshold=0.05, the 5% of the
      sequences with the lowest coverage will be filtered.
    - In dynamic mode (cov_dynamic_filt=True): A hard_threhold is computed from the distribution
      of depths by the formula: hard_threhold = median(depths) - cov_threshold * std(depths).
      The cov_threshold parameter modulates how strict the filter is from deviations from the
      observed median depth. To make the computation more robust, the 5% lowest depths and 5%
      highest depths are filtered from the distribution before the computation. Sequences with
      a depth below hard_threshold are filtered.

3. Alignment againt probes filtering: Filter the sequences whose best alignment against probes
   do not meet the minimum identity (scf_min_idt) or coverage (scf_min_cov) required.
"""

__version__ = "0.1"

# =======================================================================================
#               IMPORTS
# =======================================================================================

from pathlib import Path
import sys

import polars as pl
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction
from intervaltree import IntervalTree

from graph_utils import build_probe_trees, filt_probe_hits


# =======================================================================================
#               CONSTANTS
# =======================================================================================

REF_PAT = r".*probe([0-9]+)$"
SAUTE_PAT = r"^Contig_(?P<sample>.*?)-(?P<probe>.*?)_(?P<cluster>.+?)_(?P<seed>\d+?):(?P<component>\d+?):[^ ]+$"
METABAT_COLS = ["query", "len", "cov", "nu", "var"]

# =============================================================================
#                FUNCTIONS
# =============================================================================


def orient_scf(rec: SeqRecord, trans: bool) -> SeqRecord:
    if trans:
        rec.seq = rec.seq.reverse_complement()
    return rec


def tag_probe(rec: SeqRecord, probe: str) -> SeqRecord:
    """Add the probe information in the record name after the '-'"""

    idx = rec.id.index("-")
    rec.id = f"{rec.id[:idx+1]}{probe}_{rec.id[idx+1:]}"
    rec.name = rec.description = ""
    return rec


def filter_gc_recs(
    df: pl.DataFrame,
    recs_dict: dict[str, SeqRecord],
    max_gc: float = 0.6,
    min_gc: float = 0.2,
):
    """Given a dataframe of hits, filter the query sequences (scfs) with a higher
    gc than max_gc. Return the resulting dataframe."""

    high_gc_map = {
        name: min_gc <= gc_fraction(rec) <= max_gc for name, rec in recs_dict.items()
    }
    return df.filter(pl.col("query").replace_strict(high_gc_map))


def filter_by_cov(
    df: pl.DataFrame, cov_table: Path, threshold: int | float, dynamic: bool = True
):
    """Given a dataframe of hits, filter the query sequences (scfs) by the observed coverage
    from mapped reads."""

    meta_df = pl.read_csv(cov_table, separator="\t", new_columns=METABAT_COLS)
    covs = meta_df["cov"].to_numpy()

    if dynamic:
        low, high = np.quantile(covs, [0.05, 0.95])
        mask = (covs > low) & (covs < high)
        subarr = covs[mask]
        median = np.median(subarr)
        std = subarr.std()
        print(f"Coverage mode: dynamic ({median=:.2f}, {std=:.2f}, {threshold=})")
        limit_cov = np.median(subarr) - threshold * subarr.std()
        limit_cov = limit_cov if limit_cov > 0 else 0
    elif threshold >= 1:
        print(f"Coverage mode: absolute ({threshold=})")
        limit_cov = threshold
    elif threshold > 0:
        print(f"Coverage mode: relative ({threshold=})")
        limit_cov = np.quantile(covs, threshold)
    else:
        raise ValueError(
            f"""Error: {threshold=} has to be [0-1] for relative filtering
            or a positive integer for absolute filtering."""
        )
    print(f"Coverage threshold set at {limit_cov:.2f}")

    surviving_scfs = meta_df.filter(pl.col("cov") >= limit_cov).select("query")
    return df.join(surviving_scfs, on="query")


def compute_hits(
    df: pl.DataFrame,
    min_cov: float,
    min_idt: float,
    qpat: str = SAUTE_PAT,
    tpat: str = REF_PAT,
) -> pl.DataFrame:
    """Parse the alignments and determine the sequences that satisfy the min_cov
    and min_idt thresholds.

    If qpat is given, a pattern is expected to extract the probe version from the
    query name. If qpat is None, a best probe version is determined from the hits
    and picked. This probe information will be used to filter hits different to
    that probe.

    The remaining sequences (queries) are filtered to comply with the minimum
    coverage and identity to the best they match.
    """

    def _real_cov(tstart: list[int], tend: list[int]) -> int:
        tree = IntervalTree.from_tuples(zip(tstart, tend))
        tree.merge_overlaps(strict=False)
        return sum(inter.length() for inter in tree)

    df.replace_column(8, df["theader"].str.split(" ").list.first())

    # If query comes with probe information, filter to keep only matches of that probe.
    if qpat:
        pre_df = df.with_columns(
            pl.col("tstart") - 1,
            qprobe=pl.col("query").str.extract(qpat, 2),
            tprobe=pl.col("theader").str.extract(tpat, 1),
            cis=pl.col("qend") > pl.col("qstart"),
        ).filter(pl.col("qprobe") == pl.col("tprobe"))

    # Otherwise Determine the best probe.
    else:
        pre_df = df.with_columns(
            pl.col("tstart") - 1,
            tprobe=pl.col("theader").str.extract(tpat, 1),
            cis=pl.col("qend") > pl.col("qstart"),
        )

        probe_trees = build_probe_trees(pre_df, min_idt)
        pre_df = filt_probe_hits(pre_df, probe_trees)

    gdf = (
        pre_df.group_by(["query", "theader", "cis"])
        .agg(
            pl.sum("nident"),
            pl.sum("mismatch"),
            pl.sum("gapopen"),
            pl.first("tlen"),
            pl.first("tprobe"),
            pl.col("tstart"),
            pl.col("tend"),
            pl.first("glob_eff_cov"),
            pl.first("glob_cov"),
        )
        .with_columns(
            adj_cov=(
                pl.struct(["tstart", "tend"]).map_elements(
                    lambda x: _real_cov(**x), return_dtype=pl.Int64
                )
            )
        )
    )

    final_df = (
        gdf.with_columns(
            cov=pl.col("adj_cov") / pl.col("tlen"),
            idt=pl.col("nident") / (pl.col("nident") + pl.col("mismatch")),
        )
        .with_columns(eff_cov=(pl.col("cov") * pl.col("idt")))
        .filter((pl.col("cov") > min_cov) & (pl.col("idt") > min_idt))
        .group_by("query")
        .agg(pl.all().sort_by("eff_cov", descending=True).first())
    )

    return final_df.drop(["tstart", "tend"])


def filt_records(recs: list[SeqRecord], filt_ids: dict[str, bool]) -> list[SeqRecord]:
    return (
        orient_scf(rec, trans)
        for rec in recs
        if (trans := filt_ids.get(rec.id)) is not None
    )


def match_mmseqs_recs(
    rec_path: Path,
    table_path: Path,
    cov_path: Path,
    out_path: Path,
    out_table: Path,
    min_cov: float,
    min_idt: float,
    qpat: str,
    tpat: str,
    max_gc: float = 0.6,
    min_gc: float = 0.2,
    cov_threshold: int | float = 0.025,
    cov_dynamic_filt: bool = False,
    sep_probes: bool = False,
    tag_scfs: bool = False,
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

    recs_dict = SeqIO.to_dict(SeqIO.parse(rec_path, "fasta"))
    df = pl.read_csv(table_path, separator="\t", has_header=True)
    n0 = df["query"].unique().count()
    print(f"Starting seeds: {n0}")

    # Apply GC filter:
    df = filter_gc_recs(df, recs_dict, max_gc, min_gc)
    n1 = df["query"].unique().count()
    print(f"Seeds after GC filter: {n1} (removed {n0 - n1})")

    # Apply Coverage filter
    df = filter_by_cov(df, cov_path, cov_threshold, cov_dynamic_filt)
    n2 = df["query"].unique().count()
    print(f"Seeds after Depth filter: {n2} (removed {n1 - n2})")

    # Compute the hits from the filtered sequences
    filt_df = compute_hits(df, min_cov, min_idt, qpat, tpat)
    n3 = filt_df["query"].unique().count()
    print(f"Final set of seeds: {n3} (removed {n2 - n3})")

    # Save the computation results in the log
    if out_table:
        filt_df.write_csv(out_table, separator="\t")

    # Separate the surviving sequences by probe
    if sep_probes:
        out_path.mkdir(exist_ok=True)
        iter_df = filt_df.group_by("tprobe").all().select(["tprobe", "query", "cis"])

        for probe, ids, orient_list in iter_df.iter_rows():
            probe_recs = [
                tag_probe(orient_scf(recs_dict[rec_id], not cis), probe)
                for rec_id, cis in zip(ids, orient_list)
            ]
            SeqIO.write(probe_recs, out_path / f"{probe}.fasta", "fasta")

    # Do not separate the sequences in multiple files, but tag them.
    elif tag_scfs:
        iter_df = filt_df.select(["tprobe", "query", "cis"]).iter_rows()
        filt_scfs = [
            tag_probe(orient_scf(recs_dict[rec_id], not cis), probe)
            for probe, rec_id, cis in iter_df
        ]
        SeqIO.write(filt_scfs, out_path, "fasta")

    # Output all sequences in a single file (they are already tagged)
    else:
        filt_ids = dict(
            filt_df.select(["query", "cis"]).with_columns(~pl.col("cis")).iter_rows()
        )
        filt_scfs = filt_records(recs_dict.values(), filt_ids)
        SeqIO.write(filt_scfs, out_path, "fasta")


def snakemake_call(snakemake):
    with open(snakemake.log[0], "w") as outlog:
        sys.stdout = sys.stderr = outlog

        # Read input and outputs
        recs = Path(snakemake.input.scfs)
        hits = Path(snakemake.input.hits)
        covs = Path(snakemake.input.covs)
        seeds = Path(snakemake.output.scfs)
        out_table = Path(snakemake.output.metrics)

        # Params: GC based filtering
        max_gc = snakemake.params.max_gc
        min_gc = snakemake.params.min_gc

        # Params: Coverage filtering
        cov_threshold = snakemake.params.cov_threshold
        cov_dynamic_filt = snakemake.params.cov_dynamic_filt

        # Params: probe hits filtering
        min_idt = snakemake.params.min_idt
        min_cov = snakemake.params.min_cov
        qpat = snakemake.params.get("qpat", SAUTE_PAT)
        tpat = snakemake.params.get("tpat", REF_PAT)

        # Output params
        sep_probes = snakemake.params.separate_probes
        tag_scfs = snakemake.params.get("tag_scfs", False)

        match_mmseqs_recs(
            rec_path=recs,
            table_path=hits,
            cov_path=covs,
            out_path=seeds,
            out_table=out_table,
            min_cov=min_cov,
            min_idt=min_idt,
            qpat=qpat,
            tpat=tpat,
            max_gc=max_gc,
            min_gc=min_gc,
            cov_threshold=cov_threshold,
            cov_dynamic_filt=cov_dynamic_filt,
            sep_probes=sep_probes,
            tag_scfs=tag_scfs,
        )


def main():

    if len(sys.argv) != 4:
        print(
            "Usage: python homolog_filtering.py <scfs.fasta> <hits.tsv> <output.fasta>"
        )
        sys.exit(1)

    min_idt = 0.6
    min_cov = 0.0
    sep_probes = False

    recs = Path(sys.argv[1])
    table = Path(sys.argv[2])
    out = Path(sys.argv[3])

    qpat = None
    tpat = r"(\d+)$"

    match_mmseqs_recs(recs, table, out, min_cov, min_idt, qpat, tpat, sep_probes)


if __name__ == "__main__":
    if "snakemake" in globals():
        snakemake_call(snakemake)
    else:
        main()
