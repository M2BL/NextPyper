#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2024
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.

"""
Collect the potentially informative seeds for saute assembly based on vsearch clustering.
"""

__version__ = "0.1"

# =======================================================================================
#               IMPORTS
# =======================================================================================

from pathlib import Path
from collections import defaultdict
import re

from Bio import SeqIO
import numpy as np
import polars as pl

# =======================================================================================
#               CONSTANTS
# =======================================================================================

CLUSTER_COLS = [
    "type",
    "cluster_id",
    "size",
    "idt",
    "strand",
    "NU1",
    "NU2",
    "cigar",
    "query",
    "centroid",
]
REC_PAT = r"^(.*?)-.*_cov_([\d\.]+)$"

# =============================================================================
#                FUNCTIONS
# =============================================================================


def snakemake_call(snakemake):

    cluster_tables_dir = snakemake.input.cluster_tables
    sample_probes_dir = snakemake.input.samples
    seeds_out = snakemake.output.seeds
    cov_log_dir = Path(snakemake.log[0]).parent

    # Parse the records fo all the samples
    sample_recs = {
        sample.stem: {
            probe.stem: SeqIO.to_dict(SeqIO.parse(probe, "fasta"))
            for probe in Path(sample).glob("*.fasta")
        }
        for sample in map(Path, sample_probes_dir)
    }

    # Redistribute the seeds accroding to the clusters
    pat = re.compile(REC_PAT)
    sample_seeds = defaultdict(list)
    for probe in map(Path, cluster_tables_dir):
        if probe.stat().st_size == 0:
            continue

        table = pl.read_csv(
            probe, separator="\t", has_header=False, new_columns=CLUSTER_COLS
        ).with_columns(sample=pl.col("query").str.extract(REC_PAT))

        for sample in sample_recs:
            # Find the clusters where there are sequences of the sample
            clusters = (
                table.filter(pl.col("sample") == sample).select("cluster_id").unique()
            )

            # First, add the seeds from the current sample
            sample_seeds[sample].extend(
                sample_recs[sample].get(probe.stem, {}).values()
            )
            # Compute median coverage observed for sample sequences
            med_cov = np.median(
                [float(pat.search(rec.id)[2]) for rec in sample_seeds[sample]]
            )
            (cov_log_dir / f"{sample}.cov").write_text(f"{med_cov}")

            # Include all the sequences of the identified clusters as seeds for the sample
            for inter_sample, rec in (
                table.join(clusters, on="cluster_id", how="inner")
                .filter((pl.col("type") == "C") & (pl.col("sample") != sample))
                .select(["sample", "query"])
                .iter_rows()
            ):
                sample_seeds[sample].append(sample_recs[inter_sample][probe.stem][rec])

    # Write the seeds for all the samples
    seeds_out_dir = Path(seeds_out[0]).parent
    seeds_out_dir.mkdir(exist_ok=True, parents=True)

    for sample, recs in sample_seeds.items():
        SeqIO.write(recs, seeds_out_dir / f"{sample}.fasta", "fasta")


def main(): ...


if __name__ == "__main__":
    if "snakemake" in globals():
        snakemake_call(snakemake)
    else:
        main()
