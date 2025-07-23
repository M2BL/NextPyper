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
from collections import Counter, defaultdict
import re

from Bio import SeqIO
import numpy as np
import polars as pl

from multi_seq_probes import group_probes

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
REC_PAT = r"^(.*?)-(.*?)_.*_cov_([\d\.]+)$"

# =============================================================================
#                FUNCTIONS
# =============================================================================


def snakemake_call(snakemake):

    cluster_tables_dir = snakemake.input.cluster_tables
    sample_probes_dir = snakemake.input.samples
    probes_path = snakemake.input.probes
    seeds_out = snakemake.output.seeds
    probes_pattern = snakemake.params.pattern
    multi_probes = snakemake.params.is_multi
    min_sister_freq = snakemake.params.min_sister_freq
    cov_log_dir = Path(snakemake.log[0]).parent

    # Parse the records fo all the samples
    sample_recs = {
        sample.stem: {
            probe.stem: SeqIO.to_dict(SeqIO.parse(probe, "fasta"))
            for probe in Path(sample).glob("*.fasta")
        }
        for sample in map(Path, sample_probes_dir)
    }

    # Read the probes in order to supplement the seeds in case they are missing
    probes = list(SeqIO.parse(probes_path, "fasta"))
    if multi_probes:
        probes_dict = group_probes(probes, probes_pattern)
    else:
        probe_pat = re.compile(probes_pattern)
        probes_dict = {probe_pat.search(probe.id)[1]: [probe] for probe in probes}

    # Reformat the probes name to follow the seeds pattern
    for probe, recs in probes_dict.items():
        for rec in recs:
            rec.id = f"{rec.id}-{probe}_EDGE_0_length_{len(rec)}_cov_0.0"
            rec.name = rec.description = ""

    pat = re.compile(REC_PAT)
    probe_tables = {}
    sister_samples = defaultdict(Counter)
    for probe in map(Path, cluster_tables_dir):
        if probe.stat().st_size == 0:
            continue

        probe_tables[probe.stem] = pl.read_csv(
            probe, separator="\t", has_header=False, new_columns=CLUSTER_COLS
        ).with_columns(sample=pl.col("query").str.extract(REC_PAT))

        for sample in sample_recs:
            clusters = (
                probe_tables[probe.stem]
                .filter(pl.col("sample") == sample)
                .select("cluster_id")
                .unique()
            )
            sister_samples[sample].update(
                probe_tables[probe.stem]
                .filter(pl.col("cluster_id").is_in(clusters))["sample"]
                .unique()
            )

    # Compute the sister samples based on the votes of each probe
    sister_samples = {
        sample: {
            sister
            for sister, count in counts.most_common()
            if sister != sample and (count / counts[sample]) > min_sister_freq
        }
        for sample, counts in sister_samples.items()
    }

    # Redistribute the seeds according complementing with the centroids of clusters
    # where sister samples are present
    # This is a clusters centric approach, where we take the seeds from the current sample
    # and then we add the centroids of the clusters its sister samples matched as extra seeds.
    sample_seeds = defaultdict(list)
    for probe, table in probe_tables.items():
        for sample in sample_recs:
            # First, add the seeds from the current sample
            sample_seeds[sample].extend(sample_recs[sample].get(probe, {}).values())

            # Compute median coverage observed for sample sequences
            med_cov = np.median(
                [float(pat.search(rec.id)[3]) for rec in sample_seeds[sample]]
            )
            (cov_log_dir / f"{sample}.cov").write_text(f"{med_cov}")

            # Find the clusters that sister samples matched
            clusters = (
                table.filter(pl.col("sample").is_in(sister_samples[sample]))
                .select("cluster_id")
                .unique()
            )

            # Then, add the centroids of the clusters as extra seeds
            for inter_sample, rec in (
                table.join(clusters, on="cluster_id", how="inner")
                .filter((pl.col("type") == "C") & (pl.col("sample") != sample))
                .select(["sample", "query"])
                .iter_rows()
            ):
                sample_seeds[sample].append(sample_recs[inter_sample][probe][rec])

            # Finally, add the probes that were not present in the sample seeds.
            covered_probes = {pat.search(seed.id)[2] for seed in sample_seeds[sample]}
            for missing_probe in set(probes_dict) - covered_probes:
                sample_seeds[sample].extend(probes_dict[missing_probe])

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
