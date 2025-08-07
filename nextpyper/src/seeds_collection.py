#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2024
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.

"""
This module does two tasks in preparation for assembly with saute:

1. Collects potentially informative seeds based on vsearch clustering.
2. A Heuristic computes the kmer sizes based on the depth observed during the seeds assembly.

The seeds for a given sample are composed of:
    - The seeds assembled for that sample (intra-seeds).
    - The seeds of "sister-samples" (inter-seeds).
    - The probes sequences for those probes that were not represented already
      in the intra or inter-seeds.

#ToDo: Expand the docstring to include more details about the sister samples and inter-seeds.

The computation of the kmer sizes is based on the depth observed during the seeds assembly.
The read median depth of the sample observed during assembly is computed from the kmer depth
and kmer size found during seeds assembly. The heuristic consists in computing the primary
and secondary kmer sizes by targeting a given primary and secondary depth. The kmer sizes
are constrained to a range of values defined in the heuristic parameters.
"""

__version__ = "0.1"

# =======================================================================================
#               IMPORTS
# =======================================================================================

from pathlib import Path
from typing import TypedDict
from collections import Counter, defaultdict
import json
import re

from Bio import SeqIO
import numpy as np
import polars as pl
from more_itertools import one

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
METABAT_COLS = ("query", "len", "cov", "nu", "var")

# =============================================================================
#                CLASSES
# =============================================================================


class KmerParams(TypedDict):
    """Parameters for kmer size computation."""

    primary_target_depth: int
    secondary_target_depth: int
    k1_min: float
    k1_max: float
    k2_min: float
    k2_max: float


# =============================================================================
#                FUNCTIONS
# =============================================================================


def get_spades_kmer_size(spades_folder: Path) -> int:
    """Get the kmer size used by spades."""
    kspades = one(spades_folder.glob("K*/final_contigs.paths")).parent.name[1:]
    return int(kspades)


def get_map_coverage(cov_table: Path) -> float:
    """Get the median depth of coverage found by read mapping."""
    meta_df = pl.read_csv(cov_table, separator="\t", new_columns=METABAT_COLS)
    return meta_df["cov"].median()


def get_read_length(read_stats: Path) -> int:
    """Get the read length from the fastp report."""

    summary = json.loads(read_stats.read_text())
    return int(summary["summary"]["after_filtering"]["read2_mean_length"])


def compute_saute_kmers(
    read_med_cov: float,
    L: int,
    kmer_params: KmerParams,
) -> str:
    """Compute the kmer sizes to be used by saute based on the median coverage,
    read length, and the target depths for primary and secondary kmers. The kmer
    sizes are constrained to a range of values defined in the heuristic parameters."""

    # Compute kmer sizes for the given targets depths
    k1 = int(L * (1 - (kmer_params["primary_target_depth"] / read_med_cov)) + 1)
    k2 = int(L * (1 - (kmer_params["secondary_target_depth"] / read_med_cov)) + 1)

    # Constrain the kmer sizes to "sane" ranges (defined in the parameters):
    k1_min, k1_max = int(L * kmer_params["k1_min"]), int(L * kmer_params["k1_max"])
    k2_min, k2_max = int(L * kmer_params["k2_min"]), int(L * kmer_params["k2_max"])

    # Cap primary and secondary kmers to their respective absolute floors (Def: 49, 21)
    k1_min = max(k1_min, int(kmer_params["k1_floor"]))
    k2_min = max(k2_min, int(kmer_params["k2_floor"]))

    if k2 < k2_min:
        k2 = k2_min
    elif k2 > k2_max:
        k2 = k2_max

    if k1 < k1_min:
        k1 = k1_min
    elif k1 > k1_max:
        k1 = k1_max

    # Ensure odd kmers
    k2 = k2 - 1 if k2 % 2 == 0 else k2
    k1 = k1 - 1 if k1 % 2 == 0 else k1

    # Pack the results to log them later
    kmer_results = {
        "read_med_cov": read_med_cov,
        "primary_target_depth": kmer_params["primary_target_depth"],
        "secondary_target_depth": kmer_params["secondary_target_depth"],
        "k1_min": k1_min,
        "k1_max": k1_max,
        "k2_min": k2_min,
        "k2_max": k2_max,
        "k1": k1,
        "k2": k2,
    }

    return kmer_results


def infer_sister_samples(
    probe_tables: dict[str, dict[str, pl.DataFrame]],
) -> dict[str, dict[str, tuple[int, int]]]:
    """Use the clustering results for all the probes to determine the sister samples
    of all samples in the run.

    For a given sample and probe, look at which other samples co-occur in the same
    clusters for that probe. Those will constitute the "sister-samples" for that sample
    according to that probe.

    After consulting all the probes, keep the sister samples that were voted by more than
    min_sister_freq proportion of all probes.

    Return a dictionary of samples as key and their inferred sister-samples as values.
    """

    sister_samples = defaultdict(Counter)
    samples = pl.concat(probe_tables.values())["sample"].unique().to_list()

    for probe, table in probe_tables.items():
        for sample in samples:
            clusters = set(
                table.filter(pl.col("sample") == sample)["cluster_id"].unique()
            )
            sister_samples[sample].update(
                table.filter(pl.col("cluster_id").is_in(clusters))["sample"].unique()
            )

    # Compute the sister samples based on the votes of each probe
    sister_samples = {
        sample: {
            sister: (count, counts[sample])
            for sister, count in counts.most_common()
            if sister != sample
        }
        for sample, counts in sister_samples.items()
    }

    return sister_samples


def snakemake_call(snakemake):

    # Read inputs for seed collection
    cluster_tables_dir = snakemake.input.cluster_tables
    sample_probes_dir = snakemake.input.samples
    probes_path = snakemake.input.probes

    # Read inputs for saute kmer size computation
    spades_folders = snakemake.input.spades_folders
    read_stats = snakemake.input.read_stats
    cov_paths = snakemake.input.covs

    # Read outputs
    seeds_out = snakemake.output.seeds
    saute_params = snakemake.output.saute_params

    # Params for seed collection
    probes_pattern = snakemake.params.pattern
    multi_probes = snakemake.params.is_multi
    interseeds_use = snakemake.params.interseeds_use
    min_sister_freq = snakemake.params.min_sister_freq

    # Params for saute kmer size computation
    cov_by_mapping = snakemake.params.cov_by_mapping
    heuristic_params = snakemake.params.heuristic_params

    # Parse the records fo all the samples
    sample_recs = {
        sample.stem: {
            probe: SeqIO.to_dict(recs)
            for probe, recs in group_probes(
                list(SeqIO.parse(sample, "fasta")), REC_PAT, 2
            ).items()
        }
        for sample in map(Path, sample_probes_dir)
    }

    ## Kmer size parameters computation
    # Organize the Spades folders and fastp reports to make them accessible by sample
    spades_folders = {folder.name: folder for folder in map(Path, spades_folders)}
    read_covs = {file.stem: get_map_coverage(file) for file in map(Path, cov_paths)}
    read_lengths = {file.stem: get_read_length(file) for file in map(Path, read_stats)}
    kmer_params_out = {file.stem: file for file in map(Path, saute_params)}

    # Compute the kmer sizes to be used by saute for each sample
    pat = re.compile(REC_PAT)
    for sample, probe_recs in sample_recs.items():
        L = read_lengths[sample]

        # Compute median depth observed in base space
        # By mapping
        if cov_by_mapping:
            read_med_cov = read_covs[sample]
            kmer_log = {"read_length": L, "cov_estimation": "mapping"}

        # By spades assembly
        else:
            kspades = get_spades_kmer_size(spades_folders[sample])
            med_cov = np.median(
                [
                    float(pat.search(rec)[3])
                    for recs in probe_recs.values()
                    for rec in recs
                ]
            )
            read_med_cov = med_cov * L / (L - kspades + 1)
            kmer_log = {
                "read_length": L,
                "cov_estimation": "assembly",
                "kspades": kspades,
            }

        kmer_results = compute_saute_kmers(read_med_cov, L, heuristic_params)
        kmer_log.update(kmer_results)

        # Write the kmer results to a log file
        kmer_params_out[sample].write_text(json.dumps(kmer_log, indent=4))

    ## Seed collection
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

    # Load the clustering tables
    probe_tables = {}
    for probe in map(Path, cluster_tables_dir):
        if probe.stat().st_size == 0:
            continue

        probe_tables[probe.stem] = pl.read_csv(
            probe, separator="\t", has_header=False, new_columns=CLUSTER_COLS
        ).with_columns(sample=pl.col("query").str.extract(REC_PAT))

    # Subset the inter-seeds with the sister-samples approach
    if interseeds_use == "sister":
        sister_dir = Path(saute_params[0]).parent.parent / "sister_samples"
        sister_dir.mkdir(parents=True, exist_ok=True)

        sister_samples_counts = infer_sister_samples(probe_tables)

        # Log the sister samples counts
        for sample, sister_samples in sister_samples_counts.items():
            (sister_dir / f"{sample}.json").write_text(json.dumps(sister_samples))

        # Filter to keep the sister samples above the frequency threshold
        sister_samples = {
            sample: {
                sister
                for sister, (count, max_count) in sister_samples.items()
                if count / max_count > min_sister_freq
            }
            for sample, sister_samples in sister_samples_counts.items()
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

            # If following sister-sample strategy, subset the table to sister clusters
            if interseeds_use == "sister":
                # Find the clusters that sister samples matched
                clusters = (
                    table.filter(pl.col("sample").is_in(sister_samples[sample]))
                    .select("cluster_id")
                    .unique()
                )
                table = table.join(clusters, on="cluster_id", how="inner")

            # Then, add the centroids of the clusters as extra seeds (inter-seeds)
            if interseeds_use != "none":
                for inter_sample, rec in (
                    table.filter((pl.col("type") == "C") & (pl.col("sample") != sample))
                    .select(["sample", "query"])
                    .iter_rows()
                ):
                    sample_seeds[sample].append(sample_recs[inter_sample][probe][rec])

    # Finally, add the probes that were not present in the sample seeds.
    for sample, seeds in sample_seeds.items():
        covered_probes = {pat.search(seed.id)[2] for seed in seeds}
        for missing_probe in set(probes_dict) - covered_probes:
            seeds.extend(probes_dict[missing_probe])

    # Write the seeds for all the samples
    seeds_out_dir = Path(seeds_out[0]).parent
    seeds_out_dir.mkdir(exist_ok=True, parents=True)

    for sample, seeds in sample_seeds.items():
        SeqIO.write(seeds, seeds_out_dir / f"{sample}.fasta", "fasta")


def main(): ...


if __name__ == "__main__":
    if "snakemake" in globals():
        snakemake_call(snakemake)
    else:
        main()
