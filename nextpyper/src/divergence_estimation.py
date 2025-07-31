#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2024
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.

"""Estimate the divergence of a given sample (relative) to the probes, using the identities
computed by homologs filtering. A Gaussian is fit to the distribution of the alignments after
filtering by a minimum identity and coverage. The estimated parameters (mu, sigma) of
each sample are reported.
"""

__version__ = "0.1"

# =======================================================================================
#               IMPORTS
# =======================================================================================

from pathlib import Path
import json

import numpy as np
import polars as pl
from polars import DataFrame
from scipy.optimize import curve_fit


# =============================================================================
#                FUNCTIONS
# =============================================================================


def gaussian(x, mu, sigma, a):
    """Defines a classical Gaussian distribution."""
    return a * np.exp(-((x - mu) ** 2) / 2 / sigma**2)


def get_histogram_data(df: DataFrame, mmseqs_idt_threshold: float) -> tuple[list, list]:
    """ "
    Automatically infer the best number of bins in the histogram.
    Return x, y coordinates with y the height of the bar and y the central value of the bar.
    """

    data = df.filter((pl.col("idt") >= mmseqs_idt_threshold) & (pl.col("idt") <= 1))[
        "idt"
    ].to_numpy()
    counts, bins = np.histogram(data, bins="fd", range=(mmseqs_idt_threshold, 1))

    return (bins, counts)


def find_sim_threshold_per_acc(
    mmseqs_log: Path,
    mmseqs_idt_threshold: float = 0.6,
    min_cov: float = 0.5,
    flattening_value: float = 0.1,
) -> float:
    """
    Given a similarity log computed by mmseqs2, fit a Gaussian distribution to the distribution of similarity values (histogram).
    Compute the mean and standard deviation of the Gaussian distribution.

    :param mmseqs_log: the path to the mmseqs2 log file.
    :param mmseqs_idt_threshold: the mmseqs2 similarity threshold.
    :param min_cov: the minimum scaffold-probe coverage.
    :param flattening_value: bins that are lower than the largest bin multiplied by the flattening_value are removed before the Gaussian curve fitting.
    :return: the lower value at 2 standard deviations from the mean that should encompass 95% of the data.
    """
    mmseqs_tables = pl.read_csv(mmseqs_log, separator="\t")
    df = mmseqs_tables.filter(pl.col("cov") > min_cov).select("idt")
    bins, counts = get_histogram_data(df, mmseqs_idt_threshold)
    min_count = max(counts) * flattening_value
    filtered_bins, filtered_counts = zip(
        *[(x, y) for x, y in zip(bins[:-1], counts) if y >= min_count]
    )
    initial = [np.median(filtered_bins), np.std(filtered_bins), max(filtered_counts)]
    gau_params, gau_covariance = curve_fit(
        gaussian, filtered_bins, filtered_counts, p0=initial
    )
    mu, sigma, _ = gau_params
    return mu, sigma


def snakemake_call(snakemake):
    mmseq_tables = snakemake.input
    out_path = Path(snakemake.output[0])

    min_idt = snakemake.params.min_idt
    min_cov = snakemake.params.min_cov
    p = snakemake.params.get("flattening_prop", 0.1)

    samples_fit = {
        sample.stem: find_sim_threshold_per_acc(sample, min_idt, min_cov, p)
        for sample in map(Path, mmseq_tables)
    }
    samples_thresholds = {
        sample: mu - 2 * abs(sigma) for sample, (mu, sigma) in samples_fit.items()
    }

    with out_path.open("w") as f:
        json.dump(samples_thresholds, f, indent=4)

    with open(snakemake.log[0], "w") as outlog:
        outlog.write("sample\tmu\t\tsigma\n")
        for sample, (mu, sigma) in samples_fit.items():
            outlog.write(f"{sample}\t{mu:.4f}\t{sigma:.4f}\n")


def main(): ...


if __name__ == "__main__":
    if "snakemake" in globals():
        snakemake_call(snakemake)
    else:
        main()
