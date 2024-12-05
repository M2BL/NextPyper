#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2024
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.

# ToDo: Expand on the assumptions and usage
"""
Functions to summarize the results from a NextPyper run. The main function of the module is
summarize workflow, and given the ouput directory of a successful NextPyper run, it creates
a table with the following metrics about the run (for each sample):

    - raw_reads: Number of raw reads in the sample.
    - trimmed_reads: Number of reads after trimming.
    - cleaned_reads: Number of reads after cleaning (rRNA + [cpDNA]).
    - probe_matching_reads: Number of reads after probe matching filtering.
    - probes_in_assembly: Number of probes that had a hit in the assembly.
    - probes_at_25pct: Number of probes in the final alignments that have at least 25% of the probe length.
    - probes_at_50pct: Number of probes in the final alignments that have at least 50% of the probe length.
    - probes_at_75pct: Number of probes in the final alignments that have at least 75% of the probe length.
    - probes_at_100pct: Number of probes in the final alignments that have at least 100% of the probe length.
    - probes_at_125pct: Number of probes in the final alignments that have at least 125% of the probe length.

Example of usage:

>>> results_dir = Path("my_directory/nextpyper_output")
>>> out_stats = Path("my_directory/stats_nextpyper.csv")
>>> stats = summarize_workflow(results_dir)
>>> stats.to_csv(out_stats, index=False)

Will create a csv table "stats_nextpyper.csv" summarizing the run.
"""

__version__ = "0.1"

# =======================================================================================
#               IMPORTS
# =======================================================================================

from dataclasses import dataclass, field
from itertools import groupby
from functools import reduce
from pathlib import Path
import json
import re
import sys

from Bio import SeqIO
import pandas as pd
import numpy as np


# =======================================================================================
#               CONSTANTS
# =======================================================================================

MATCHING_PATTERN = r"(Input|Contaminants):\s+(\d+)"
CLEANING_PATTERN = r"(Input|Result):\s+(\d+)"

SCF_PATTERN = re.compile(
    r"^(?P<sample>.*?)-(?P<probe>.*?)_(?P<cluster>\d+?)_(?P<seed>\d+?)", re.VERBOSE
)

# =======================================================================================
#               CLASSES
# =======================================================================================


@dataclass
class AlnTable:
    samples: list[str]
    probe_lens: dict[str, int]
    mat: pd.DataFrame = field(init=False, repr=False)

    def __post_init__(self):
        seq_lens = np.zeros((len(self.samples), len(self.probe_lens)))
        self.mat = pd.DataFrame(
            seq_lens, index=self.samples, columns=self.probe_lens, dtype=int
        )

    def update_counts(self, counts: dict[str, int], probe: str) -> None:

        if not counts:
            return

        if probe not in self.mat.columns:
            raise ValueError(f"{probe=} not in probe set")

        rows, vals = zip(*counts.items())
        self.mat.loc[rows, probe] = vals

    def make_normalized_table(self) -> pd.DataFrame:
        return self.mat / (np.array(list(self.probe_lens.values())) * 3)


# =======================================================================================
#               FUNCTIONS
# =======================================================================================


def _get_raw_reads(sample: str, rootdir: Path) -> int:
    with (rootdir / f"logs/preprocessing/fastp/{sample}.json").open() as file:
        return json.load(file)["summary"]["before_filtering"]["total_reads"]


def _get_cleaned_reads(sample: str, rootdir: Path, pat: str) -> tuple[int, int]:
    """Extract from bbduk logs the number of input and output reads. Depending on the
    pattern, the output would be tha matching (contaminants) or non matching (clean reads).
    """

    read_ext = re.compile(pat)
    with (rootdir / f"logs/preprocessing/bbduk_cleaning/{sample}.log").open() as file:
        return tuple(int(match[2]) for line in file if (match := read_ext.match(line)))


def _get_matched_reads(sample: str, rootdir: Path, pat: str) -> tuple[int, int]:
    """Extract from bbduk logs the number of input and output reads. Depending on the
    pattern, the output would be tha matching (contaminants) or non matching (clean reads).
    """

    read_ext = re.compile(pat)
    with (
        rootdir / f"logs/preprocessing/bbduk_probe_matching/{sample}.log"
    ).open() as file:
        return tuple(int(match[2]) for line in file if (match := read_ext.match(line)))


def _get_matched_probes(sample: str, rootdir: Path) -> int:
    """Return the number of probes found for a given sample after spades
    assembly and graph splitting."""

    return len(
        {
            file.stem.rsplit("_", 1)[0]
            for file in (rootdir / f"assembled/split_components/{sample}").glob(
                "*.fasta"
            )
        }
    )


def _get_len_probes(file: Path) -> int:
    "Given a file with multiple versions of a probe, return the size of the longest."

    return (
        np.fromiter((len(rec) for rec in SeqIO.parse(file, "fasta")), dtype=int)
    ).max()


def _get_seqs_lens(aln: Path, pat: re.Pattern) -> dict[str, int]:
    """Given a an alignment file with multiple sequences, group them by sample using
    the provided pattern and take the longest for each sample. Return a dictionary
    sample: length pairs."""

    get_sample = lambda rec: pat.search(rec.id)["sample"]
    return {
        sample: max(map(len, recs))
        for sample, recs in groupby(
            sorted(SeqIO.parse(aln, "fasta"), key=get_sample), key=get_sample
        )
    }


def _condense_counts(mul_counts: list[dict[str, int]]) -> dict[str, int]:
    "Given a list of dictionaries of sample:count with overlapping samples, returned a unified dictionary with the sum of the counts."

    return {
        key: sum(counts.get(key, 0) for counts in mul_counts)
        for key in reduce(set.union, mul_counts, set())
    }


def summarize_workflow(results_dir: Path) -> pd.DataFrame:
    """Given the root path of a NextPyper output, collect and summarize the results of the run.
    Returns a dataframe with a row per sample and multiple statistics.

    Current columns include:

    - raw_reads: Number of raw reads in the sample.
    - trimmed_reads: Number of reads after trimming.
    - cleaned_reads: Number of reads after cleaning (rRNA + [cpDNA]).
    - probe_matching_reads: Number of reads after probe matching filtering.
    - probes_in_assembly: Number of probes that had a hit in the assembly.
    - probes_at_25pct: Number of probes in the final alignments that have at least 25% of the probe length.
    - probes_at_50pct: Number of probes in the final alignments that have at least 50% of the probe length.
    - probes_at_75pct: Number of probes in the final alignments that have at least 75% of the probe length.
    - probes_at_100pct: Number of probes in the final alignments that have at least 100% of the probe length.
    - probes_at_125pct: Number of probes in the final alignments that have at least 125% of the probe length.
    """

    ## Get the samples processed
    samples = sorted(
        sample.name for sample in (Path(results_dir) / "assembled/spades").iterdir()
    )
    stats = pd.DataFrame(index=samples)

    ## Get the raw reads
    stats["raw_reads"] = {
        sample: _get_raw_reads(sample, results_dir) for sample in samples
    }

    ## Get Preprocessing stats:
    clean_stats = {
        sample: _get_cleaned_reads(sample, results_dir, CLEANING_PATTERN)
        for sample in samples
    }
    stats["trimmed_reads"] = {
        sample: trimmed for sample, (trimmed, cleaned) in clean_stats.items()
    }
    stats["cleaned_reads"] = {
        sample: cleaned for sample, (trimmed, cleaned) in clean_stats.items()
    }
    if (results_dir / "logs/preprocessing/bbduk_probe_matching").exists():
        stats["probe_matching_reads"] = {
            sample: _get_matched_reads(sample, results_dir, MATCHING_PATTERN)[1]
            for sample in samples
        }

    ## Get matched probes
    stats["probes_in_assembly"] = {
        sample: _get_matched_probes(sample, results_dir) for sample in samples
    }

    ## Get lens in aa
    probe_lens = {
        file.stem: _get_len_probes(file)
        for file in sorted(
            (results_dir / "homolog_prospection/region_separation/input_probes").glob(
                "*.fasta"
            )
        )
    }

    ## Make an auxiliary matrix of [Samples x probes]
    table = AlnTable(samples, probe_lens)

    for probe_dir in (
        results_dir / "homolog_prospection/region_separation/separation_output/scfs"
    ).iterdir():
        if probe_dir.is_file():
            continue

        if counts := _condense_counts(
            [_get_seqs_lens(aln, SCF_PATTERN) for aln in probe_dir.glob("*.fasta")]
        ):
            table.update_counts(counts, probe=probe_dir.name)

    norm_table = table.make_normalized_table()

    for level in (25, 50, 75, 100, 125):
        stats[f"probes_at_{level}pct"] = (norm_table > level / 100).sum(axis=1)

    return stats.reset_index().rename({"index": "samples"}, axis=1), norm_table


def snakemake_call(snakemake):
    with open(snakemake.log[0], "w") as outlog:
        sys.stdout = sys.stderr = outlog

        wf_dir = Path(snakemake.input.workflowdir)
        out_stats = Path(snakemake.output[0])

        df, _ = summarize_workflow(wf_dir)
        df.to_csv(out_stats, index=False)


def main():
    results_dir = sys.argv[1]
    out_file = sys.argv[2]

    df, table = summarize_workflow(Path(results_dir))
    df.to_csv(out_file, index=False)

    try:
        tab_file = sys.argv[3]
        table.T.to_csv(tab_file, float_format="%.2f")
    except IndexError:
        pass


if __name__ == "__main__":
    if "snakemake" in globals():
        snakemake_call(snakemake)
    else:
        main()
