#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2024
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.

"""Basic stats computation tribbles in a nextpyper run."""

# =======================================================================================
#               IMPORTS
# =======================================================================================

from pathlib import Path
import argparse
import yaml
import sys

import dpath
from Bio import SeqIO
import polars as pl

sys.path.append(str((Path(__file__) / "../../../src").resolve()))

from var_asm_parser import collapse_alleles_df, collapse_variants_df, query2df

# =============================================================================
#                CONSTANTS
# =============================================================================

from var_asm_parser import TARGET_PAT, FIELDS

TRIBBLE_LIM_DPATH = "nextpyper/pipeline/saute/reassembly/explosive_limit"


# =============================================================================
#                FUNCTIONS
# =============================================================================


def read_target_vars(
    vars: Path, tribble_lim: int, header_pat: str = TARGET_PAT
) -> pl.DataFrame:
    df = query2df(SeqIO.parse(vars, "fasta"), header_pat, FIELDS)
    return (
        collapse_variants_df(collapse_alleles_df(df))
        .with_columns(pl.col("probe").cast(pl.String))
        .filter(pl.col("var_count") > tribble_lim)
        .group_by("probe")
        .agg(tribble=pl.len())
    )


def find_tribbles(run_dir: Path, sample: str, tribble_lim: int = 30) -> pl.DataFrame:
    """Compute the tribbles for a single sample"""

    before_vars = run_dir / f"saute/target_assembly/{sample}/target_vars.fasta"
    after_vars = run_dir / f"saute/final/collected/{sample}.fasta"

    tribble_df = read_target_vars(before_vars, tribble_lim)

    if after_vars.exists():
        tribble_df2 = read_target_vars(after_vars, tribble_lim)
        merged = tribble_df.join(tribble_df2, on="probe", how="left", suffix="_2nd")
        merged = merged.select(
            pl.col("probe").len(),
            pl.sum("tribble"),
            (pl.col("tribble_2nd") > 0).sum().alias("probe_2nd"),
            pl.sum("tribble_2nd"),
        )
    else:
        merged = tribble_df.select(
            pl.col("probe").len(), pl.sum("tribble")
        ).with_columns(probe_2nd=None, tribble_2nd=None)

    return merged.insert_column(0, pl.lit(sample).alias("sample"))


def summarize_tribbles(run_dir: Path, tribble_lim: int | None = None) -> None:
    """Compute and summarize tribble stats for a whole run. Print the results table to stdout.
    The table includes the counts of probes with tribbles, and number of tribbles for both
    the first assembly and second assembly (secondary)."""

    with (run_dir / "config.yaml").open() as file:
        config = yaml.load(file, yaml.BaseLoader)
        try:
            reasm = dpath.get(config, "nextpyper/args/reasm")
        except KeyError:
            reasm = True

    if tribble_lim is None:
        tribble_lim = int(dpath.get(config, TRIBBLE_LIM_DPATH))

    df = pl.concat(
        find_tribbles(run_dir, sample.stem, tribble_lim)
        for sample in (run_dir / "saute/target_assembly").iterdir()
    )

    print(f"# {tribble_lim=}, {reasm=}")
    df.sort("sample").write_csv(sys.stdout, separator="\t")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute tribble stats of a NextPyper run to stdout as a tsv"
    )
    parser.add_argument("rundir", type=Path, help="Path to rundir")
    parser.add_argument(
        "-t",
        "--tribble-lim",
        type=int,
        default=None,
        help="""Maximum number of variants in a component to consider it a tribble. By default,
                read the limit used in the run.""",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    summarize_tribbles(args.rundir, args.tribble_lim)


if __name__ == "__main__":
    main()
