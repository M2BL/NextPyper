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
from more_itertools import one


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


def collapse_alleles_df(df: pl.DataFrame) -> pl.DataFrame:
    """Group the sequences from the same component with the same length and pick the one
    with the highest kmer count. It is assumed that these sequences are different alleles.

    Return the dataframe with the collapsed alleles."""

    return df.group_by("probe", "seed", "seed_id", "comp", "len").agg(
        pl.col("query").gather(pl.col("kmers").arg_max()).first(),
        pl.col("kmers").max(),
        count=pl.len(),
    )


def collapse_variants_df(df: pl.DataFrame, max_vars: int = 1) -> pl.DataFrame:
    """Given a dataframe with already collapsed alleles, further collapse the variants
    in the same component.

    Return a DataFrame with the collapsed variants
    """
    return (
        df.group_by("probe", "seed", "seed_id", "comp")
        .agg(
            pl.all().sort_by("kmers", descending=True).head(max_vars),
            var_count=pl.len(),
            allele_count=pl.sum("count"),
        )
        .explode("len", "query", "kmers", "count")
    )


def collapse_records(recs: list[SeqRecord], df: pl.DataFrame) -> list[SeqRecord]:
    """Group the sequences from the same component with the same length and pick the one
    with the highest kmer count. It is assumed that these sequences are different alleles.

    Return a reduced set of records with the selected alleles."""

    return [rec for rec in recs if rec.id in set(df["query"])]


def split_explosive_probes(
    df: pl.DataFrame, explosive_limit: int
) -> tuple[list[SeqRecord], list[SeqRecord]]:

    allele_df = collapse_alleles_df(df)
    var_df = collapse_variants_df(allele_df)

    explosive_probes = (
        var_df.filter(pl.col("var_count") >= explosive_limit).select("probe").unique()
    )

    print(f"Total probes: {var_df["probe"].n_unique()}. Setting {explosive_limit=}")
    print(f"Normal probes: {var_df["probe"].n_unique() - len(explosive_probes)}")
    print(f"Explosive probes: {len(explosive_probes)}")

    normal_df = df.join(explosive_probes, on="probe", how="anti")
    expl_df = df.join(explosive_probes, on="probe")

    return normal_df, expl_df


def snakemake_call(snakemake):
    with open(snakemake.log[0], "w") as outlog:
        sys.stdout = sys.stderr = outlog

        records_path = Path(snakemake.input[0])
        out_path = Path(snakemake.output.normal)
        expl_out_path = snakemake.output.get("expl")

        # Mandatory parameters
        pattern = snakemake.params.get("pattern", TARGET_PAT)
        mode = snakemake.params.get("mode", "collapse")

        # Optional parameters (required according to mode)
        empty_ok = snakemake.params.get("empty_ok", False)
        collapse_vars = snakemake.params.get("collapse_vars", False)
        max_vars = one(snakemake.params.get("max_vars", [10]))

        recs = list(SeqIO.parse(records_path, "fasta"))

        print(f"Starting sequences: {len(recs)}")
        print(f"{mode=}")

        # If the input is empty and we allow for it, just touch the output
        if empty_ok and len(recs) == 0:
            print("Exiting without error due to empty input.")
            out_path.touch()
            sys.exit(0)

        df = query2df(recs, pattern, FIELDS)

        match mode:
            case "collapse":
                print(f"{collapse_vars=}")
                allele_df = collapse_alleles_df(df)
                print(f"Sequences after allele collapsing: {len(allele_df)}")

                if collapse_vars:
                    var_df = collapse_variants_df(allele_df, max_vars)
                    new_recs = collapse_records(recs, var_df)
                    print(f"Sequences after variant collapsing: {len(var_df)}")
                else:
                    new_recs = collapse_records(recs, allele_df)

                SeqIO.write(new_recs, out_path, "fasta")

            case "split" | "cap":
                if mode == "split" and expl_out_path is None:
                    raise ValueError(f"For {mode=}, two outputs have to be specified.")

                print("Splitting into normal and explosive sequence sets ({max_vars})")
                normal_df, explosive_df = split_explosive_probes(df, max_vars)

                print(f"Sequences in normal set: {len(normal_df)}")
                print(f"Sequences in explosive set: {len(explosive_df)}")
                expl_allele_df = collapse_alleles_df(explosive_df)
                print(
                    f"Sequences in explosive set (allele collapsed): {len(expl_allele_df)}"
                )

                if mode == "cap":
                    expl_var_df = collapse_variants_df(expl_allele_df, max_vars)
                    print(
                        f"Sequences in explosive set (after capping to {max_vars=}): {len(expl_var_df)}"
                    )

                    final_df = normal_df.select("query").vstack(
                        expl_var_df.select("query")
                    )
                    final_recs = collapse_records(recs, final_df)

                    SeqIO.write(final_recs, out_path, "fasta")
                else:

                    normal_recs = collapse_records(recs, normal_df)
                    explosive_recs = collapse_records(recs, expl_allele_df)

                    SeqIO.write(normal_recs, out_path, "fasta")
                    SeqIO.write(explosive_recs, Path(expl_out_path), "fasta")

            case _:
                raise NotImplementedError(f"{mode=} not implemented.")


def main():
    if len(sys.argv) != 3:
        print("Usage: python var_asm_parser.py <saute_asm.fasta> <output.fasta>")
        sys.exit(1)

    records_path = Path(sys.argv[1])
    out_path = Path(sys.argv[2])

    recs = list(SeqIO.parse(records_path, "fasta"))
    df = query2df(recs, TARGET_PAT, FIELDS)
    new_recs = collapse_records(recs, df)
    SeqIO.write(new_recs, out_path, "fasta")


if __name__ == "__main__":
    if "snakemake" in globals():
        snakemake_call(snakemake)
    else:
        main()
