#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2024
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.

# ToDo: Add docstring explaining the classification of simulated target sequences into categories.
"""
Category codes are defined as:

0: Complete single copy
1: Complete and duplicated
2: Fragmented (single or multiple copies all fragmented)
3: Missing
4: Recovered chimera (The query is chimeric but still meets acceptance criteria)
5: Failed chimera (The query is chimeric and does not meet acceptance criteria)

An additional pseudo-category is computed as:
6: Noise sequences (count of all the query sequences that were not matched to any targets)

"""

# =======================================================================================
#               IMPORTS
# =======================================================================================

from __future__ import annotations
from collections import Counter
from itertools import chain, groupby, pairwise
from operator import itemgetter
from pathlib import Path
import sys
from typing import NamedTuple
import argparse

from intervaltree import Interval, IntervalTree
import polars as pl
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# =======================================================================================
#               CONSTANTS
# =======================================================================================

# Columns parsed from magicblast output table
MAGIC_COLS = {
    "# Fields: query acc.": "query",
    "reference acc.": "target",
    "% identity": "fident",
    "query start": "qstart",
    "query end": "qend",
    "reference start": "tstart",
    "reference end": "tend",
    "score": "score",
    "query length": "qlen",
}

# Columns used for building the target trees
TREE_COLS = [
    "query",
    "target",
    "tstart",
    "tend",
    "fident",
    "qlen",
    "qcov",
    "tlen",
    "tcov",
]

# =======================================================================================
#               CLASSES
# =======================================================================================


class QueryHit(NamedTuple):
    query: str
    idt: float
    qlen: int
    qcov: float
    tlen: int
    tcov: float


class LongestPath:
    """Result of the get_longest_path recursion"""

    path: list[Interval]

    def __init__(self, *args: Interval):
        self.path = list(args)

    def __getitem__(
        self, index: int | slice | list[int | slice]
    ) -> Interval | LongestPath:
        match index:
            case int():
                return self.path[index]
            case slice():
                return LongestPath(*self.path[index])
            case [*idxs] if all(isinstance(idx, (int, slice)) for idx in idxs):
                return LongestPath(*chain.from_iterable(self.path[sli] for sli in idxs))
            case _:
                raise NotImplementedError(f"Given {index}")

    def append(self, inter: Interval) -> None:
        self.path.append(inter)

    @property
    def path_length(self) -> int:
        return sum(inter.length() for inter in self.path) - sum(
            left.overlap_size(right) for left, right in pairwise(self.path)
        )

    @property
    def global_idt(self) -> float:
        """Calculate weighted average identity based on interval coverage."""  # ToDo: Correct for the idt of overlaps

        if not self.path:
            return 0.0

        return (
            sum(inter.data.idt * inter.length() for inter in self.path)
            / self.path_length
        )


# =======================================================================================
#               CLASSES
# =======================================================================================


def build_target_trees(
    df: pl.DataFrame, min_idt: float = 0.7
) -> dict[str, IntervalTree]:
    """Given a Dataframe of target hits, for each target build a target tree with the
    coverage of all the hits, while keeping only the best hit in each region. Hits
    with similarity below min_idt are ignored.
    Returns a dict mapping target names to IntervalTrees."""

    target_hits = df.sort("target").select(TREE_COLS)

    target_trees = {}
    for target, phits in groupby(target_hits.iter_rows(), itemgetter(1)):
        target_trees[target] = IntervalTree(
            (
                Interval(tstart - 1, tend, info)
                if tstart < tend
                else Interval(tend - 1, tstart, info)
            )
            for query, _, tstart, tend, idt, qlen, qcov, tlen, tcov in phits
            if (info := QueryHit(query, idt, qlen, qcov, tlen, tcov)).idt >= min_idt
        )

    return target_trees


def get_longest_path(itree: IntervalTree, max_ovlp: int, min_idt: float = 0.0) -> int:
    """
    Perform a dfs on the interval in order to find the combination that yields the longest path
    Parameters
    ----------
    itree
    max_ovlp: max overlap in nucleotides allowed between two fragments to fuse the paths

    Returns the length of the longest path
    -------

    """
    longest_path = LongestPath()

    def dfs_util(
        current_fragment: Interval,
        remaining_fragments: list[Interval],
        longest_path: LongestPath,
        current_path: LongestPath,
    ) -> None:

        if (
            current_path.path_length > longest_path.path_length
            and current_path.global_idt >= min_idt
        ):
            longest_path.path = current_path.path

        for next_fragment in remaining_fragments:
            if current_fragment.overlap_size(next_fragment) <= max_ovlp:
                current_path.append(next_fragment)

            no_ovlp_fragments = sorted(
                set(itree)
                - set(itree.overlap(0 + max_ovlp, next_fragment.end - max_ovlp))
            )

            dfs_util(
                next_fragment,
                no_ovlp_fragments,
                longest_path,
                current_path[:],
            )

    for start_frg in sorted(itree):
        remaining_fragments = sorted(
            set(itree) - set(itree.overlap(0 + max_ovlp, start_frg.end - max_ovlp))
        )

        dfs_util(
            start_frg,
            remaining_fragments,
            longest_path,
            LongestPath(start_frg),
        )

    return longest_path.path_length


def find_busco_category(
    itree: IntervalTree,
    min_length_percent: float,
    max_overlap_percent: float = 0,
    min_idt: float = 0.0,
) -> int:
    """
    Category codes are defined as:
    0: Complete single copy
    1: Complete and duplicated
    2: Fragmented success (single or multiple copies all fragmented, fulfill acceptance criteria)
    3: Fragmented failure (single or multiple fragments detected, that do not fulfill acceptance criteria)
    4: Missing
    Parameters
    ----------
    itree: intput interval tree
    min_length_percent: percentage of the target length covered by the query to qualify for completeness
    max_overlap_percent: percentage of the target overlap covered by the query allowed during fragment fusion
    ----------
    Returns the category code
    -------

    """
    itree = itree.copy()  # avoid modifying the original tree
    if itree.is_empty():
        return 4

    interval_0 = list(itree.all_intervals)[0]
    length_threshold = interval_0.data.tlen * min_length_percent
    overlap_threshold = interval_0.data.qlen * max_overlap_percent

    # Case complete single copy
    if len(itree) == 1:
        if interval_0.data.tcov >= min_length_percent:
            return 0
        return 3

    longest = max(itree, key=lambda i: (i.data.tcov, i.data.idt))

    # Case either complete single copy or multiple copies
    if longest.data.tcov >= min_length_percent:
        itree.remove(longest)
        if get_longest_path(itree, overlap_threshold, min_idt) >= length_threshold:
            return 1
        return 0

    # Case check for fragmented
    if get_longest_path(itree, overlap_threshold, min_idt) >= length_threshold:
        return 2

    return 3


def compute_target_covs(target_trees: dict[str, IntervalTree]) -> pl.DataFrame:
    """Compute the covered bases of the targets based on their interval trees."""

    def tcov(tree: IntervalTree) -> int:
        t = tree.copy()
        t.merge_overlaps()
        return sum(inter.length() for inter in t)

    return pl.DataFrame(
        [(probe, tcov(tree)) for probe, tree in target_trees.items()],
        {"target": str, "cov": int},
        orient="row",
    )


def categorize_sample(
    hits: pl.DataFrame,
    chimera_df: pl.DataFrame,
    targets: dict[str, SeqRecord],
    min_hit_tcov: float = 0.0,
    min_length_cov: float = 0.7,
    min_idt: float = 0.99,
    include_borderline: bool = False,
) -> tuple[dict[str, int], int, int, pl.DataFrame]:
    """
    Computes target-centric categories for a sample based on magicblast hits from a given assembly
    and the reference targets.

    Category codes are defined as:
    0: Complete single copy
    1: Complete and duplicated
    2: Fragmented success (single or multiple copies all fragmented)
    3: Fragmented failure (single or multiple fragments detected, that do not fulfill acceptance criteria)
    4: Missing

    Assembly-centric categories are also computed: Number of chimeric sequences found by
    vsearch and the number of noise sequences (sequences that do not hit any target or do not meet the idt criteria).

    Parameters
    ----------
    hits: mmseqs hits dataframe.
    chimera_df: dataframe with vsearch chimeric results.
    targets: list of all the targets (including the missing ones).
    min_hit_tcov: minimum target coverage of the QUERY to consider a hit.
    min_idt: minimum identity of the hits on the target to accept it.
    min_length_cov: minimum total (considering multiple hits if needed) target coverage to consider the target recovered.
    include_borderline: whether to include borderline chimeric cases as chimeras (see vsearch documentation)
    ----------
    Returns the category code dictionary, number of chimeras and number of noise sequences.
    -------
    """

    # Get the length of the targets
    target_lens = pl.DataFrame(
        [(target, len(rec)) for target, rec in targets.items()],
        {"target": str, "tlen": int},
        orient="row",
    )

    # Add tlen, tcov and qcov to the hits table
    predf = (
        hits.select(*MAGIC_COLS.keys())
        .rename(MAGIC_COLS)
        .join(target_lens, on="target")
        .with_columns(
            fident=pl.col("fident") / 100,
            tcov=(pl.col("tend") - pl.col("tstart")).abs() / pl.col("tlen"),
            qcov=(pl.col("qend") - pl.col("qstart")).abs() / pl.col("qlen"),
        )
    )

    df2 = predf.filter(pl.col("tcov") >= min_hit_tcov).filter(
        pl.col("score") == pl.col("score").max().over("query")
    )
    chimera_subset_df = chimera_df.select(
        pl.col("column_2"), pl.col("column_18")
    ).rename({"column_2": "query", "column_18": "chimera"})

    if include_borderline:
        final_chimera_df = chimera_subset_df.with_columns_seq(
            chimera=pl.col("chimera") != "N"
        ).filter(pl.col("chimera"))
    else:
        final_chimera_df = chimera_subset_df.with_columns_seq(
            chimera=pl.col("chimera") == "Y"
        ).filter(pl.col("chimera"))

    df3 = df2.join(final_chimera_df, on="query", how="anti")
    n_chimeras = len(final_chimera_df)

    # Build target trees with and without chimeras
    no_chim_trees = build_target_trees(df3, 0.0)

    # Compute the number of "noise" sequences
    # All unaligned sequences count
    noise = len(hits[:, 0].unique()) - len(predf["query"].unique())
    # Add sequences that won't be assigned to any target tree
    noise += len(df3.filter(pl.col("fident") < 0.0))

    # Compute the categories from the trees
    no_chim_cat_dict = {
        target: find_busco_category(tree, min_length_cov, 0, min_idt)
        for target, tree in no_chim_trees.items()
    }

    final_cat = {target: no_chim_cat_dict.get(target, 4) for target in targets}

    cov_df = compute_target_covs(no_chim_trees)
    final_cov_df = target_lens.join(cov_df, on="target", how="left").fill_null(0)
    final_cov_df = final_cov_df.select(pl.sum("tlen"), pl.sum("cov")).with_columns(
        tcov=pl.col("cov") / pl.col("tlen")
    )

    return final_cat, n_chimeras, noise, final_cov_df


def parse_args():
    parser = argparse.ArgumentParser(description="Target benchmarking script")
    parser.add_argument("hits", type=Path, help="Path to MMseqs hits table (.tsv)")
    parser.add_argument(
        "chimeras", type=Path, help="Path to vsearch chimera table (.tsv)"
    )
    parser.add_argument("targets", type=Path, help="Path to targets (.fasta)")
    parser.add_argument(
        "prefix", type=Path, help="Prefix for output tables (categories and coverages)"
    )
    parser.add_argument(
        "--min_tcov",
        type=float,
        default=0.0,
        help="Minimum target coverage of the query to consider a hit (default: 0.0)",
    )
    parser.add_argument(
        "--min_tot_cov",
        type=float,
        default=0.7,
        help="Minimum total target coverage to consider the target recovered (default: 0.7)",
    )
    parser.add_argument(
        "--min_idt",
        type=float,
        default=0.99,
        help="Minimum identity of the hits on the target to accept it (default: 0.99)",
    )
    parser.add_argument(
        "--include_borderline",
        action="store_true",
        help="Include borderline chimeric cases as chimeras (see vsearch documentation)",
    )
    parser.add_argument(
        "--batch",
        action="store_true",
        help="Process a batch of samples. Inputs (hits, chimeras, targets) are expected to be folders containing files with the same names.",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    hits_path = args.hits
    chimera_path = args.chimeras
    targets_path = args.targets
    prefix = args.prefix

    output_path = prefix.with_stem(f"{prefix.stem}_cat.tsv")
    cov_path = prefix.with_stem(f"{prefix.stem}_cov.tsv")

    min_tcov = args.min_tcov
    min_tot_cov = args.min_tot_cov
    min_idt = args.min_idt
    include_borderline = args.include_borderline
    batch = args.batch

    if not hits_path.exists() or not chimera_path.exists() or not targets_path.exists():
        print("One or more input files do not exist.")
        sys.exit(1)

    # Single sample processing
    if not batch:
        if (
            not hits_path.is_file()
            or not chimera_path.is_file()
            or not targets_path.is_file()
        ):
            print("One or more input paths are not files.")
            sys.exit(1)
        if not output_path.parent.exists():
            output_path.parent.mkdir(parents=True)

        df = pl.read_csv(
            hits_path,
            separator="\t",
            has_header=True,
            skip_rows=2,
            schema_overrides={"% identity": pl.Float64},
            infer_schema_length=10000,
        )
        chimera_df = pl.read_csv(
            chimera_path, separator="\t", has_header=False, infer_schema_length=10000
        )
        targets = SeqIO.to_dict(SeqIO.parse(targets_path, "fasta"))
        categories, n_chimeras, noise, cov_df = categorize_sample(
            df,
            chimera_df,
            targets,
            min_hit_tcov=min_tcov,
            min_length_cov=min_tot_cov,
            min_idt=min_idt,
            include_borderline=include_borderline,
        )

        with open(output_path, "w") as out_file:
            out_file.write("target\tcategory\n")
            for target, category in categories.items():
                out_file.write(f"{target}\t{category}\n")

    # Batch processing (unified table)
    else:
        if (
            not hits_path.is_dir()
            or not chimera_path.is_dir()
            or not targets_path.is_dir()
        ):
            print("One or more input paths are not directories.")
            sys.exit(1)

        with open(output_path, "w") as out_file, open(cov_path, "w") as cov_out:
            for targets_file in sorted(targets_path.glob("*.fasta")):
                sample_name = targets_file.stem
                hits_file = hits_path / f"{sample_name}.tsv"
                chimera_file = chimera_path / f"{sample_name}.tsv"

                if not hits_file.exists() or not chimera_file.exists():
                    print(f"Missing files for sample {sample_name}. Skipping.")
                    continue

                df = pl.read_csv(
                    hits_file,
                    separator="\t",
                    has_header=True,
                    skip_rows=2,
                    schema_overrides={"% identity": pl.Float64},
                    infer_schema_length=10000,
                )
                chimera_df = pl.read_csv(
                    chimera_file,
                    separator="\t",
                    has_header=False,
                    infer_schema_length=10000,
                )
                targets = SeqIO.to_dict(SeqIO.parse(targets_file, "fasta"))
                categories, n_chimeras, noise, cov_df = categorize_sample(
                    df,
                    chimera_df,
                    targets,
                    min_hit_tcov=min_tcov,
                    min_length_cov=min_tot_cov,
                    min_idt=min_idt,
                    include_borderline=include_borderline,
                )

                # Summarize categories and write to output
                category_counts = Counter(categories.values())

                # Add number of chimeras as extra assembly-centric category
                category_counts[5] = n_chimeras

                # Add a last noise assembly-centric category with the counts of
                # extra sequences that do not hit any target or meet idt criteria.
                category_counts[6] = noise

                # Finally, add the total number of sequences in the assembly
                category_counts[7] = len(chimera_df)

                categories = "\t".join(
                    str(category_counts.get(i, 0))
                    for i in range(max(category_counts) + 1)
                )
                out_file.write(f"{sample_name}\t{categories}\n")

                # If requested, write also the length based recovery
                cov_df.insert_column(0, pl.lit(sample_name)).write_csv(
                    cov_out, separator="\t", include_header=False, float_precision=5
                )


if __name__ == "__main__":
    main()
