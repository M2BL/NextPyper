from dataclasses import dataclass
from itertools import groupby
from operator import itemgetter
from pathlib import Path
import sys
from typing import NamedTuple
import argparse

from intervaltree import Interval, IntervalTree
import polars as pl
from Bio import SeqIO


class QueryHit(NamedTuple):
    query: str
    idt: float
    qlen: int
    qcov: float
    tlen: int
    tcov: float


def build_target_trees(
    df: pl.DataFrame, min_idt: float = 0.7
) -> dict[str, IntervalTree]:
    """Given a Dataframe of target hits, for each target build a target tree with the
    coverage of all the hits, while keeping only the best hit in each region. Hits
    with similarity below min_idt are ignored.
    Returns a dict mapping target names to IntervalTrees."""
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


@dataclass
class LongestPath:
    """
    Result of the get_longest_path recursion
    """

    path: int


def get_longest_path(itree: IntervalTree, min_overlap_length_int: int) -> int:
    """
    Perform a dfs on the interval in order to find the combination that yields the longest path
    Parameters
    ----------
    itree
    min_overlap_length_int: max overlap in nucleotides allowed between two fragments to fuse the paths

    Returns the length of the longest path
    -------

    """
    longest_path = LongestPath(0)
    all_intervals = sorted(itree, key=lambda i: i.begin)

    def dfs_util(
        current_fragment: Interval,
        remaining_fragments: list[Interval],
        longest_path: LongestPath,
        current_length: int = 0,
    ) -> None:
        if not remaining_fragments:
            if current_length > longest_path.path:
                longest_path.path = current_length
            return

        for next_fragment in remaining_fragments:
            if (
                overlap := current_fragment.overlap_size(next_fragment)
            ) <= min_overlap_length_int:
                current_length += next_fragment.data.tcov - overlap
                # ToDo: Check if it is working as expected (hint: tcov here is suspicious)

            next_fragment_idx = remaining_fragments.index(next_fragment)
            dfs_util(
                next_fragment,
                remaining_fragments[next_fragment_idx + 1 :],
                longest_path,
                current_length,
            )

    for start_frg in sorted(itree):
        start_frg_idx = all_intervals.index(start_frg)
        dfs_util(
            start_frg,
            all_intervals.copy()[start_frg_idx + 1 :],
            longest_path,
            start_frg.data.tcov,
        )

    return longest_path.path


def find_busco_category(
    itree: IntervalTree, min_length_percent: float, max_overlap_percent: float = 0
) -> int:
    """
    Category codes are defined as:
    0: Complete single copy
    1: Complete and duplicated
    2: Fragmented (single or multiple copies all fragmented)
    3: Missing
    Parameters
    ----------
    itree: intput interval tree
    min_length_percent: percentage of the target length covered by the query to qualify for completeness
    max_overlap_percent: percentage of the target overlap covered by the query allowed during fragment fusion
    ----------
    Returns the category code
    -------

    """
    if itree.is_empty():
        return 3

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
        if get_longest_path(itree, overlap_threshold) >= length_threshold:
            return 1
        return 0

    # Case check for fragmented
    if get_longest_path(itree, overlap_threshold) >= length_threshold:
        return 2

    return 3


def categorize_sample(
    hits: pl.DataFrame,
    chimera_df: pl.DataFrame,
    targets: list[str],
    min_hit_tcov: float = 0.25,
    min_length_cov: float = 0.7,
    min_idt: float = 0.99,
    include_borderline: bool = False,
) -> dict[str, int]:
    """
    Category codes are defined as:
    0: Complete single copy
    1: Complete and duplicated
    2: Fragmented (single or multiple copies all fragmented)
    3: Missing
    4: Recovered chimera (The query is chimeric but still meets acceptance criteria)
    5: Failed chimera (The query is chimeric and does not meet acceptance criteria)
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
    Returns the category code
    -------
    """

    df2 = (
        hits.filter(pl.col("tcov") >= min_hit_tcov)
        .sort(["query", "bits"], descending=True)
        .group_by("query")
        .agg(pl.all().first())
    )
    chimera_subset_df = chimera_df.select(
        pl.col("column_2"), pl.col("column_18")
    ).rename({"column_2": "query", "column_18": "chimera"})

    if include_borderline:
        df3 = df2.join(
            chimera_subset_df.with_columns_seq(chimera=pl.col("chimera") != "N"),
            on="query",
            how="left",
        ).fill_null(False)
    else:
        df3 = df2.join(
            chimera_subset_df.with_columns_seq(chimera=pl.col("chimera") == "Y"),
            on="query",
            how="left",
        ).fill_null(False)

    # Build target trees with and without chimeras
    chim_trees = build_target_trees(df3, min_idt)
    no_chim_trees = build_target_trees(df3.filter(~pl.col("chimera")), min_idt)

    # Compute the categories with both sets of trees
    chim_cat_dict = {
        target: find_busco_category(tree, min_length_cov, 0)
        for target, tree in chim_trees.items()
    }
    no_chim_cat_dict = {
        target: find_busco_category(tree, min_length_cov, 0)
        for target, tree in no_chim_trees.items()
    }

    final_cat = {}
    for target in targets:
        no_chim_cat = no_chim_cat_dict.get(target, None)
        chim_cat = chim_cat_dict.get(target, None)

        if no_chim_cat is None and chim_cat is None:
            final_cat[target] = 3
        elif no_chim_cat == chim_cat:
            final_cat[target] = no_chim_cat
        elif chim_cat != 3:
            final_cat[target] = 4
        else:
            final_cat[target] = 5

    return final_cat


def parse_args():
    parser = argparse.ArgumentParser(description="Target benchmarking script")
    parser.add_argument("hits", type=Path, help="Path to MMseqs hits table (.tsv)")
    parser.add_argument(
        "chimeras", type=Path, help="Path to vsearch chimera table (.tsv)"
    )
    parser.add_argument("targets", type=Path, help="Path to targets (.fasta)")
    parser.add_argument("output", type=Path, help="Path to output table (.tsv)")
    parser.add_argument(
        "--min_tcov",
        type=float,
        default=0.25,
        help="Minimum target coverage of the query to consider a hit (default: 0.25)",
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
    return parser.parse_args()


def main():
    args = parse_args()

    hits_file = args.hits
    chimera_file = args.chimeras
    targets_file = args.targets
    output_file = args.output

    min_tcov = args.min_tcov
    min_tot_cov = args.min_tot_cov
    min_idt = args.min_idt
    include_borderline = args.include_borderline

    if not hits_file.exists() or not chimera_file.exists() or not targets_file.exists():
        print("One or more input files do not exist.")
        sys.exit(1)
    if (
        not hits_file.is_file()
        or not chimera_file.is_file()
        or not targets_file.is_file()
    ):
        print("One or more input paths are not files.")
        sys.exit(1)
    if not output_file.parent.exists():
        output_file.parent.mkdir(parents=True)

    df = pl.read_csv(hits_file, separator="\t", has_header=True)
    chimera_df = pl.read_csv(chimera_file, separator="\t", has_header=False)
    targets = [rec.id for rec in SeqIO.parse(targets_file, "fasta")]

    with open(output_file, "w") as out_file:
        out_file.write("target\tcategory\n")
        categories = categorize_sample(
            df,
            chimera_df,
            targets,
            min_hit_tcov=min_tcov,
            min_length_cov=min_tot_cov,
            min_idt=min_idt,
            include_borderline=include_borderline,
        )
        for target, category in categories.items():
            out_file.write(f"{target}\t{category}\n")


if __name__ == "__main__":
    main()
