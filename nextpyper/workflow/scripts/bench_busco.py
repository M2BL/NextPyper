from dataclasses import dataclass
from itertools import groupby
from operator import itemgetter
from pathlib import Path
import sys
from typing import NamedTuple

from intervaltree import Interval, IntervalTree
import polars as pl

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
    print("building target trees...")
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
    path:int

def get_longest_path(itree:IntervalTree, min_overlap_length_int) -> int:
    longest_path = LongestPath(0)
    all_intervals = sorted(list(itree.all_intervals), key=lambda i: i.begin)

    def dfs_util(
            current_fragment: Interval, remaining_fragments: list[Interval], longest_path:LongestPath, current_length: int = 0,
    ) -> None:
        if not remaining_fragments:
            if current_length > longest_path.path:
                longest_path.path = current_length
            return

        for next_fragment in remaining_fragments:
            if (overlap := current_fragment.overlap_size(next_fragment)) <= min_overlap_length_int:
                current_length += next_fragment.data.tcov - overlap

            next_fragment_idx = remaining_fragments.index(next_fragment)
            dfs_util(next_fragment, remaining_fragments[next_fragment_idx+1:], longest_path, current_length)
    for start_frg in sorted(itree):
        start_frg_idx = all_intervals.index(start_frg)
        dfs_util(start_frg, all_intervals.copy()[start_frg_idx+1:], longest_path,  start_frg.data.tcov)

    return longest_path.path



def find_busco_category(itree: IntervalTree, min_length_percent:float, max_overlap_percent:float=0) -> int:
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
    interval_0 = list(itree.all_intervals)[0]
    length_threshold = interval_0.data.tlen*min_length_percent
    overlap_threshold = interval_0.data.qlen*max_overlap_percent

    # Case complete single copy
    if len(itree) == 1:
        if interval_0.data.tcov >= length_threshold:
            return 0
        return 3

    longest = max(list(itree.all_intervals), key=lambda i: i.data.tcov)

    # Case either complete single copy or multiple copies
    if longest.data.tcov >=  length_threshold:
        itree.remove(longest)
        if get_longest_path(itree, overlap_threshold) >= length_threshold:
            return 1
        return 0

    # Case check for fragmented
    if get_longest_path(itree, overlap_threshold) >= length_threshold:
        return 2

    return 3


class testHit(NamedTuple):
    qlen: int


def main():
    # itree = IntervalTree(
    #         [Interval(0, 5, testHit(5)),Interval(6, 10, testHit(4)), Interval(11, 20, testHit(9))])
    #
    # print(get_longest_path(itree))
    # return 0
    table_path = Path("data/test_bench/test_out_Hexa_E_1.tsv")
    table_path = Path('/home/yjkbertrand/Documents/projects/Nextpyper/nextpyper/data/test_data/busco_bench/test_out_Hexa_E_1.tsv')
    min_tcov = 0.25

    df = pl.read_csv(table_path, separator="\t", has_header=True)
    print(df.head())
    df2 = (
        df.filter(pl.col("tcov") >= min_tcov)
        .sort(["query", "bits"], descending=True)
        .group_by("query")
        .agg(pl.all().first())
    )  ## Filter by target coverage and Keep only the best hit per query
    target_trees = build_target_trees(df2)
    i = 0
    for target, itree in target_trees.items():
        if i == 15:
            break
        print(target, len(itree))
        print("score", find_busco_category(itree, 0.7, 0))
        i+=1

if __name__ == "__main__":
    main()