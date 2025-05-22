from pathlib import Path
import polars as pl
from intervaltree import Interval, IntervalTree
from itertools import groupby
from operator import itemgetter
from typing import NamedTuple


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
    with similarity below min_idt are ignored."""

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


def main():

    table_path = Path("data/test_bench/test_out_Hexa_E_1.tsv")
    min_tcov = 0.25

    df = pl.read_csv(table_path, separator="\t", has_header=True)
    df2 = (
        df.filter(pl.col("tcov") >= min_tcov)
        .sort(["query", "bits"], descending=True)
        .group_by("query")
        .agg(pl.all().first())
    )  ## Filter by target coverage and Keep only the best hit per query
    target_trees = build_target_trees(df2)
