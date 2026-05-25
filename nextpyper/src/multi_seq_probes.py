#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2024
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.


__version__ = "0.1"

# =======================================================================================
#               IMPORTS
# =======================================================================================

from pathlib import Path
from operator import itemgetter
from itertools import groupby
from collections import defaultdict
from functools import partial
from typing import TextIO
import re
import sys

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from snakemake.utils import validate
import polars as pl
from polars.exceptions import NoDataError
import pandas as pd

# =============================================================================
#                FUNCTIONS
# =============================================================================

PROBE_SCHEMA = (Path(__file__).parent / "../workflow/schemas/probes.yaml").resolve()


class NoGrouping(Exception):
    "Exception raised when a pattern (RegEx) does not yield any grouping"


class NoMatch(Exception):
    "Exception raised when a pattern (RegEx) does not match a probre name"


def group_probes(
    recs: list[SeqRecord],
    pattern: str,
    match_group: int | str = 1,
    strict: bool = True,
) -> dict[str, list[SeqRecord]]:
    """Given a list of records with multiple sequences per probe, use the given
    pattern (a RegEx) to group the sequences using their ID. The pattern must
    have at least one capture group. By default, the first capture group will be used
    to group the probes, although other captures group can be specified with match_group.
    If strict, raise a NoGrouping exception if the pattern does not group any records.
    """

    def get_probe_generic(rec: SeqRecord, pattern: re.Pattern, match_group: int | str):
        if (match := pattern.search(rec.id)) is None:
            raise NoMatch(f"Probe {rec.id} does not match pattern {pattern.pattern}")
        else:
            return match[match_group]

    pat = re.compile(pattern)
    if pat.groups == 0:
        raise ValueError(
            "RegEx must have at least one capture group that groups the sequences."
        )
    elif isinstance(match_group, int):
        if match_group > pat.groups:
            raise IndexError(
                f"{match_group=} bigger than the number of capture groups ({pat.groups}) in RegEx."
            )
    elif isinstance(match_group, str):
        if match_group not in pat.groupindex:
            raise ValueError(f"{match_group=} not defined in RegEx {pattern}.")

    get_probe = partial(get_probe_generic, pattern=pat, match_group=match_group)
    probe_recs = {
        probe: list(recs)
        for probe, recs in groupby(sorted(recs, key=get_probe), key=get_probe)
    }

    if len(probe_recs) == len(recs) and strict:
        raise NoGrouping(f"Pattern {pattern} yielded no grouping of the sequences.")

    return probe_recs


def write_summary(probe_counts: dict[str, int], output: Path | TextIO) -> None:
    if isinstance(output, Path):
        with output.open("w") as out:
            for probe, count in probe_counts.items():
                out.write(f"{probe}\t{count}\n")
    else:
        for probe, count in probe_counts.items():
            output.write(f"{probe}\t{count}\n")


def write_hierarchy(
    probe_counts: dict[str, list[str]], output: Path, sep: str = ","
) -> None:
    with output.open("w") as out:
        for probe, names in probe_counts.items():
            out.write(f"{probe}\t{sep.join(names)}\n")


def check_probes(
    probes_path: Path,
    pattern: str,
    out_summary: Path | None,
    out_hierarchy: Path | None,
):
    """Validate that the given probes file follows the naming convention of the schema
    and that the given pattern does group the probes.
    """  # ToDo: Write proper docstring

    probes = list(SeqIO.parse(probes_path, "fasta"))
    probe_ids = [rec.id for rec in probes]

    # Test if the probe names comply the naming convention
    validate(pd.DataFrame({"probe_name": probe_ids}), schema=PROBE_SCHEMA)

    print("Probe sequence names comply with naming convention.", file=sys.stderr)
    try:
        probes_dict = group_probes(probes, pattern)
    except NoMatch as err:
        print(err, file=sys.stderr)
        print(
            "At least one of probe does not match the given pattern",
            file=sys.stderr,
        )
        return
    except NoGrouping as err:
        print(err, file=sys.stderr)
        print(
            "Either this is a single-probe set or the pattern is not appropiate for multi-probe mode",
            file=sys.stderr,
        )
        return

    probes_counts = {probe: len(recs) for probe, recs in probes_dict.items()}

    print(
        f"The pattern yields {len(probes_dict)} groups from {len(probes)} probe sequences",
        file=sys.stderr,
    )
    # print("The probes have the following number of members: \n", file=sys.stderr)
    # write_summary(probes_counts, sys.stdout)

    if out_summary:
        write_summary(probes_counts, out_summary)

    if out_hierarchy:
        probes_hier = {
            probe: [rec.id for rec in recs] for probe, recs in probes_dict.items()
        }
        write_hierarchy(probes_hier, out_hierarchy)


def add_chim_tag(
    rec: SeqRecord,
    chimera_set: set[str],
    tribble_set: set[tuple[str]],
    pat: re.Pattern,
):
    """Add chimera tag to the record based on the Vsearch de novo detection and tribble provenance.
    If the record comes from a tribble it is tagged as 'putative'. If it is additionally detected as
    chimeric by vsearch it is tagged as 'chimeric'. Otherwise, the record is considered normal.
    """

    get_comp = itemgetter("probe", "sample2", "seed_id", "comp")

    in_tribble = get_comp(pat.match(rec.id)) in tribble_set
    in_chimeras = rec.id in chimera_set

    if in_tribble and in_chimeras:
        tag = "chimeric"
    elif in_tribble:
        tag = "putative"
    else:
        tag = "normal"

    rec.description += f" [chimeric={tag}]"


# Main execution for rule "per_probe_scaffold_grouping" in workflow
def snakemake_call(snakemake):
    with open(snakemake.log[0], "w") as outlog:
        sys.stdout = sys.stderr = outlog

        inputs = snakemake.input
        outfolder = Path(snakemake.output[0]).parent
        pattern = snakemake.params.pattern
        probes_list = snakemake.params.get("probes")
        mode = snakemake.params.mode
        pat = re.compile(pattern, re.VERBOSE)

        # Todo: Refactor this block as a function
        match mode:
            case "exons" | "genetigs" | "supercontigs":
                all_recs = [
                    rec
                    for file in Path(inputs.scfs[0]).parent.rglob(f"*_{mode}.fasta")
                    for rec in SeqIO.parse(file, "fasta")
                ]
                grouped_recs = group_probes(all_recs, pat, match_group="sample1")

                ## Chimera tagging
                tag_sets = defaultdict(set)
                tribble_sets = defaultdict(set)

                # Load de novo chimera detection results by vsearch
                # if inputs.get("chimera_tags"):
                #     for table in map(Path, inputs.chimera_tags):
                #         try:
                #             tag_sets[table.stem].update(
                #                 pl.read_csv(table, separator="\t", has_header=False)[
                #                     "column_2"
                #                 ]
                #             )
                #         except NoDataError:
                #             pass

                if inputs.get("tribbles"):
                    # Load tribble tables
                    tribble_sets = {
                        tribble.stem: set(
                            pl.read_csv(tribble, separator="\t").iter_rows()
                        )
                        for tribble in map(Path, inputs.tribbles)
                    }

                # Add chimera tags to the sequences
                if inputs.get("tribbles") or inputs.get("chimera_tags"):
                    for sample, recs in grouped_recs.items():
                        for rec in recs:
                            add_chim_tag(
                                rec, tag_sets[sample], tribble_sets[sample], pat
                            )

            case "scfs":
                all_recs = [
                    rec for file in inputs for rec in SeqIO.parse(file, "fasta")
                ]
                grouped_recs = group_probes(all_recs, pat, match_group="probe")
            case "single_probes":
                grouped_recs = {
                    pat.search(rec.id)[1]: rec
                    for rec in SeqIO.parse(inputs.probes, "fasta")
                }
            case "multi_probes":
                tables_dir = Path(inputs.tables[0]).parent
                matched_probes = set(
                    pd.concat(
                        pd.read_csv(table, sep="\t") for table in tables_dir.iterdir()
                    )["theader"].unique()
                )
                all_recs = [
                    rec
                    for rec in SeqIO.parse(inputs.probes, "fasta")
                    if rec.id in matched_probes
                ]
                grouped_recs = group_probes(all_recs, pat)
                # It could happen that in a very ideal scenario where a single probe
                # version per probe is the only surviving probe, this would raise a NoGrouping exception.
            case _:
                raise ValueError(
                    f"{mode=} not recognized. Use scfs, supercontigs, single_probes or multi_probes."
                )

        # Writing output is common to all uses
        # In "supercontigs" mode, "probe" is actually samples.
        outfolder.mkdir(exist_ok=True)
        for probe, recs in grouped_recs.items():
            SeqIO.write(recs, outfolder / f"{probe}.fasta", "fasta")

        # In all cases, we may have missing probes, so we need to at least touch them
        # In supercontigs mode this is not needed, because the outputs are per sample,
        # which are guaranteed to exist.
        if mode in ("scfs", "single_probes", "multi_probes"):
            for probe in probes_list:
                (outfolder / f"{probe}.fasta").touch(exist_ok=True)


if __name__ == "__main__":
    if "snakemake" in globals():
        snakemake_call(snakemake)
