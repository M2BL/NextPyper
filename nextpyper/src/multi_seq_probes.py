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
from itertools import groupby
from functools import partial
import re
import sys

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from snakemake.utils import validate
import pandas as pd

# =============================================================================
#                FUNCTIONS
# =============================================================================

PROBE_SCHEMA = (Path(__file__).parent / "../schemes/probes.yaml").resolve()


class NoGrouping(Exception):
    "Exception raised when a pattern (RegEx) does not yield any grouping"


class NoMatch(Exception):
    "Exception raised when a pattern (RegEx) does not match a probre name"


def group_probes(
    recs: list[SeqRecord], pattern: str, match_group: int | str = 1
) -> dict[str, list[SeqRecord]]:
    """Given a list of records with multiple sequences per probe, use the given
    pattern (a RegEx) to group the sequences using their ID. The pattern must
    have at least one capture group. The first capture group will be used
    to group the probes
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

    if len(probe_recs) == len(recs):
        raise NoGrouping(f"Pattern {pattern} yielded no grouping of the sequences.")

    return probe_recs


def write_summary(probe_counts: dict[str, int], output: Path) -> None:
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
            "Either this is a single-probe set or the pattern is not appropiate for the multi-probe",
            file=sys.stderr,
        )
        return

    probes_counts = {probe: len(recs) for probe, recs in probes_dict.items()}

    print("The probes have the following number of members: \n", file=sys.stderr)
    write_summary(probes_counts, sys.stdout)

    if out_summary:
        write_summary(probes_counts, out_summary)

    if out_hierarchy:
        probes_hier = {
            probe: [rec.id for rec in recs] for probe, recs in probes_dict.items()
        }
        write_hierarchy(probes_hier, out_hierarchy)


# Main execution for rule "per_probe_scaffold_grouping" in workflow
def snakemake_call(snakemake):
    with open(snakemake.log[0], "w") as outlog:
        sys.stdout = sys.stderr = outlog

        inputs = snakemake.input
        outfolder = Path(snakemake.output[0]).parent
        pattern = snakemake.params.pattern
        probes_list = snakemake.params.probes
        pat = re.compile(pattern, re.VERBOSE)

        all_recs = [rec for file in inputs for rec in SeqIO.parse(file, "fasta")]
        grouped_scfs = group_probes(all_recs, pat, match_group="probe")

        outfolder.mkdir(exist_ok=True)
        for probe, recs in grouped_scfs.items():
            SeqIO.write(recs, outfolder / f"{probe}.fasta", "fasta")

        for probe in probes_list:
            (outfolder / f"{probe}.fasta").touch(exist_ok=True)


if __name__ == "__main__":
    if "snakemake" in globals():
        snakemake_call(snakemake)
