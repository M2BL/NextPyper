#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2024
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.
"""Fix Saute sequences headers by updating the length and cov values and add the sample name"""

__version__ = "0.1"

# =======================================================================================
#               IMPORTS
# =======================================================================================

from pathlib import Path
from itertools import pairwise
import re
import sys

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from more_itertools import interleave_longest


# =============================================================================
#                FUNCTIONS
# =============================================================================

SAUTE_PAT = r"^(?P<sample>.*?)-(?P<probe>.*?)_EDGE_(?P<seed_id>\d+)_length_(?P<len>\d+)_cov_(?P<cov>[\w.]+):[^ ]+:(?P<kmers>\d+)$"


def update_specs(rec: SeqRecord, pat: re.Pattern) -> SeqRecord:
    """Mach the record name with the provided pattern and replace
    the length and cov matching groups.

    The pattern must have len, cov and kmers as capture groups.
    """

    if (match := pat.search(rec.id)) is None:
        raise ValueError(f"{rec.id} does not match {pat.pattern}")

    info = match.groupdict()
    info["len"] = str(len(rec))
    info["cov"] = str(round(int(info["kmers"]) / len(rec), 4))

    prefix = rec.id[0 : match.regs[1][0]]
    suffix = rec.id[match.regs[-1][1] : len(rec.id)]

    gaps = [rec.id[a:b] for (_, a), (b, _) in pairwise(match.regs[1:])]
    name = prefix + "".join(interleave_longest(info.values(), gaps)) + suffix

    rec.name = rec.description = ""
    rec.id = name
    return rec


def fix_saute_headers(
    recs: list[SeqRecord], sample: str, pattern: str
) -> list[SeqRecord]:
    """Given a Saute assembly, fix the sequence headers by prepending the
    provided sample to the name, updating the specs (length, cov) and
    removing the "Contig_" at the start.
    """

    pat = re.compile(pattern, re.VERBOSE)
    new_recs = []
    for rec in recs:
        new_rec = update_specs(rec, pat)
        new_rec.id = f"{sample}|{new_rec.id.removeprefix("Contig_")}"
        new_recs.append(new_rec)

    return new_recs


def snakemake_call(snakemake):
    records_path = Path(snakemake.input[0])
    out_path = Path(snakemake.output[0])
    sample = snakemake.params.sample
    pattern = snakemake.params.get("pattern", SAUTE_PAT)

    recs = SeqIO.parse(records_path, "fasta")
    new_recs = fix_saute_headers(recs, sample, pattern)
    SeqIO.write(new_recs, out_path, "fasta")


def main():
    if len(sys.argv) != 4:
        print("Usage: python fix_headers.py <saute_asm.fasta> <output.fasta> <sample>")
        sys.exit(1)

    records_path = Path(sys.argv[1])
    out_path = Path(sys.argv[2])
    sample = sys.argv[3]

    recs = SeqIO.parse(records_path, "fasta")
    new_recs = fix_saute_headers(recs, sample, SAUTE_PAT)
    SeqIO.write(new_recs, out_path, "fasta")


if __name__ == "__main__":
    if "snakemake" in globals():
        snakemake_call(snakemake)
    else:
        main()
