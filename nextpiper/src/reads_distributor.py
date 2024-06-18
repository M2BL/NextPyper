#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2024
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.
"""Functions to assign and distribute mapped reads to the corresponding
probes where they mapped."""

__version__ = "0.1"

# =============================================================================
#               IMPORTS
# =============================================================================

from pathlib import Path
from operator import attrgetter
from collections import defaultdict
from itertools import groupby
import pysam

# =============================================================================
#                FUNCTIONS
# =============================================================================


same_map = lambda aln: aln.reference_name == aln.next_reference_name
which_mapq = lambda aln1, aln2: (
    aln1 if aln1.mapping_quality > aln2.mapping_quality else aln2
)


def assign_alns(alns: tuple[pysam.AlignedSegment]) -> str | None:
    if len(alns) > 2:
        raise ValueError(f"Multi mapping not implemented yet {alns}")
    else:
        match (alns[0].is_unmapped, alns[1].is_unmapped):
            case (True, True):
                return None
            case (False, True):
                return alns[0].reference_name
            case (True, False):
                return alns[1].reference_name
            case (False, False):
                chosen = alns[0] if same_map(alns[0]) else which_mapq(*alns)
                return chosen.reference_name


def distribute_reads(inbam: Path, outdir: Path):
    """Given a input bam distribute the read pairs into the different
    targets they mapped to."""

    outdir = Path(outdir)
    handle = pysam.AlignmentFile(inbam, "rb")
    assigned_dict = defaultdict(list)

    # Assign read pairs to a target/probe (best mapping)
    for _, g_alns in groupby(
        (aln for aln in handle if not aln.is_supplementary),
        key=attrgetter("query_name"),
    ):
        alns = tuple(g_alns)
        if (ref := assign_alns(alns)) is not None:
            assigned_dict[ref].extend(alns)

    # Write uncompressed bams for each target/probe in the outdir
    outdir.mkdir(exist_ok=True, parents=True)
    for probe, alns in assigned_dict.items():
        outdir / f"{probe}.bam"
        with pysam.AlignmentFile(
            outdir / f"{probe}.bam", "wbu", template=handle
        ) as outbam:
            for aln in alns:
                outbam.write(aln)


if __name__ == "__main__":
    # Snakemake rule execution by the "script:" directive
    if "snakemake" in globals():
        distribute_reads(str(snakemake.input), str(snakemake.output))
