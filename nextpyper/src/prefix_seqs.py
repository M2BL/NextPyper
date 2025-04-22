#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2024
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.
"""Prefixes the given string to all the sequences of the given fasta file
"""

__version__ = "0.1"

# =======================================================================================
#               IMPORTS
# =======================================================================================

from pathlib import Path
from itertools import repeat
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys

# =============================================================================
#                FUNCTIONS
# =============================================================================


def prefix_fasta(input: Path, output: Path, p: str) -> None:
    """Prefix all the records in the input fasta with p and write to output."""
    SeqIO.write(map(pref_rec, SeqIO.parse(input, "fasta"), repeat(p)), output, "fasta")


def pref_rec(record: SeqRecord, prefix: str) -> SeqRecord:
    record.id = record.name = prefix + record.name
    record.description = ""
    return record


def pref_path(path: str, p: str) -> str:
    return f"P\t{p}{path[2:]}"


def prefix_graph(graph_lines: list[str], p: str) -> list[str]:
    return [pref_path(line) if line.startswith("P") else line for line in graph_lines]


def prefix_gfa(input: Path, output: Path, p: str) -> None:
    """Prefix all the paths in the input gfa with p and write to output."""

    with Path(input).open() as graph_in, Path(output).open("w") as graph_out:
        new_graph = prefix_graph(graph_in, p)
        graph_out.write("".join(new_graph))


def prefix_hybseq_dir(root_dir: Path, out_dir: Path, sep="_") -> None:
    """Given a root_folder with hybseq samples and assemblies
    Prefix the sequences in their assemblies with sample_probe
    and output them in out_dir

    Note: Temporal solution while implementing Snakemake pipeline
    """

    for file in root_dir.glob("*/*/*/contigs.fasta"):
        sample, probe = file.parts[-4:-2]
        output = out_dir / sample / probe / file.name
        output.parent.mkdir(parents=True, exist_ok=True)
        prefix_fasta(file, output, f"{sample}{sep}{probe}{sep}")


if __name__ == "__main__":
    sys.exit(prefix_hybseq_dir(Path(sys.argv[1]), Path(sys.argv[2])))
