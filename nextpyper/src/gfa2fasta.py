#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2024
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.

"""Simple script to retrieve the paths of an assembly graph and output it as sequences in
fasta format. The functions are tested to work for standard De Brujn graphs, and assume a
constant overlap between segments. Heterogeneous overlaps are supported in J-lines.
Blunted graphs (e.g. variation graphs) should work too, but have not been tested.

#Usage example:
    gfa_file = Path("path/to/my_graph.gfa")
    out_fasta = Path("my_path_sequences.fasta")
    paths2fasta(gfa_file, out_fasta)

The headers of the sequences in the fasta files correspond to the names in the P lines.
"""

__version__ = "0.1"

# =======================================================================================
#               IMPORTS
# =======================================================================================

from pathlib import Path
import sys

from Bio import SeqIO
from gfa_graph import Assembly_graph


# =============================================================================
#                FUNCTIONS
# =============================================================================


def paths2fasta(graph_path: Path, output_path: Path) -> None:
    """Retrieve the paths from the assembly graph and write them to
    the specified output file in fasta format."""

    graph = Assembly_graph(graph_path)
    recs = (graph.retrieve_path(path) for path in graph.paths.values())
    SeqIO.write(recs, output_path, "fasta")


def snakemake_call(snakemake):
    graph_path = Path(snakemake.input[0])
    output_path = Path(snakemake.output[0])
    paths2fasta(graph_path, output_path)


def main():
    if len(sys.argv) != 3:
        print("Usage: python gfa2fasta.py <graph.gfa> <output.fasta>")
        sys.exit(1)

    graph_path = Path(sys.argv[1])
    output_path = Path(sys.argv[2])
    paths2fasta(graph_path, output_path)


if __name__ == "__main__":
    if "snakemake" in globals():
        snakemake_call(snakemake)
    else:
        main()
