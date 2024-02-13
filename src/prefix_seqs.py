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
