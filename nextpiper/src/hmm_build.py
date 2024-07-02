#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2024
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.

"""
Function used to build HMM profiles from alignment and generate consensus sequences from profiles
"""
__version__ = "0.1"

# =======================================================================================
#               IMPORTS
# =======================================================================================

from pathlib import Path
import pyhmmer
from pyhmmer.easel import MSAFile
from pyhmmer.plan7 import HMMFile

alphabet = pyhmmer.easel.Alphabet.amino()

# =============================================================================
#                FUNCTIONS
# =============================================================================


def hmm_build(fasta_path: Path, hmm_output_path: Path) -> None:
    """
    Build HMM profile from alignment and save it in hmm format.
    :param fasta_path:
    :param hmm_output_path:
    :return: save the hmm profile into a file.
    """
    with MSAFile(fasta_path, digital=True, alphabet=alphabet, format="afa") as msa_file:
        msa = msa_file.read()
        msa.name = str.encode(fasta_path.name.removesuffix(".fasta"))
        builder = pyhmmer.plan7.Builder(alphabet)
        background = pyhmmer.plan7.Background(alphabet)
        hmm, _, _ = builder.build_msa(msa, background)
        with open(hmm_output_path, "wb") as output_file:
            hmm.write(output_file)


def hmm_consensus(hmm_path: Path, output_path: Path) -> None:
    with HMMFile(hmm_path) as hmm_file:
        while True:
            try:
                # iterate over all hmms in profile
                hmm = next(hmm_file)
                seq1 = pyhmmer.easel.TextSequence(name=hmm.name, sequence=hmm.consensus)
                with output_path.open("wb") as f:
                    seq1.write(f)
            except StopIteration:
                break


if __name__ == "__main__":
    ...
