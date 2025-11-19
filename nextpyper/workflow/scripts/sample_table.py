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
import os

# =============================================================================
#                FUNCTIONS
# =============================================================================

COLS = ["sample", "path_forward", "path_reverse", "ploidy"]


# =============================================================================
#                FUNCTIONS
# =============================================================================


def make_table(data_dir: Path, out_table: Path, sep: str = "\t") -> None:
    """Parse the directory with paired-end reads data and write a sample table
    specifying paths for each sample. The sample name is inferred from the common
    prefix of each pair of reads.

    The output table has 3 columns: sample_name, path_forward_reads, path_reverse_reads.
    If extra is given, a fourth column will be included for all the samples with extra.
    """

    with out_table.open("w") as table:
        files = sorted(data_dir.iterdir())

        if len(files) % 2 != 0:
            raise ValueError(
                f"There is an odd number of files in {data_dir}. All samples must have forward and reverse files."
            )

        # Add the header to the table
        table.write(f"{sep.join(COLS)}\n")

        for file1, file2 in zip(files[::2], files[1::2]):
            # Get sample name
            sample = os.path.commonprefix([file1.name, file2.name])

            ## Do name processing magic
            if sample.endswith("R"):
                sample = sample.removesuffix("R")
            sample = sample.rstrip("._-")
            sample = sample.replace("-", "_")

            row = f"{sample}{sep}{file1.resolve()}{sep}{file2.resolve()}{sep}0\n"
            table.write(row)
