#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2024
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.
# Cleaning alignments with TrimAl and TAPER (https://github.com/chaoszhang/TAPER)
#  pyjulia issue:Your Python interpreter "/usr/bin/python3.7"
# is statically linked to libpython.  Currently, PyJulia does not fully
# support such Python interpreter.
# solutions: https://stackoverflow.com/questions/64486932/compile-and-use-custom-system-image-for-pyjulia
# https://github.com/JuliaPy/pyjulia/issues/310
# or solve it in snakemake

__version__ = "0.1"

# =======================================================================================
#               IMPORTS
# =======================================================================================

import os
from pathlib import Path
import subprocess
import shutil


# =============================================================================
#                FUNCTIONS
# =============================================================================
def run_subprocess(cmd: str):
    try:
        process = subprocess.run(
            cmd.split(),
            timeout=100,
            check=True,
            capture_output=False,
            encoding="utf-8",
        )
    except FileNotFoundError as exc:
        print(f"Process failed because the executable could not be found.\n{exc}")
        raise
    except subprocess.CalledProcessError as exc:
        print(
            f"Process failed because did not return a successful return code. "
            f"Returned {exc.returncode}\n{exc}"
        )
        raise
    except subprocess.TimeoutExpired as exc:
        print(f"Process timed out.\n{exc}")
        raise


def align_clean(
    alignment: Path,
    output_folder: Path,
    taper_parameters: Path,
    nbr_repetitions: int = 2,
    trimal_gt: float = 0.2,
):
    run_nbr = 0

    if not (Path.cwd() / "temp").exists():
        os.mkdir(Path.cwd() / "temp")
    shutil.copy(alignment, Path.cwd() / "temp")
    os.chdir(Path.cwd() / "temp")

    while run_nbr < nbr_repetitions:
        suffix = alignment.stem
        run_nbr += 1
        taper_cmd = f"julia /home/yjkbertrand/Documents/projects/radseq/bin/TAPER/correction_multi.jl -m N -a N -c 1 {alignment} > taper_{run_nbr}.fasta"
        trimal_cmd = f"trimal -gt {trimal_gt} -in taper_{run_nbr}.fasta  -out {suffix}_{run_nbr}.fasta  "
        os.system(taper_cmd)
        run_subprocess(trimal_cmd)


def main(): ...


if __name__ == "__main__":
    os.chdir(
        "/home/yjkbertrand/Documents/projects/nextpyper/test_data/batrachium/exonerate/clusters_2"
    )
    alignment = Path.cwd() / "probe_3_aa_1.fasta"
    taper_parameters = Path.cwd() / "taper_parameters.txt"
    align_clean(alignment, Path.cwd(), taper_parameters)
    # main()
