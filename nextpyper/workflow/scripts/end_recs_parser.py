import argparse
from pathlib import Path
from collections import defaultdict
import re

from Bio import SeqIO

"""This script parses the final output of different pipelines (NextPyper, Captus, HybPiper) to restructure
the results from being organized per probe to being organized per sample."""


CAPTUS_HEADER = r"\[query=(.*)\]\s\[hit=(\d+)\]"


def parse_nextpyper(results_path: Path, parsed_dir: Path, kind: str = "supercontigs"):
    """
    Parse the end recs from the NextPyper results_path (separation_output folder) structured
    per probe, restructure them per sample and save them in parsed_dir, in a file per sample.
    """

    parsed_dir.mkdir(parents=True, exist_ok=True)
    if not results_path.exists():
        raise FileNotFoundError(f"Results path {results_path} does not exist.")
    if not results_path.is_dir():
        raise NotADirectoryError(f"Results path {results_path} is not a directory.")

    sample_recs = defaultdict(list)
    for file in results_path.rglob(f"*{kind}.fasta"):
        for rec in SeqIO.parse(file, "fasta"):
            sample = rec.id.split("|")[0]
            sample_recs[sample].append(rec)

    for sample, recs in sample_recs.items():
        SeqIO.write(recs, parsed_dir / f"{sample}.fasta", "fasta")


def parse_captus(results_path: Path, parsed_dir: Path):
    """
    Parse the end recs from the captus results_path structured per probe,
    restructure them per sample and save them in parsed_dir, in a file per sample.
    Sequences are uppercased and dashes are removed.

    The probe and hit fields in the description are added as a suffix to the sequence
    name to eliminate name duplicates.
    """

    head_pat = re.compile(CAPTUS_HEADER)

    sample_recs = defaultdict(list)
    for file in results_path.glob("*.fna"):
        for rec in SeqIO.parse(file, "fasta"):
            sample = rec.id.split("__")[0]
            rec.description = rec.description.removeprefix(f"{rec.id} ")
            match = head_pat.search(rec.description)
            rec.id = rec.name = f"{rec.id}|{match.group(1)}_hit{match.group(2)}"
            rec.seq = rec.seq.upper().replace("-", "")
            sample_recs[sample].append(rec)

    for sample, recs in sample_recs.items():
        SeqIO.write(recs, parsed_dir / f"{sample}.fasta", "fasta")


def parse_hybpiper(results_path: Path, parsed_dir: Path):
    """
    Parse the end recs from the HybPiper results_path where two folders are expected paralogs
    and supercontigs, structured per probe. The sequences are restructured per sample. For a
    given probe, if there are more sequences in the paralogs, the paralogs are used, otherwise
    the supercontigs are used. The sequences are then saved in parsed_dir, in a file per sample.

    The file name (expected to be the probe) is added as a suffix to the sequence name to
    eliminate name duplicates.
    """

    sample_recs = defaultdict(list)
    for sctg_path in (results_path / "supercontigs").glob("*.fasta"):
        paralogs_path = results_path / f"paralogs/{sctg_path.stem}_paralogs_all.fasta"

        sctg = list(SeqIO.parse(sctg_path, "fasta"))
        paralogs = list(SeqIO.parse(paralogs_path, "fasta"))

        recs = paralogs if len(paralogs) > len(sctg) else sctg
        for rec in recs:
            sample = rec.id.split(".")[0]
            rec.id = rec.name = f"{rec.id}|{sctg_path.stem}"
            rec.seq = rec.seq.upper().replace("-", "").replace("N", "")
            sample_recs[sample].append(rec)

    for sample, recs in sample_recs.items():
        SeqIO.write(recs, parsed_dir / f"{sample}.fasta", "fasta")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Parse NextPyper final output structure from per probe to per sample."
    )
    parser.add_argument(
        "pipeline_results", type=Path, help="Path to per probe sequence output."
    )
    parser.add_argument(
        "output_folder", type=Path, help="Path where to write the parsed results."
    )
    parser.add_argument(
        "pipeline",
        type=str,
        choices=["nextpyper", "captus", "hybpiper"],
        help="Pipelined used (affects the expected file structure).",
    )

    parser.add_argument(
        "-k",
        "--kind",
        type=str,
        default="supercontigs",
        choices=["supercontigs", "genotigs", "exons"],
        help="For NextPyper, which kind of sequences to retrieve.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    args.output_folder.mkdir(parents=True, exist_ok=True)
    if not args.pipeline_results.exists():
        raise FileNotFoundError(f"Results path {args.pipeline_results} does not exist.")
    if not args.pipeline_results.is_dir():
        raise NotADirectoryError(
            f"Results path {args.pipeline_results} is not a directory."
        )
    if not args.output_folder.is_dir():
        raise NotADirectoryError(
            f"Parsed directory {args.output_folder} is not a directory."
        )

    match args.pipeline:
        case "nextpyper":
            parse_nextpyper(args.pipeline_results, args.output_folder, args.kind)
        case "captus":
            parse_captus(args.pipeline_results, args.output_folder)
        case "hybpiper":
            parse_hybpiper(args.pipeline_results, args.output_folder)
        case _:
            raise ValueError(f"Unknown pipeline: {args.pipeline}")


if __name__ == "__main__":
    main()
