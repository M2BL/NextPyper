#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2024
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.
"""
Functions and classes for running miniprot in different steps of the pipeline.

The class ComponentFilter is used after the SAUTE assembly to discard the scaffolds that do not match
the selected probes. At this step, the probe set should have been downsized as the complexity of the search is
quadratic. Filtered scaffolds are then used for a vsearch clustering. The scaffolds file can be per sample or per probe
per sample. Adjust the probe file accordingly.
#  Usage example:
    probe_fasta = ".../test_data/test_clustering/probes_aa.fasta"
    scaffolds_fasta = "../test_data/test_clustering/sp_1.fasta
    # load the data, perform the computation:
    cpf = ComponentFilter(probe_fasta, scaffolds_fasta)
    # save the results in fasta format:
    cpf.save("../test_data/test_clustering/sp_1_filtered.fasta")

The class OverlappingCds is used after the vsearch clustering of SAUTE scaffolds and consensus estimation.
It maps all probe versions against the consensus sequences of vsearch clusters. It seeks
to find the best mapping probe version over all the consensus sequences and use it as
a common reference in order to determine which sequences are not overlapping.
Strictly non-overlapping sequences (there are no bridging scaffolds) are saved in separate
files.
#  Usage example:
    probe_fasta = ".../test_data/test_clustering/probe_3_aa.fasta"
    consensus_contig_fasta = "../test_data/test_clustering/gene_3_consensus.fasta"
    # load the data, perform the computation:
    olc = OverlappingCds(probe_fasta, consensus_contig_fasta)
    # specify a folder where each non-overlapping group of sequences is saved in fasta format.
    olc.save_scaffolds("../test_data/test_clustering/non_overlapping/")
    # save the sequence of the best overall mapping probe.
    olc.save_best_probe("../test_data/test_clustering/gene_3_best_probe.fasta")
"""

__version__ = "0.1"

# =======================================================================================
#               IMPORTS
# =======================================================================================

from dataclasses import dataclass, field
from collections import defaultdict
from importlib import reload
from io import StringIO
from operator import attrgetter
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import textwrap
from typing import Optional, Self

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# import gff_parser
# gff_parser = reload(gff_parser)
from gff_parser import Fragment, Cds
from interval_tree import IntervalST, Interval

# Parse the header of SAUTE scaffolds
saute_pattern = re.compile(
    r"^Contig_(?P<name>.*?):(?P<component>\d+?):[^ ]+$",
    re.VERBOSE,
)

# =======================================================================================
#               EXCEPTIOMS
# =======================================================================================


# =======================================================================================
#               FUNCTIONS
# =======================================================================================


def get_nuc_coordinates(
    probe_start: int, probe_end: int, fragments: list[Fragment]
) -> list[Optional[int]]:
    """
    Given a start and end position on a probe in AA space, return the nucleotide coordinates of start and end.
    :param probe_start:
    :param probe_end:
    :param fragments: list of fragments
    :return: [None, None] if for some reason there are no AA-nt correspondence, otherwise [start, end]
    """
    correspondences = {}
    for fragment in fragments:
        correspondences.update(fragment.get_correspondence())
    while probe_start < probe_end:
        if (first := correspondences.get(probe_start)) is None:
            probe_start += 1
        if (last := correspondences.get(probe_end)) is None:
            probe_end -= 1
        if None not in [first, last]:
            return sorted(
                [
                    first,
                    last,
                ]
            )
    return [None, None]


def merge_intervals(arr: list["OverlappingSeqs"]):
    """
    Combine overlapping intervals of sequences and cluster their names.
    """
    # Sorting based on the increasing order of the start intervals
    arr.sort(key=attrgetter("probe_start"))
    index = 0
    index_to_del = []
    for i in range(1, len(arr)):
        if arr[index].probe_end > arr[i].probe_start:
            arr[index].probe_end = max(arr[index].probe_end, arr[i].probe_end)
            arr[index].combine(arr[i].seq_names)
            index_to_del.append(i)
        else:
            index = index + 1
            arr[index] = arr[i]
    return [arr[index] for index in range(len(arr)) if index not in index_to_del]


def run_miniprot(
    probe_path: Path,
    contig_path: Path,
    treads=2,
    min_similarity=0.85,
    min_coverage=0.01,
) -> StringIO:
    """
    Wrapper for running miniprot.
    :return:
    """

    miniprot_cmd = f"miniprot -t {treads} --gff --aln --outn 1 --outs {min_similarity} --outc {min_coverage} {contig_path} {probe_path}".split()

    try:
        miniprot = subprocess.run(
            miniprot_cmd,
            timeout=1000,
            check=True,
            capture_output=True,
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
    return StringIO(miniprot.stdout)


# =======================================================================================
#               CLASSES
# =======================================================================================


@dataclass
class RankCoverage:
    """
    Container for keeping track of the rank of a given probe out of all probes ranks.
    Attributes
    ----------
    -rank_score: list of ranks. Each int represents how a well given probe aligns to a scaffolds
        compared to the other probes. Smaller ranks denote better mapping.
    -coverage: list of scaffolds to which this probe aligns.
    """

    rank_score: list[int] = field(default_factory=list)
    coverage: list[str] = field(default_factory=list)

    def append_to(self, item_0, item_1):
        self.rank_score.append(item_0)
        self.coverage.append(item_1)

    def get_total_score(self):
        return sum(self.rank_score) / len(self.coverage)


@dataclass
class OverlappingSeqs:
    """
    Container for groups of sequences that overlap a common region of the probe.
    Attributes
    ----------
    -probe_start: the smallest index on the probe (AA space) that has a matching sequence.
    -probe_end: the largest index on the probe (AA space) that has a matching sequence.
    -seq_names: a list of sequence names that overlap this region.
    """

    probe_start: int
    probe_end: int
    seq_names: list[str]

    def combine(self, other: list[str]) -> Self:
        """
        Add items to seq_names
        """
        self.seq_names.extend(other)
        return self

    def get_seq_names(self) -> list[str]:
        return self.seq_names

    def get_probe_interval(self) -> str:
        """
        Get the boundaries of the probe covered by the group of sequences.
        """
        return f"start_{self.probe_start}_end_{self.probe_end}"


@dataclass
class MiniprotInit:
    """
    Base class for miniprot initialization.
    Attributes
    ----------
    -probes_fasta: path to the amino acid probes fasta file (multiprobe)
    -contigs_fasta: path to the nucleotide contigs fasta file.
    -treads: for miniprot.
    -min_probe_contig_sim: mapping similarity for miniprot.
    -min_fragment_cov: fraction of the probe covered by a contig for miniprot.
    -min_contig_length: minimum length of a contig after trimming, if the sequences are saved.

    Post Init
    -probes_path: probe_fasta string converted to Path.
    -contigs_path: contigs_fasta string converted to Path.
    -probes_dict: dict of probe name as key and SeqRecord as value.
    -contigs_dict: dict of contig name as key and SeqRecord as value.

    """

    probes_fasta: str
    contigs_fasta: str
    treads: int = field(default=8)
    min_probe_contig_sim: float = field(default=0.85)
    min_fragment_cov: float = field(default=0.05)
    min_contig_length: int = field(default=300)
    probes_path: Path = field(init=False, repr=False)
    contigs_path: Path = field(init=False, repr=False)
    probes_dict: dict[str, SeqRecord] = field(
        init=False, repr=False, default_factory=dict
    )
    contigs_dict: dict[str, SeqRecord] = field(
        init=False, repr=False, default_factory=dict
    )

    def __post_init__(self):
        self.probes_path = Path(self.probes_fasta)
        self.contigs_path = Path(self.contigs_fasta)
        assert self.probes_path.exists(), f"{self.probes_fasta} does not exist"
        assert self.contigs_path.exists(), f"{self.contigs_fasta} does not exist"
        self._parse_fasta()

    def _parse_fasta(self) -> None:
        """
        Load the probe fasta file and the contig fasta file.
        :return:
        """
        try:
            self.probes_dict = SeqIO.to_dict(SeqIO.parse(self.probes_path, "fasta"))
        except Exception as err:
            sys.exit(f"[ERROR] {err}")
        try:
            self.contigs_dict = SeqIO.to_dict(SeqIO.parse(self.contigs_path, "fasta"))
        except Exception as err:
            sys.exit(f"[ERROR] {err}")

    def _run_miniprot(self, record, probe_path_temp, contig_path_temp) -> Cds:
        SeqIO.write(record, contig_path_temp, "fasta")
        miniprot_out = run_miniprot(
            probe_path_temp,
            contig_path_temp,
            self.treads,
            self.min_probe_contig_sim,
            self.min_fragment_cov,
        )
        return Cds(miniprot_out)


@dataclass
class ComponentFilter(MiniprotInit):
    """
    Scaffolds generated by Saute are filtered for mapping quality. The filtering is performed at the component level.
    Probes are mapped to scaffolds with miniprot, a valid hit per component is sufficient.
    Poorly aligned contigs are filtered out.
    We don't search for the best probe, any hit is sufficient for re-orienting the probe.
    Attributes
    ----------
    -min_global_sim:
    Post Init
    -cds_dict: dict of contig name as key and Cds object as value.
    -used_probes: search first the list of probe variants that have already matched a scaffold.
    """

    min_global_sim: float = field(default=0.80)
    filtered_scaffolds: list[str] = field(init=False, repr=False, default_factory=list)

    def __post_init__(self):
        super().__post_init__()

        # run miniprot on each scaffold/probe combination.
        # When a hit is found, add to filtered_scaffolds.
        print("Running miniprot")
        used_probes: list[str] = []
        valid_components: set[str] = set()
        rejected_components: set[str] = set()
        length_sorted = sorted(
            [(ctg, rec) for ctg, rec in self.contigs_dict.items()],
            key=lambda x: len(x[1].seq),
            reverse=True,
        )
        with tempfile.TemporaryDirectory() as tmpdirname:
            tmp_probe_paths = {}
            for probe_name, record in list(self.probes_dict.items())[:]:
                fasta_name = f"{probe_name}.fas"
                contig_fasta = Path(tmpdirname) / fasta_name
                SeqIO.write(record, contig_fasta, "fasta")
                tmp_probe_paths[probe_name] = contig_fasta

            for contig, record in length_sorted:
                print(f"working on contig {contig}")
                #  Check whether this component has been analyzed previously
                if (match := saute_pattern.match(contig)) is not None:
                    component = f"{match.group('name')}_{match.group('component')}"
                else:
                    sys.exit(f"[ERROR] {contig} does not match saute  header pattern")
                if component in valid_components:
                    self.filtered_scaffolds.append(contig)
                    continue
                if component in rejected_components:
                    continue

                if used_probes:
                    for probe_name in used_probes:
                        probe_path = tmp_probe_paths[probe_name]
                        contig_path = Path(tmpdirname) / f"{contig}.fas"
                        cds = self._run_miniprot(record, probe_path, contig_path)
                        if (
                            not cds.is_empty()
                            and cds.get_global_sim() > self.min_global_sim
                        ):
                            self.filtered_scaffolds.append(contig)
                            valid_components.add(component)
                            break
                    else:
                        for probe_name in set(self.probes_dict.keys()) - set(
                            used_probes
                        ):
                            probe_path = tmp_probe_paths[probe_name]
                            contig_path = Path(tmpdirname) / f"{contig}.fas"
                            cds = self._run_miniprot(record, probe_path, contig_path)
                            if (
                                not cds.is_empty()
                                and cds.get_global_sim() > self.min_global_sim
                            ):
                                self.filtered_scaffolds.append(contig)
                                valid_components.add(component)
                                used_probes.append(probe_name)
                                break
                        else:
                            rejected_components.add(component)
                else:
                    for probe_name in self.probes_dict.keys():
                        probe_path = tmp_probe_paths[probe_name]
                        contig_path = Path(tmpdirname) / f"{contig}.fas"
                        cds = self._run_miniprot(record, probe_path, contig_path)
                        if (
                            not cds.is_empty()
                            and cds.get_global_sim() > self.min_global_sim
                        ):
                            self.filtered_scaffolds.append(contig)
                            valid_components.add(component)
                            used_probes.append(probe_name)
                            break
                    else:
                        rejected_components.add(component)

    def save(self, fasta_file: str) -> None:
        """Save the trimmed sequences in a fasta file"""
        print(f"Saving filtered scaffolds into {fasta_file}")
        records = [self.contigs_dict[rec] for rec in self.filtered_scaffolds]
        SeqIO.write(records, fasta_file, "fasta")


@dataclass
class OverlappingCds(MiniprotInit):
    """
    Data structure for selecting the best probe over the entire set of scaffolds.
    The probe is then used as reference for separating non overlapping scaffolds.
    Attributes
    ----------
    Post Init
    -cds_dict: dict of contig name as key and list of Cds object as value.
    -filtered_cds_dict: dict of contig name as key and the Cds that corresponds to the best overall probe,
        discarding sequences without matches
    -non_overlapping: list of OverlappingSeqs objects that keep track of overlapping sequences and position on probe.
    -best_probe: name of the probe version that has the overall best mapping score.
    """

    cds_dict: defaultdict[str, list[Cds]] = field(
        init=False, repr=False, default_factory=lambda: defaultdict(list)
    )
    filtered_cds_dict: dict[str, Cds] = field(
        init=False, repr=False, default_factory=dict
    )
    non_overlapping: list[OverlappingSeqs] = field(init=False, default_factory=list)
    best_probe: str = field(init=False, default=None)

    def __post_init__(self):
        super().__post_init__()
        # run miniprot on each scaffold/probe combination
        with tempfile.TemporaryDirectory() as tmpdirname:
            tmp_probe_paths = {}
            for probe_name, record in list(self.probes_dict.items())[:]:
                fasta_name = f"{probe_name}.fas"
                contig_fasta = Path(tmpdirname) / fasta_name
                SeqIO.write(record, contig_fasta, "fasta")
                tmp_probe_paths[probe_name] = contig_fasta

            for contig, record in self.contigs_dict.items():
                print(f"working on contig {contig}")
                for probe_name, probe_path in tmp_probe_paths.items():
                    print(f"working on probe {probe_name}")
                    fasta_name = f"{contig}.fas"
                    contig_path = Path(tmpdirname) / fasta_name
                    cds = self._run_miniprot(record, probe_path, contig_path)
                    cds.probe_name = probe_name
                    if not cds.is_empty():
                        self.cds_dict[contig].append(cds)
                # print("cds",self.cds_dict.get(contig))
            self._compute_rank_score()
            self._sorting_overlapping()
            print("non_overlapping", self.non_overlapping)

    def _compute_rank_score(self) -> Self:
        """
        Find the probe that has the best overall match using the miniprot score.
        Populate the 'best_probe' attribute.
        """
        if not self.cds_dict.values():
            return self
        probe_ranks = defaultdict(RankCoverage)  # score per contig
        for contig, list_cds in self.cds_dict.items():
            probe_scores = sorted(
                [(cds.probe_name, cds.get_global_score()) for cds in list_cds],
                key=lambda x: x[1],
                reverse=True,
            )
            for idx, probe_score in enumerate(probe_scores):
                probe_name = probe_score[0]
                probe_ranks[probe_name].append_to(idx, contig)

        best_combination = sorted(
            probe_ranks.items(), key=lambda x: x[1].get_total_score()
        )[0]
        self.best_probe = best_combination[0]
        print(f"Best overall probe: {self.best_probe}")
        for contig in best_combination[1].coverage:
            self.filtered_cds_dict[contig] = [
                cds
                for cds in self.cds_dict[contig]
                if cds.probe_name == self.best_probe
            ][0]
        return self

    def _sorting_overlapping(self) -> Self:
        """
        Separate the Cds by overlap with the probe sequence.
        Populate the non_overlapping attribute
        :return:
        """
        intervals = []
        for contig, cds in self.filtered_cds_dict.items():
            if cds.is_empty():
                continue
            start = cds.probe_start
            end = cds.probe_end
            intervals.append(OverlappingSeqs(start, end, [contig]))
        self.non_overlapping = merge_intervals(intervals)
        return self

    def save_scaffolds(self, fasta_folder: str) -> None:
        """
        Save the scaffolds from each OverlappingSeqs object into a separate fasta file.
        Scaffolds are trimmed in order to fit the regions covered by the probe.
        Files are labelled by the aa coordinates covered on the probe
        :param fasta_folder: folder where the fasta files are saved.
        :return:
        """
        assert (
            self.non_overlapping
        ), f"[ERROR] No overlap detected in {self.contigs_fasta}"

        print(f"Saving clusters to fasta folder {fasta_folder}")
        Path(fasta_folder).mkdir(parents=True, exist_ok=True)
        fasta_prefix = self.probes_path.stem
        idx = 0
        for merged in self.non_overlapping:
            probe_interval = merged.get_probe_interval()
            trimmed_records = self._trim_msa(merged.get_seq_names())
            if trimmed_records:
                name = f"{fasta_prefix}_{probe_interval}_{idx}.fasta"
                name_fasta = Path(fasta_folder) / name
                SeqIO.write(trimmed_records, name_fasta, "fasta")
                print(f"saving fasta {name_fasta}")
                idx += 1

    def save_best_probe(self, fasta_file: str) -> None:
        fasta_path = Path(fasta_file)
        fasta_path.parent.mkdir(parents=True, exist_ok=True)
        print(f"Saving best probe to fasta file {fasta_file}")
        SeqIO.write([self.probes_dict[self.best_probe]], fasta_path, "fasta")

    def _trim_msa(self, cluster_names: list[str]) -> Optional[list[SeqRecord]]:
        """
        Given a list of records, trim them to the smallest region of the probe shared by at least two sequences.
        In case, the group contains a single member, the sequence is trimmed to fit the probe's boundaries.
        :param cluster_names:
        :return: None is there is no match with a probe on this cluster, otherwise returns the trim sequence(s).
        """
        if len(cluster_names) == 1:
            contig = cluster_names[0]
            cds = self.filtered_cds_dict[contig]
            record = self.contigs_dict[contig]
            if cds.is_empty():
                return
            if cds.get_global_sim() < self.min_probe_contig_sim:
                return
            strand = cds.fragments[0].get_strand()
            if strand == 1:
                trim_start, trim_end = cds.mRNA_start, cds.mRNA_end
            else:
                trim_start, trim_end = cds.mRNA_end, cds.mRNA_start
            if abs(trim_end - trim_start) < self.min_contig_length:
                return

            if strand == -1:
                new_record = SeqRecord(
                    record.seq[trim_end:trim_start].reverse_complement(),
                    id=contig,
                    name="",
                    description="",
                )
            else:
                new_record = SeqRecord(
                    record.seq[trim_start:trim_end],
                    id=contig,
                    name="",
                    description="",
                )
            return [new_record]

        trimmed_records = []
        interval_tree = IntervalST()
        #  Find overlapping edges.
        for contig in cluster_names:
            cds = self.filtered_cds_dict[contig]
            if cds.is_empty():
                continue
            fragments = cds.get_fragments()
            start = list(fragments[0].get_correspondence().keys())[0]
            end = list(fragments[-1].get_correspondence().keys())[-1]
            interval_tree.put(Interval(start, end), contig)
        common_intervals = interval_tree.get_all_intersections()
        probe_start = min([overlap.get_start() for overlap in common_intervals])
        probe_end = max([overlap.get_end() for overlap in common_intervals])

        #  Trim the sequences
        for contig in cluster_names:
            cds = self.filtered_cds_dict[contig]
            if cds.is_empty():
                continue
            record = self.contigs_dict[contig]
            fragments = cds.get_fragments()
            strand = cds.fragments[0].get_strand()
            contig_start = list(fragments[0].get_correspondence().keys())[0]
            contig_end = list(fragments[-1].get_correspondence().keys())[-1]
            if strand == 1:
                trim_start, trim_end = cds.mRNA_start, cds.mRNA_end
            else:
                trim_start, trim_end = cds.mRNA_end, cds.mRNA_start
            # Case contig within core region
            if probe_start < contig_start < contig_end < probe_end:
                pass
            if contig_start <= probe_start:
                first, last = get_nuc_coordinates(probe_start, contig_end, fragments)
                if None in [first, last]:
                    continue
                trim_start = first if strand == 1 else last

            if probe_end <= contig_end:
                first, last = get_nuc_coordinates(
                    contig_start, probe_end - 1, fragments
                )
                if None in [first, last]:
                    continue
                trim_end = last if strand == 1 else first

            if abs(trim_end - trim_start) < self.min_contig_length:
                continue

            if strand == -1:
                new_record = SeqRecord(
                    record.seq[trim_end:trim_start].reverse_complement(),
                    id=contig,
                    name="",
                    description="",
                )
            else:
                new_record = SeqRecord(
                    record.seq[trim_start:trim_end],
                    id=contig,
                    name="",
                    description="",
                )
            trimmed_records.append(new_record)
        return trimmed_records


def snakemake_call(snakemake):
    with open(snakemake.log[0], "w") as outlog:
        sys.stdout = sys.stderr = outlog

        outdir = Path(snakemake.output[0])
        outdir.mkdir(exist_ok=True)
        probes_dir = Path(snakemake.input.probes_dir)

        for scfs in Path(snakemake.input.scfs_dir).glob("*.fasta"):
            probes = probes_dir / scfs.name
            olc = OverlappingCds(probes, scfs, **snakemake.params)

            if not olc.non_overlapping:
                print(f"Scaffolds from {probes.stem} did not yield any overlap")
                (outdir / "scfs" / scfs.stem).touch(exist_ok=True)
                continue

                ## Save scfs
            olc.save_scaffolds(outdir / "scfs" / scfs.stem)

            ## Save Probes
            olc.save_best_probe(outdir / f"probes/{probes.name}")


def main(): ...


if __name__ == "__main__":
    if "snakemake" in globals():
        snakemake_call(snakemake)
    else:
        main()
