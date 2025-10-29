#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2025
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.
"""
Functions and classes for running miniprot.

The class OverlappingCds is used after the vsearch clustering of SAUTE scaffolds and consensus estimation.
It maps all probe versions against the consensus sequences of vsearch clusters. It seeks
to find the best mapping probe version over all the consensus sequences and use it as
a common reference in order to determine which sequences are not overlapping.
Strictly non-overlapping sequences (there are no bridging scaffolds) are saved in separate files.
Sequences that do not meet miniprot length and similarity thresholds are discarded.
Exon-intron boundaries are then inferred with miniprot-boundary-scored (https://anaconda.org/bioconda/miniprot-boundary-scorer)
that rely on a AA substitution matrix.
The minimum similarity information for miniprot to keep a scaffold is passed via a dictionary with accession names
as keys and similarity thresholds as values.
#  Usage example:
    probe_fasta = ".../test_data/test_clustering/probe_3_aa.fasta"
    consensus_contig_fasta = "../test_data/test_clustering/gene_3_consensus.fasta"
    matrix = "../test_data/test_paralogy_2/blosum62.csv"
    thresholds_dict = {"accession_0":0.8, accession_1":0.9, "accession_2":0.75, "accession_3":0.65}}
    parameters = [8, 0.1, 10, 0.7] # See attribute definitions in the MiniprotInit class.
    # load the data, perform the computation in order to separate non-overlapping sequences:
    olc = OverlappingCds(probe_fasta, consensus_contig_fasta, matrix, thresholds_dict, *parameters)
    # Non overlapping sets of sequences are saved in a separate file.
    out_dir = "../test_data/test_clustering/nonoverlapping"
    min_exon_length = 10
    olc.save_records(out_dir, 10)
    # save the sequence of the best overall mapping probe.
    olc.save_best_probe("../test_data/test_clustering/gene_3_best_probe.fasta")

TODO: refactor find_global_boundaries() so that the alignments limits are inferred from the column's occupancy
instead of the exon's boundaries.
"""

__version__ = "0.1"

# =======================================================================================
#               IMPORTS
# =======================================================================================

from dataclasses import dataclass, field
from collections import defaultdict, namedtuple
from io import StringIO
import json
from operator import attrgetter, add
from pathlib import Path
import subprocess
import sys
import tempfile
from typing import Optional, Self, Literal

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
from intervaltree import Interval, IntervalTree

from gff_parser import Cds
from exon_intron import Exon

Region = namedtuple("Region", ["start", "end"])

MAX_EXPANSION_INTERVAL = (
    10  # on scaffold size of flanking region on each side of the probe hits,
)
# used to find the common region that matches a given probe across multiple scaffolds.


# =======================================================================================
#               FUNCTIONS
# =======================================================================================


def cluster_intervals(intervals: list[Interval]) -> list[Interval]:
    """
    Cluster overlapping Intervals using an interval tree.
    The data field of intervals is made out of a list of EndPoint objects.
    """
    it = IntervalTree.from_tuples(intervals)
    it.merge_overlaps(strict=True, data_reducer=add)
    return it


def select_best_limit(
    intervals: list[Interval], direction: Literal["left", "right"]
) -> int:
    """
    Select the interval that has the most number of sequences supporting it.
    If several scaffolds share this limit, return either the smallest if 'left' is set otherwise the largest.
    """
    sorted_intervals = (
        sorted(intervals, key=lambda i: i.begin)
        if direction == "left"
        else sorted(intervals, key=lambda i: i.begin, reverse=True)
    )
    best_limit = sorted_intervals[0]
    covered_sequences = len(best_limit.data)
    for interval in sorted_intervals[1:]:
        if covered_sequences > 3:
            break
        covered_sequences += len(interval.data)
        best_limit = interval
    else:
        best_limit = sorted_intervals[0]
    if direction == "left":
        return sorted(best_limit.data, key=lambda i: i.idx)[0].idx
    return sorted(best_limit.data, key=lambda i: i.idx)[-1].idx


def merge_intervals(arr: list["OverlappingSeqs"]) -> list["OverlappingSeqs"]:
    """
    Combine overlapping intervals of sequences and cluster their names.
    We can thus separate the group of sequences according to non-overlapping intervals.
    Used to cluster overlapping Exons from different scaffolds.
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
            index += 1
            arr[index] = arr[i]
    return [arr[index] for index in range(len(arr)) if index not in index_to_del]


def remove_gff(miniprot_out) -> bytes:
    """
    Transform a PAF file generated by 'miniprot --gff --aln' that is used to create Cds objects,
    into a file that would be obtained with the option 'miniprot --aln' that can be used as
    input for 'miniprot --gff --aln'.
    """
    new_lines = []
    for line in miniprot_out:
        if line.startswith("##PAF"):
            new_lines.append(line.replace("##PAF	", ""))
        if (
            line.startswith("##ATN")
            or line.startswith("##ATA")
            or line.startswith("##AAS")
            or line.startswith("##AQA")
        ):
            new_lines.append(line)
    return "".join(new_lines).encode()


def run_miniprot(
    probe_path: Path,
    scaffold_path: Path,
    threads=2,
    min_coverage=0.01,
) -> StringIO:
    """
    Wrapper for running miniprot.
    :return:
    """
    miniprot_cmd = f"miniprot -t {threads} --gff --aln --outn 1 -J 50 -E 3 --outc {min_coverage} {scaffold_path} {probe_path}".split()
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


def run_miniprot_boundary_scorer(
    miniprot: bytes, boundary_scorer_out: Path, matrix_path: Path
):
    """
    Wrapper for running miniprot_boundary_scorer.
    :return:
    """
    boundary_scorer_cmd = (
        f"miniprot_boundary_scorer -o '{boundary_scorer_out}' -s {matrix_path}"
    )
    try:
        subprocess.run(
            boundary_scorer_cmd,
            timeout=1000,
            shell=True,
            input=miniprot,
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


# =======================================================================================
#               CLASSES
# =======================================================================================


@dataclass
class RankCoverage:
    """
    Container for keeping track of the rank of a given probe out of all probes ranks.
    Probes are ranked by similarity to the scaffolds. The most similar is ranked first, the second most similar
    second, etc. The sum of ranks over all scaffolds gives the final rank.
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

    def get_mean_score(self, nbr_scaffolds: int):
        """The  score used for the global ranking is based on probe positions, the total number of
        scaffolds and the scaffolds without a match.
        -nbr_scaffolds: the total number of scaffolds that are subjected to miniprot.
        Missing matches get a penalty equal to the maximum rank,
         e.g. rank_score=[0, 3, 1, 0], coverage=['contig_0', 'contig_1', 'contig_2', 'contig_3']
         and nbr_contigs=6 (i.e. missing 'contig_4', 'contig_6')  produces a score of (3+1+6+6)/6
        """
        nbr_missing_scaffolds = nbr_scaffolds - len(self.coverage)
        missing_scores = [nbr_scaffolds] * nbr_missing_scaffolds
        return (sum(self.rank_score) + sum(missing_scores)) / nbr_scaffolds


@dataclass(frozen=True, slots=True)
class ExonCorrespondence:
    """
    Probe to scaffold coordinates
    Attributes
    ----------
    exon_probe: start-end coordinates on the probe.
    exon_scaffold: start-end coordinates on the scaffold.
    length_on_scaffold: nbr of AAs that are aligned between the probe and the scaffold.
    """

    exon_probe: Exon
    exon_scaffold: Exon
    length_on_scaffold: int


@dataclass(slots=True)
class Boundary:
    """
    coordinate on scaffold
    """

    scaffold_name: str
    scaffold_start: int
    scaffold_end: int


@dataclass(slots=True)
class EndPoint:
    """
    Container for start or end of an exon.
    """

    idx: int
    name: str


@dataclass
class ExtendedCds(Cds):
    """
    Cds object with additional features to handle paralogy information
    Post Init
    -correspondence: scaffold coordinates to probe coordinates
    -rev_correspondence: probe coordinates to scaffold coordinates
    -exon_correspondences: list of ExonCorrespondence(exon_probe=Exon(start=71, end=110),
        exon_scaffold=Exon(start=17, end=132), length_on_scaffold=115) that contains the
        coordinates of exons on probe and scaffolds.
    #-accepted_exons: Most common exons as identified on the probes.
    """

    correspondence: dict[int, int] = field(default_factory=dict, init=False, repr=True)
    rev_correspondence: dict[int, int] = field(
        default_factory=dict, init=False, repr=True
    )
    exon_correspondences: list[ExonCorrespondence] = field(
        default_factory=list, init=False, repr=True
    )

    def __post_init__(self):
        super().__post_init__()
        self._reverse_correspondence()

    def _reverse_correspondence(self) -> Self:
        """
        Populate the correspondence, rev_correspondence, start_on_probe and end_on_probe attributes.
        """
        for fragment in self.fragments:
            for k, v in fragment.correspondence.items():
                self.correspondence[k] = v
                self.rev_correspondence[v] = k
        self.start_on_probe = min(self.rev_correspondence)
        self.end_on_probe = max(self.rev_correspondence)
        return self

    def _find_probe_exons(self, scaff_exons: list[Exon]) -> Self:
        """
        Find coordinates pairs on the probe where exons are located from scaffold coordinates.
        Because there might be boundaries slight variation due to  miniprot alignment,
        we search a few AA around the probe edge to find a correspondence on the scaffold.
        """
        for scaff_exon in scaff_exons:
            probe_start = None
            probe_end = None
            for idx in [0, 1, -1, 2, -2]:
                if (
                    probe_start := self.rev_correspondence.get(scaff_exon.start + idx)
                ) is not None:
                    break
            for idx in [0, 1, -1, 2, -2]:
                if (
                    probe_end := self.rev_correspondence.get(scaff_exon.end + idx)
                ) is not None:
                    break
            if None not in [probe_start, probe_end]:
                probe_exon = Exon(probe_start, probe_end)
                length = scaff_exon.end - scaff_exon.start
                self.exon_correspondences.append(
                    ExonCorrespondence(probe_exon, scaff_exon, length)
                )
        return self


@dataclass
class OverlappingSeqs:
    """
    Container for groups of sequences that overlap a common region of the probe.
    Attributes
    ----------
    -probe_start: the smallest index on the probe (AA space) that has a matching sequence.
    -probe_end: the largest index on the probe (AA space) that has a matching sequence.
    -seq_names: a list of sequence names that overlap this region.
    -cds_dict
    -common_exons: exons found in a majority of scaffolds
    -paralogs: list of scaffolds that are deemed to be paralogs.
    To do: in the final paralog selection add the sequences that are present in OverlappingCds.cds_dict
    but disappear in filtered_cds_dict.
    -scaffold_exons: for each exon region keep a list of scaffold Boundaries that correspond to this exon.
    -scaffold_introns: for each intron region keep a list of scaffold Boundaries that correspond to this intron.
    """

    probe_start: int
    probe_end: int
    seq_names: list[str]
    extended_cds_dict: dict[str, ExtendedCds] = field(
        init=False, repr=False, default_factory=dict
    )
    global_start: Optional[int] = field(init=False)
    global_end: Optional[int] = field(init=False)

    def find_global_boundaries(self):
        """
        For each non overlapping object, find the most common endpoint at the 5' and 3'
        exon endpoints.
        """
        start_intervals = []
        end_intervals = []
        if not self.extended_cds_dict:
            self.global_start, self.global_end = None, None
            return
        for scaffold_name, cds in self.extended_cds_dict.items():
            start_intervals.append(
                Interval(
                    cds.probe_start - MAX_EXPANSION_INTERVAL,
                    cds.probe_start + MAX_EXPANSION_INTERVAL,
                    [EndPoint(cds.probe_start, scaffold_name)],
                )
            )
            end_intervals.append(
                Interval(
                    cds.probe_end - MAX_EXPANSION_INTERVAL,
                    cds.probe_end + MAX_EXPANSION_INTERVAL,
                    [EndPoint(cds.probe_end, scaffold_name)],
                )
            )
        self.global_start = select_best_limit(
            cluster_intervals(start_intervals), "left"
        )
        self.global_end = select_best_limit(cluster_intervals(end_intervals), "right")
        print(f"{self.global_start=}\t{self.global_end =}")

    def combine(self, other: list[str]) -> Self:
        """
        Add items to seq_names
        """
        self.seq_names.extend(other)
        return self

    def __repr__(self):
        return (
            f"probe_start={self.probe_start}, probe_end={self.probe_end},"
            f"{sorted(self.seq_names)}"
        )


@dataclass
class MiniprotInit:
    """
    Base class for miniprot initialization.
    Attributes
    ----------
    -probes_fasta: path to the amino acid probes fasta file (multiprobe)
    -scaffold_fasta: path to the nucleotide scaffold fasta file.
    -threads: for miniprot.
    -min_fragment_cov: fraction of the probe covered by a contig for miniprot.
    -min_exonic_length: minimum length of concatenated exons after trimming used for saving.

    Post Init
    -probes_path: probe_fasta string converted to Path.
    -scaffold_path: scaffold_fasta string converted to Path.
    -probes_dict: dict of probe name as key and SeqRecord as value.
    -scaffold_dict: dict of scaffold name as key and SeqRecord as value.

    """

    probes_fasta: str
    scaffold_fasta: str
    substitution_matrix: str
    threads: int = field(default=8)
    min_fragment_cov: float = field(default=0.05)
    min_exonic_length: int = field(default=200)
    probes_path: Path = field(init=False, repr=False)
    scaffold_path: Path = field(init=False, repr=False)
    probes_dict: dict[str, SeqRecord] = field(
        init=False, repr=False, default_factory=dict
    )
    scaffold_dict: dict[str, SeqRecord] = field(
        init=False, repr=False, default_factory=dict
    )

    def __post_init__(self):
        self.probes_path = Path(self.probes_fasta)
        self.scaffold_path = Path(self.scaffold_fasta)
        self.substitution_matrix_path = Path(self.substitution_matrix)
        assert self.probes_path.exists(), f"{self.probes_fasta} does not exist"
        assert self.scaffold_path.exists(), f"{self.scaffold_fasta} does not exist"
        assert (
            self.substitution_matrix_path.exists()
        ), f"{self.substitution_matrix} does not exist"
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
            self.scaffold_dict = SeqIO.to_dict(SeqIO.parse(self.scaffold_path, "fasta"))
        except Exception as err:
            sys.exit(f"[ERROR] {err}")

    def _run_miniprot(
        self, record, probe_path_temp, scaffold_path_temp, exon_scorer=False
    ) -> StringIO:
        SeqIO.write(record, scaffold_path_temp, "fasta")
        miniprot_out = run_miniprot(
            probe_path_temp,
            scaffold_path_temp,
            self.threads,
            self.min_fragment_cov,
        )
        #print(list(miniprot_out))
        return miniprot_out


@dataclass
class OverlappingCds(MiniprotInit):
    """
    Data structure for selecting the best probe over the entire set of scaffolds.
    The probe is then used as reference for separating non overlapping scaffolds.
    Attributes
    ----------
    -user_probe: if specified, the probe with the best matching score to the scaffold is not computed with miniprot.
    -min_global_identity_dict: per sample dictionary of identity thresholds to filter scaffolds, using
        the global identity found over several fragments.
    -min_global_identity: identity value to use if sample is not present in min_global_identity_dict
    -max_intron_dict: dictionary with per sample values for maximum intron length allowed to keep a scaffold.
    -max_intron_length: max intron length to use if a sample is not present in max_intron_dict.
    #-min_overlapping: proportion of the exon length that is used to find an overlap with another exon.
        The length that is explored on each endpoint of the exon is capped with hard coded boundary.
    Post Init
    -cds_dict: dict of contig name as key and list of Cds object as value.
    -filtered_cds_dict: dict of scaffolds names as key and the Cds that corresponds to the best overall probe,
        discarding sequences without matches
    -non_overlapping: list of OverlappingSeqs objects that keep track of overlapping sequences and position on probe.
    -best_probe: name of the probe version that has the overall best mapping score.
    """

    user_probe: str = field(default=None)
    min_global_identity_dict: dict[str:float] = field(default_factory=dict)
    min_global_identity: float = field(default=0.85)
    max_intron_length: int = field(default=1000)
    max_intron_dict: dict[str, int] = field(default_factory=dict)

    miniprot_out: Optional[defaultdict[str, dict[str, StringIO]]] = field(
        init=False, repr=False, default_factory=lambda: defaultdict(dict)
    )
    cds_dict: Optional[defaultdict[str, list[Cds]]] = field(
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
        print(f"{self.user_probe=}")
        with tempfile.TemporaryDirectory() as tmpdirname:
            tmp_probe_paths = {}
            # write all probe records in separate file
            if (
                self.user_probe is not None
                and self.user_probe in self.probes_dict.keys()
            ):
                print(f"User selected probe is {self.user_probe}")
                fasta_name = f"{self.user_probe}.fas"
                probe_fasta = Path(tmpdirname) / fasta_name
                SeqIO.write(self.probes_dict[self.user_probe], probe_fasta, "fasta")
                tmp_probe_paths[self.user_probe] = probe_fasta

            else:
                print("creating probe list")
                for probe_name, record in list(self.probes_dict.items())[:]:
                    fasta_name = f"{probe_name}.fas"
                    probe_fasta = Path(tmpdirname) / fasta_name
                    SeqIO.write(record, probe_fasta, "fasta")
                    tmp_probe_paths[probe_name] = probe_fasta
            for scaffold, record in self.scaffold_dict.items():
                accession = scaffold.split("|")[0]

                print(f"working on scaffold {scaffold}")
                for probe_name, probe_path in tmp_probe_paths.items():
                    print(f"working on probe {probe_name}")
                    fasta_name = f"{scaffold}.fas"
                    scaffold_path = Path(tmpdirname) / fasta_name
                    miniprot_result = self._run_miniprot(
                        record, probe_path, scaffold_path
                    )
                    if scaffold in self.miniprot_out:
                        self.miniprot_out[scaffold].update(
                            {probe_name: miniprot_result}
                        )
                    else:
                        self.miniprot_out[scaffold] = {probe_name: miniprot_result}
                    cds = Cds(miniprot_result)
                    cds.probe_name = probe_name
                    # print(f"{accession=} with identity {self.min_global_identity_dict.get(accession, 0.85)}")
                    if (
                        not cds.is_empty()
                        and cds.global_identity
                        >= self.min_global_identity_dict.get(
                            accession, self.min_global_identity
                        )
                        and cds.get_longest_intron()
                        < self.max_intron_dict.get(accession, self.max_intron_length)
                    ):
                        self.cds_dict[scaffold].append(cds)
                print("cds", self.cds_dict.get(scaffold))
            scaffold_with_cds = sorted(
                [
                    scf
                    for scf in self.scaffold_dict
                    if self.cds_dict.get(scf) is not None
                ],
            )
            scaffold_without_cds = sorted(
                [scf for scf in self.scaffold_dict if scf not in scaffold_with_cds]
            )
            print(f"scaffold with cds:\n{"\n".join(scaffold_with_cds)}")
            print(f"scaffold  without cds:\n{"\n".join(scaffold_without_cds)}")
            print("computing rank score")
            self._compute_rank_score()
            print("initial sequences:", len(self.cds_dict))
            print("filtered sequences:", len(self.filtered_cds_dict))
            print("sorting overlapping")
            self._sorting_overlapping()
            print("done sorting overlapping")
            print("exon search")
            self._exon_search()
            print("exon extraction")
            # free memory
            self.miniprot_out = None
            self.cds_dict = None

    def _compute_rank_score(self) -> Self:
        """
        Find the probe that has the best overall match using the miniprot score.
        Populate the 'best_probe' attribute.
        """
        if not self.cds_dict.values():
            return self
        probe_ranks = defaultdict(RankCoverage)  # score per contig
        for scaffold, list_cds in self.cds_dict.items():
            probe_scores = sorted(
                [(cds.probe_name, cds.get_global_score()) for cds in list_cds],
                key=lambda x: x[1],
                reverse=True,
            )
            for idx, probe_score in enumerate(probe_scores):
                probe_name = probe_score[0]
                probe_ranks[probe_name].append_to(idx, scaffold)

        best_combination = sorted(
            probe_ranks.items(), key=lambda x: x[1].get_mean_score(len(self.cds_dict))
        )[0]

        self.best_probe = best_combination[0]
        print(f"Best overall probe: {self.best_probe}")
        # Keep only scaffolds that are covered by the best probe
        for scaffold in best_combination[1].coverage:
            self.filtered_cds_dict[scaffold] = [
                cds
                for cds in self.cds_dict[scaffold]
                if cds.probe_name == self.best_probe
            ][0]
        return self

    def _sorting_overlapping(self) -> Self:
        """
        Separate the Cds by overlap with the probe sequence.
        Populate the non_overlapping attribute
        """
        interval_list = []
        for scaffold, cds in self.filtered_cds_dict.items():
            if cds.is_empty():
                continue
            start = cds.probe_start
            end = cds.probe_end
            interval_list.append(OverlappingSeqs(start, end, [scaffold]))
        self.non_overlapping = merge_intervals(interval_list)
        return self

    def _boundary_scorer_parser(self, boundary_scorer_path: Path) -> list[Exon]:
        """
        Parser of the miniprot_boundary_scorer output.
        Create for each record a list of Exon objects.
        """
        lines = boundary_scorer_path.read_text().split("\n")
        exons = []
        for line in lines:
            if not line:
                continue
            splt = line.split("\t")
            if splt[2] == "CDS" and splt[6] == "+":
                exons.append(Exon(int(splt[3]) - 1, int(splt[4]) - 3))
        return exons

    def _exon_search(self):
        """run miniprot boundary scorer"""
        if not self.non_overlapping:
            print("No overlapping group")
            return
        with tempfile.TemporaryDirectory() as tmpdirname:
            for overlapping_gp in self.non_overlapping:
                for scaffold_name in overlapping_gp.seq_names:
                    miniprot_result = self.miniprot_out[scaffold_name][self.best_probe]
                    miniprot_result.seek(0)
                    cds = ExtendedCds(miniprot_result)
                    boundary_scorer_out = (
                        Path(tmpdirname) / f"{scaffold_name}_boundary.gfa"
                    )
                    miniprot_result.seek(0)
                    # run miniprot_boundary scorer
                    try:
                        run_miniprot_boundary_scorer(
                            remove_gff(miniprot_result),
                            boundary_scorer_out,
                            self.substitution_matrix_path,
                        )
                    except:
                        continue
                    cds._find_probe_exons(
                        self._boundary_scorer_parser(boundary_scorer_out)
                    )
                    if not cds.exon_correspondences:
                        print(
                            f"Something went terribly wrong with {scaffold_name} as no Exon was inferred"
                        )
                        continue
                    overlapping_gp.extended_cds_dict[scaffold_name] = cds
                overlapping_gp.find_global_boundaries()

    def save_records(self, outdir: Path, min_exon_size: int) -> None:
        """
        Save exonic regions, supercontigs (the whole scaffold), and genetigs that include the sequence
        within the probe boundaries.
        """
        stem_name = self.probes_path.stem
        print(f"saving {stem_name}")
        for overlapping_gp in self.non_overlapping:
            prefix = (
                f"{stem_name}_{overlapping_gp.probe_start}_{overlapping_gp.probe_end:}"
            )

            if overlapping_gp.extended_cds_dict:
                if not overlapping_gp.global_start and not overlapping_gp.global_end:
                    sys.exit("[ERROR] 'global_start' and 'global_end' are undefined")

                def _build_description()->dict:
                    """
                    Building the description part of the header, inspired by the way captus reports the information.
                    """
                    name_description_dict = {}
                    for scaffold_name, extended_cds in overlapping_gp.extended_cds_dict.items():
                        probe_cov = sum([frag.query_end -frag.query_start for frag in extended_cds.fragments])
                        probe_length = len(self.probes_dict[self.best_probe].seq)
                        identity=extended_cds.global_identity
                        score=extended_cds.get_global_score()
                        description = f"[query={self.best_probe}] [cover={probe_cov/probe_length:.2f}] [ident={identity}] [score={score}]"
                        print(f"{scaffold_name},  {probe_cov/probe_length} {identity} {score}")
                        name_description_dict[scaffold_name] = description
                    return name_description_dict

                #  Save exons
                name = f"{prefix}_exons.fasta"
                out_path_exon = outdir / name
                global_start = overlapping_gp.global_start
                global_end = overlapping_gp.global_end
                exon_records: list[SeqRecord] = []
                description_dict = _build_description()

                for (
                    scaffold_name,
                    extended_cds,
                ) in overlapping_gp.extended_cds_dict.items():
                    seq = str(self.scaffold_dict[scaffold_name].seq)
                    new_seq = ""
                    for exon_correspondence in extended_cds.exon_correspondences:
                        probe_start = exon_correspondence.exon_probe.start
                        probe_end = exon_correspondence.exon_probe.end
                        exon_start = exon_correspondence.exon_scaffold.start
                        exon_end = exon_correspondence.exon_scaffold.end
                        # Adjusting the exon boundaries in order to match the most common exon
                        if global_start > probe_start:
                            for idx in [0, 1, -1, 2, -2]:
                                if (
                                    tmp_start := extended_cds.correspondence.get(
                                        global_start + idx
                                    )
                                ) is not None:
                                    exon_start = tmp_start
                                    break
                        if probe_end > global_end:
                            for idx in [0, -1, 1, -2, 2]:
                                if (
                                    tmp_end := extended_cds.correspondence.get(
                                        global_end + idx
                                    )
                                ) is not None:
                                    exon_end = tmp_end
                                    break
                        new_seq += seq[exon_start:exon_end]
                    if len(new_seq) > min_exon_size:
                        exon_records.append(
                            SeqRecord(
                                Seq(new_seq),
                                name="",
                                description=description_dict[scaffold_name] + f" [length={len(new_seq)}]",
                                id=scaffold_name,
                            )
                        )
                SeqIO.write(
                    sorted(exon_records, key=lambda x: x.id),
                    out_path_exon,
                    "fasta",
                )
                #  Save supercontigs
                super_name = f"{prefix}_supercontigs.fasta"
                out_path_supercontigs = outdir / super_name
                super_records = []
                for scaffold_name in overlapping_gp.extended_cds_dict:
                    record = self.scaffold_dict[scaffold_name]
                    record.description = description_dict[scaffold_name] + f" [length={len(record.seq)}]"
                    super_records.append(record)

                if super_records:
                    SeqIO.write(
                        sorted(super_records, key=lambda x: x.id),
                        out_path_supercontigs,
                        "fasta",
                    )
                # save genetigs
                genetigs_name = f"{prefix}_genetigs.fasta"
                out_path_genetigs = outdir / genetigs_name
                genetigs_records = []

                for (
                    scaffold_name,
                    extended_cds,
                ) in overlapping_gp.extended_cds_dict.items():
                    seq = str(self.scaffold_dict[scaffold_name].seq)
                    new_scaffold = ""
                    first_exon = extended_cds.exon_correspondences[0]
                    probe_start = first_exon.exon_probe.start
                    scaffold_start = first_exon.exon_scaffold.start
                    last_exon = extended_cds.exon_correspondences[-1]
                    probe_end = last_exon.exon_probe.end
                    scaffold_end = last_exon.exon_scaffold.end
                    if global_start > probe_start:
                        for idx in [0, 1, -1, 2, -2]:
                            if (
                                tmp_start := extended_cds.correspondence.get(
                                    global_start + idx
                                )
                            ) is not None:
                                scaffold_start = tmp_start
                                break
                    if probe_end > global_end:
                        for idx in [0, -1, 1, -2, 2]:
                            if (
                                tmp_end := extended_cds.correspondence.get(
                                    global_end + idx
                                )
                            ) is not None:
                                scaffold_end = tmp_end
                                break
                    new_scaffold += seq[scaffold_start:scaffold_end]
                    if len(new_scaffold) > min_exon_size:
                        genetigs_records.append(
                            SeqRecord(
                                Seq(new_scaffold),
                                name="",
                                description=description_dict[scaffold_name] + f" [length={len(new_scaffold)}]",
                                id=scaffold_name,
                            )
                        )

                if genetigs_records:
                    SeqIO.write(
                        sorted(genetigs_records, key=lambda x: x.id),
                        out_path_genetigs,
                        "fasta",
                    )

    def save_best_probe(self, fasta_path: Path) -> None:
        fasta_path.parent.mkdir(parents=True, exist_ok=True)
        if not self.best_probe:
            return
        print(f"Saving best probe to fasta file {fasta_path}")
        SeqIO.write([self.probes_dict[self.best_probe]], fasta_path, "fasta")


def snakemake_call(snakemake):
    with open(snakemake.log[0], "w") as outlog:
        sys.stdout = sys.stderr = outlog

        threads = snakemake.threads
        outdir = Path(snakemake.output[0])
        scfs = Path(snakemake.input.scfs)
        probes = Path(snakemake.input.probes)
        div_map = json.loads(Path(snakemake.input.div_map).read_bytes())
        max_intron_map = json.loads(Path(snakemake.input.max_intron_map).read_bytes())

        outdir.mkdir(parents=True, exist_ok=True)
        if scfs.stat().st_size > 0:
            olc = OverlappingCds(
                probes,
                scfs,
                min_global_identity_dict=div_map,
                max_intron_dict=max_intron_map,
                threads=threads,
                **snakemake.params,
            )

            ## Save scfs
            outdir.mkdir(parents=True, exist_ok=True)
            olc.save_records(outdir, snakemake.params.min_exonic_length)

            ## Save Probes
            olc.save_best_probe(outdir.parent.parent / f"probes/{probes.name}")


def main():

    if len(sys.argv) != 8:
        print(
            "Usage: python miniprot.py <probes.fasta> <scfs.fasta> <matrix.csv> <div_map.json> <intron_map.json> <outdir> <miniprot.log>"
        )
        sys.exit(1)

    class Run:
        def __init__(self, **kwargs):
            setattr(self, "_dict", kwargs)
            for key, val in kwargs.items():
                setattr(self, key, val)

        def keys(self):
            return self._dict.keys()

        def __getitem__(self, key):
            return self._dict[key]

    # Mock the snakemake object
    snakemake = Run(
        input=Run(
            probes=sys.argv[1],
            scfs=sys.argv[2],
            div_map=sys.argv[4],
            max_intron_map=sys.argv[5],
        ),
        output=[sys.argv[6]],
        log=[sys.argv[7]],
        threads=1,
        params=Run(
            substitution_matrix=sys.argv[3],
            min_fragment_cov=0.1,
            min_exonic_length=10,
        ),
    )

    snakemake_call(snakemake)


def debug():
    from random import uniform

    def mk_threshold_dict(fasta: str):
        records = SeqIO.parse(fasta, "fasta")
        return {rec.id.split("|")[0]: uniform(0.6, 1) for rec in records}

    matrix = (
        "/home/yjkbertrand/Documents/projects/Nextpyper/nextpyper/data/blosum62.csv"
    )
    probe_file = "/home/yjkbertrand/Documents/projects/nextpiper/debug/miniprot/refactor_7_30_25/tmp_multi/6544_probe.fasta"
    scfs = "/home/yjkbertrand/Documents/projects/nextpiper/debug/miniprot/refactor_7_30_25/tmp_multi/6544_scfs.fasta"
    parameters = [

        8,
        0.1,
        10,
        None,
        mk_threshold_dict(scfs),
        0.5
    ]
    olc = OverlappingCds(probe_file, scfs, matrix, *parameters)
    print("overlapping:", olc)
    outdir = Path(
        "/home/yjkbertrand/Documents/projects/nextpiper/debug/miniprot/key_error/outdir"
    )
    olc.save_records(outdir, 10)


if __name__ == "__main__":
    debug()
    # if "snakemake" in globals():
    #     snakemake_call(snakemake)
    # else:
    #     main()
