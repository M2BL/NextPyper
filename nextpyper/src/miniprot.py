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
Strictly non-overlapping sequences (there are no bridging scaffolds) are saved in separate
files.
Sequences that do not meet miniprot length and similarity thresholds are discarded.
Exon-intron boundaries are then inferred with miniprot-boundary-scored (https://anaconda.org/bioconda/miniprot-boundary-scorer)
that rely on a AA substitution matrix.
#  Usage example:
    probe_fasta = ".../test_data/test_clustering/probe_3_aa.fasta"
    consensus_contig_fasta = "../test_data/test_clustering/gene_3_consensus.fasta"
    matrix = "../test_data/test_paralogy_2/blosum62.csv"
    parameters = [8, 0.85, 0.1, 10, 0.7] # See attribute definitions in the MiniprotInit class.
    # load the data, perform the computation in order to separate non-overlapping sequences:
    olc = OverlappingCds(probe_fasta, consensus_contig_fasta, matrix, *parameters)
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
from itertools import chain
from operator import attrgetter
from pathlib import Path
import subprocess
import sys
import tempfile
import textwrap
from typing import Optional, Self, Literal

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
from intervaltree import Interval, IntervalTree

from gff_parser import Cds
from exon_intron import Exon

Region = namedtuple("Region", ["start", "end"])

MAX_EXPANSION_INTERVAL = 10


# =======================================================================================
#               FUNCTIONS
# =======================================================================================
def fuse_intervals(intervals: list[Interval]) -> Interval:
    """Fuse several Interval objects into a single Interval.
    Data attributes are fused into a single list, lower and higher bounds
    are set to min and max values respectively."""
    min_value = sorted(intervals, key=lambda i: i.begin)[0].begin
    max_value = sorted(intervals, key=lambda i: i.end)[-1].end
    data = list(chain(*[x.data for x in intervals]))
    return Interval(min_value, max_value, data)


def cluster_intervals(intervals: list[Interval]) -> list[Interval]:
    """
    Cluster overlapping Intervals using an interval tree.
    The data field of intervals is made out of a list of EndPoint objects.
    """
    it = IntervalTree.from_tuples(intervals)
    used_idxs = []
    for interval in intervals:
        idx = interval.data[0].idx
        if idx in used_idxs:
            continue
        centered_intervals = sorted(it[idx])
        new_interval = fuse_intervals(centered_intervals)
        del it[idx]
        it.add(new_interval)
        used_idxs.extend(x.idx for x in new_interval.data)
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
            index = index + 1
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
    treads=2,
    min_similarity=0.80,
    min_coverage=0.01,
) -> StringIO:
    """
    Wrapper for running miniprot.
    :return:
    """
    miniprot_cmd = f"miniprot -t {treads} --gff --aln --outn 1 -J 50 --outs {min_similarity} --outc {min_coverage} {scaffold_path} {probe_path}".split()
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
    print(boundary_scorer_cmd)
    try:
        subprocess.run(
            boundary_scorer_cmd,
            timeout=100,
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


@dataclass(frozen=True)
class ExonCorrespondence:
    """
    Probe to scaffold coordinates
    Attributes
    ----------
    exon_probe: start-end coordinates on the probe.
    exon_scaffold: start-end coordinates on the scaffold.
    length_on_scaffold: nbr of AAs that are aligned between the probe and the scaffold.
    """

    __slots__ = ["exon_probe", "exon_scaffold", "length_on_scaffold"]
    exon_probe: Exon
    exon_scaffold: Exon
    length_on_scaffold: int


@dataclass
class Boundary:
    """
    coordinate on scaffold
    """

    __slots__ = ["scaffold_name", "scaffold_start", "scaffold_end"]
    scaffold_name: str
    scaffold_start: int
    scaffold_end: int


@dataclass
class EndPoint:
    """
    Container for start or end of an exon.
    """

    __slots__ = ["idx", "name"]
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
    -treads: for miniprot.
    -min_probe_scaffold_sim: mapping similarity for miniprot.
    -min_fragment_cov: fraction of the probe covered by a contig for miniprot.
    -min_exonic_length: minimum length of concatenated exons after trimming used for saving.
    -min_global_identity: identity over all several fragments. Should be slightly lower than min_probe_contig_sim.

    Post Init
    -probes_path: probe_fasta string converted to Path.
    -scaffold_path: cscaffold_fasta string converted to Path.
    -probes_dict: dict of probe name as key and SeqRecord as value.
    -scaffold_dict: dict ofscaffold name as key and SeqRecord as value.

    """

    probes_fasta: str
    scaffold_fasta: str
    substitution_matrix: str
    treads: int = field(default=8)
    min_probe_scaffold_sim: float = field(default=0.80)
    min_fragment_cov: float = field(default=0.05)
    min_exonic_length: int = field(default=200)
    min_global_identity: float = field(default=0.00)
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
            self.treads,
            self.min_probe_scaffold_sim,
            self.min_fragment_cov,
        )
        return miniprot_out


@dataclass
class OverlappingCds(MiniprotInit):
    """
    Data structure for selecting the best probe over the entire set of scaffolds.
    The probe is then used as reference for separating non overlapping scaffolds.
    Attributes
    ----------
    -user_probe: if specified, the probe with the best matching score to the scaffold is not computed with miniprot.
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
    # min_overlapping = 0.1
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

                    if (
                        not cds.is_empty()
                        and cds.global_identity >= self.min_global_identity
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
                    run_miniprot_boundary_scorer(
                        remove_gff(miniprot_result),
                        boundary_scorer_out,
                        self.substitution_matrix_path,
                    )
                    cds._find_probe_exons(
                        self._boundary_scorer_parser(boundary_scorer_out)
                    )
                    if not cds.exon_correspondences:
                        print(
                            f"Something went terribly wrong with {scaffold_name} as no Exon was inferred"
                        )
                        continue
                    # assert (
                    #     cds.exon_correspondences
                    # ), f"Something went terribly wrong with {scaffold_name} as no Exon was inferred"

                    overlapping_gp.extended_cds_dict[scaffold_name] = cds
                overlapping_gp.find_global_boundaries()

    def save_records(self, outdir: Path, min_exon_size: int) -> None:
        stem_name = self.probes_path.stem
        print(f"saving {stem_name}")
        for overlapping_gp in self.non_overlapping:
            prefix = (
                f"{stem_name}_{overlapping_gp.probe_start}_{overlapping_gp.probe_end:}"
            )

            if overlapping_gp.extended_cds_dict:
                if not overlapping_gp.global_start and not overlapping_gp.global_end:
                    sys.exit(f"[ERROR] 'global_start' and 'global_end' are undefined")
                #  Save exons
                name = f"{prefix}_exons.fasta"
                out_path_exon = outdir / name
                global_start = overlapping_gp.global_start
                global_end = overlapping_gp.global_end
                exon_records: list[SeqRecord] = []
                super_records = [
                    self.scaffold_dict[scaffold_name]
                    for scaffold_name in overlapping_gp.extended_cds_dict
                ]
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
                                description="",
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
                if super_records:
                    SeqIO.write(super_records, out_path_supercontigs, "fasta")

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

        outdir.mkdir(parents=True, exist_ok=True)
        if scfs.stat().st_size > 0:
            olc = OverlappingCds(probes, scfs, treads=threads, **snakemake.params)

            ## Save scfs
            outdir.mkdir(parents=True, exist_ok=True)
            olc.save_records(outdir, snakemake.params.min_exonic_length)

            ## Save Probes
            olc.save_best_probe(outdir.parent.parent / f"probes/{probes.name}")


def main():
    # probe_dir = Path(
    #     '/home/yjkbertrand/Documents/projects/nextpiper/debug/NextPyper_hieracium/Merged_run/First_rna_targeted/results_full2/aster_kew_rna/homolog_prospection/region_separation/input_probes')
    # scaffolds_dir = Path(
    #     '/home/yjkbertrand/Documents/projects/nextpiper/debug/NextPyper_hieracium/Merged_run/First_rna_targeted/results_full2/aster_kew_rna/homolog_prospection/region_separation/input_scfs')
    # out_dir = Path(
    #     '/home/yjkbertrand/Documents/projects/nextpiper/debug/NextPyper_hieracium/Merged_run/First_rna_targeted/results_full2/aster_kew_rna/homolog_prospection/region_separation/output_scfs')
    # probe_dir = Path(
    #     "/home/yjkbertrand/Documents/projects/nextpiper/debug/centroids_noHMM2/homolog_prospection/region_separation/input_probes"
    # )
    # probe_dir = Path(
    #     "/home/yjkbertrand/Documents/projects/nextpiper/debug/centroids_noHMM2/homolog_prospection/region_separation/input_probes"
    # )
    #
    # scaffolds_dir = Path(
    #     "/home/yjkbertrand/Documents/projects/nextpiper/debug/centroids_noHMM2/homolog_prospection/region_separation/input_scfs"
    # )
    #
    # scaffolds_dir = Path(
    #     "/home/yjkbertrand/Documents/projects/nextpiper/debug/centroids_noHMM2/bug_instances/scaffolds"
    # )
    # scaffolds_dir = Path(
    #     "/home/yjkbertrand/Documents/projects/nextpiper/debug/centroids_noHMM2/homolog_prospection/region_separation/input_scfs"
    # )
    # out_dir = Path(
    #     "/home/yjkbertrand/Documents/projects/nextpiper/debug/centroids_noHMM2/homolog_prospection/test"
    # )
    # out_dir = Path(
    #     "/home/yjkbertrand/Documents/projects/nextpiper/debug/centroids_noHMM2/bug_instances/out_dir"
    # )
    # parameters = [8, 0.85, 0.1, 10, 0.7, "TIUZ_probe4471"]
    #
    # parameters = [8, 0.85, 0.1, 10, 0.7, "HLJG_probe5551"]
    # parameters = [8, 0.85, 0.1, 10, 0.7]
    #
    # problematic = ["5865.fasta"]
    # for scfs in scaffolds_dir.glob("*.fasta"):
    #     if scfs.name in problematic:
    #         continue
    #     # if not scfs.name == "5865.fasta":
    #     #     continue
    #     # if not scfs.name == "6221.fasta":
    #     #     continue
    #     # # probes_name = f"probes_{scfs.name.split('_')[1]}"
    #     # probes = probe_dir / probes_name
    #     probes = probe_dir / scfs.name
    #     matrix = "/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_paralogy_2/blosum62.csv"
    #     print(f"{scfs=},  {probes=}")
    #     olc = OverlappingCds(str(probes), str(scfs), matrix, *parameters)
    #     print("overlapping:", olc)
    #     # matrix = "/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_paralogy_2/blosum62.csv"
    #     # olc.paralogy_search(matrix)
    #     # print("saving")
    #     olc.save_records(out_dir, 10)

    # break

    # scfs = "/home/yjkbertrand/Documents/projects/nextpiper/debug/5899_test.fasta"
    # probes = probe_dir / "5899.fasta"
    # olc = OverlappingCds(str(probes), str(scfs), *parameters)
    # matrix = "/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_paralogy_2/blosum62.csv"
    # olc.paralogy_search(matrix)
    scfs = "/home/yjkbertrand/Documents/projects/nextpiper/test_data/temp/scfs_At4g32140.fasta"
    probe = "/home/yjkbertrand/Documents/projects/nextpiper/test_data/temp/probe_At4g32140.fasta"
    matrix = "/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_paralogy_2/blosum62.csv"
    parameters = [8, 0.85, 0.1, 10, 0.7]
    olc = OverlappingCds(str(probe), str(scfs), matrix, *parameters)
    out_dir = Path(
        "/home/yjkbertrand/Documents/projects/nextpiper/debug/centroids_noHMM2/bug_instances/out_dir"
    )
    olc.save_records(out_dir, 10)


if __name__ == "__main__":
    if "snakemake" in globals():
        snakemake_call(snakemake)
    else:
        main()
