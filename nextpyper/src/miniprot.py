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
from collections import defaultdict, namedtuple, deque
from importlib import reload
from io import StringIO
from operator import attrgetter
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import textwrap
from typing import Optional, Self, NamedTuple

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq

# import gff_parser
# gff_parser = reload(gff_parser)
from gff_parser import Fragment, Cds
from interval_tree import IntervalST, Interval, Node

# # Parse the header of SAUTE scaffolds
# saute_pattern = re.compile(
#     r"^Contig_(?P<name>.*?):(?P<component>\d+?):[^ ]+$"
# )

Region = namedtuple("Region", ["start", "end"])

# =======================================================================================
#               FUNCTIONS
# =======================================================================================
#
# def get_nuc_coordinates(
#         probe_start: int, probe_end: int, fragments: list[Fragment]
# ) -> list[Optional[int]]:
#     """
#     Given a start and end position on a probe in AA space, return the nucleotide coordinates of start and end.
#     :param probe_start:
#     :param probe_end:
#     :param fragments: list of fragments
#     :return: [None, None] if for some reason there are no AA-nt correspondence, otherwise [start, end]
#     """
#     correspondences = {}
#     for fragment in fragments:
#         correspondences.update(fragment.get_correspondence())
#     while probe_start < probe_end:
#         if (first := correspondences.get(probe_start)) is None:
#             probe_start += 1
#         if (last := correspondences.get(probe_end)) is None:
#             probe_end -= 1
#         if None not in [first, last]:
#             return sorted(
#                 [
#                     first,
#                     last,
#                 ]
#             )
#     return [None, None]


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
        if (line.startswith("##ATN")
                or line.startswith("##ATA")
                or line.startswith("##AAS")
                or line.startswith("##AQA")):
            new_lines.append(line)
    return "".join(new_lines).encode()


def run_miniprot(
        probe_path: Path,
        contig_path: Path,
        treads=2,
        min_similarity=0.85,
        min_coverage=0.01,
) -> StringIO | str:
    """
    Wrapper for running miniprot.
    :return:
    """
    miniprot_cmd = f"miniprot -t {treads} --gff --aln --outn 1 -J 50 --outs {min_similarity} --outc {min_coverage} {contig_path} {probe_path}".split()
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
    # print(miniprot.stdout)
    return StringIO(miniprot.stdout)


def run_miniprot_boundary_scorer(miniprot: bytes, boundary_scorer_out: Path, matrix_path: Path):
    """
    Wrapper for running miniprot_boundary_scorer.
    :return:
    """
    boundary_scorer_cmd = f"miniprot_boundary_scorer -o {boundary_scorer_out} -s {matrix_path}"
    try:
        subprocess.run(
            boundary_scorer_cmd,
            timeout=10,
            shell=True,
            input=miniprot, )
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

    def get_mean_score(self, nbr_contigs:int):
        """The  score used for the global ranking is based on probe positions, the total number of
        scaffolds and the scaffolds without a match. Missing matches get a penalty equal to the maximum rank,
         e.g. rank_score=[0, 3, 1, 0], coverage=['contig_0', 'contig_1', 'contig_2', 'contig_3']
         and nbr_contigs=6  produces a score of (4+6+6)/6"""
        missing_contigs = nbr_contigs - len(self.coverage)
        missing_scores = [nbr_contigs for _ in range(missing_contigs)]
        return (sum(self.rank_score) + sum(missing_scores)) / nbr_contigs


@dataclass
class Exon:
    __slots__ = ['start', 'end']
    start: int
    end: int


@dataclass
class Exon_correspondence:
    """
    Probe to scaffold coordinates
    """
    __slots__ = ['exon_probe', 'exon_scaffold', 'length_on_scaffold']
    exon_probe: Exon
    exon_scaffold: Exon
    length_on_scaffold: int

    def check_identity(self, other:Exon, sensitivity: int)->bool:
        """
        Given another probe exon coordinates, check if this exon is compatible:
        whether it is the same set of coordinates or whether this exon comes from
        the fusion of two exons and the query is one of them.
        Parameters
        ----------
        other
        sensitivity: the range around the query coordinates that should be searched for a match.

        Returns
        -------

        """
        probe_start = self.exon_probe.start
        probe_end = self.exon_probe.end
        other_probe_start = other.start
        other_probe_end = other.end
        print(f"{probe_start=},{probe_end=},{other_probe_start=},{other_probe_end=}")
        if self.exon_probe == other:
            return True
        # case two exon fusion left part
        if ((probe_start - sensitivity <= other_probe_start <= probe_start + sensitivity) and (
                other_probe_end - sensitivity <= probe_end )):
            return True
        # case two exon fusion right part
        if ((probe_start <= other_probe_start + sensitivity)
                and (probe_end - sensitivity <= other_probe_end <= probe_end + sensitivity)):
            return True
        return False
    # exon discovery score?


@dataclass
class Boundary:
    """
    coordinate on scaffole
    """
    __slots__ = ['scaffold_name', 'scaffold_start', 'scaffold_end']
    scaffold_name: str
    scaffold_start: int
    scaffold_end: int


@dataclass
class ParalogyCds(Cds):
    """
    Cds object with additional features to handle paralogy information
    """
    exon_probe_scaffold: list = field(init=False, default_factory=list)
    correspondence: dict[int, int] = field(default_factory=dict, init=False, repr=True)
    rev_correspondence: dict[int, int] = field(default_factory=dict, init=False, repr=True)
    exon_correspondences: list[Exon_correspondence] = field(default_factory=list, init=False, repr=True)
    accepted_exons: list[Exon] = field(default_factory=list, init=False, repr=True)

    def __post_init__(self):
        super().__post_init__()
        self._reverse_correspondence()

    def _reverse_correspondence(self) -> Self:
        for fragment in self.fragments:
            for k, v in fragment.correspondence.items():
                self.correspondence[k] = v
                self.rev_correspondence[v] = k
        self.start_on_probe = min(self.rev_correspondence)
        self.end_on_probe = max(self.rev_correspondence)
        return self

    def find_probe_exons(self, exons: list[Exon]) -> Self:
        """
        Find coordinates pairs on the probe where exons are located
        """
        for scaff_exon in exons:
            probe_start = None
            probe_end = None
            for idx in [0, 1, 2, -1, -2]:
                if (probe_start := self.rev_correspondence.get(scaff_exon.start + idx)) is not None:
                    break
            for idx in [0, 1, 2, -1, -2]:
                if (probe_end := self.rev_correspondence.get(scaff_exon.end + idx)) is not None:
                    break
            if None not in [probe_start, probe_end]:
                probe_exon = Exon(probe_start, probe_end)
                length = scaff_exon.end - scaff_exon.start
                self.exon_correspondences.append(Exon_correspondence(probe_exon, scaff_exon, length))
                # print(
            #     f"{scaff_exon.start=}\t{scaff_exon.end=}, {probe_start=}, {probe_end=}, length={scaff_exon.end - scaff_exon.start}")
        return self

    def find_scaffold_exon(self, probe_exon:Exon, sensitivity: int=10) -> Optional[Exon]:
        """
        Given a pair of coordinates on the probe, retrieve the coordinates on the scaffold
        Parameters
        ----------
        probe_exon
        sensitivity

        Returns
        -------

        """
        if sensitivity < 2:
            sensitivity = 2
        if sensitivity > 10:
            sensitivity = 10
        search_range = [0]
        for (i, j) in zip(range(1, sensitivity), range(-1,-sensitivity,-1)):
            search_range.extend([i,j])
        scaffold_start = None
        scaffold_end = None
        for idx in search_range:
            if (scaffold_start  := self.correspondence.get(probe_exon.start + idx)) is not None:
                break
        for idx in search_range:
            if (scaffold_end := self.correspondence.get(probe_exon.end + idx)) is not None:
                break
        print(f"{scaffold_start=}, {scaffold_end=}")
        if None in [scaffold_start, scaffold_end]:
            return
        return Exon(scaffold_start, scaffold_end)

    def find_scaffold_intron(self, probe_start, probe_end, sensitivity: int=10) -> Optional[Exon]:
        if sensitivity < 2:
            sensitivity = 2
        if sensitivity > 10:
            sensitivity = 10
        search_range = [0]
        for (i, j) in zip(range(1, sensitivity), range(-1,-sensitivity,-1)):
            search_range.extend([i,j])


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
    cds_dict: dict[str, ParalogyCds] = field(init=False, repr=False, default_factory=dict)
    common_exons: list[Node] = field(init=False, repr=False, default_factory=list)
    paralogs: list[str] = field(init=False, repr=False, default_factory=list)
    orthologs: list[str] = field(init=False, repr=False, default_factory=list)
    flank_left: list[Boundary] = field(init=False, repr=False, default_factory=list)
    flank_right: list[Boundary] = field(init=False, repr=False, default_factory=list)
    exons: defaultdict[Region, list[Boundary]] = field(default_factory=lambda: defaultdict(list))
    introns: defaultdict[int, list[Boundary]] = field(default_factory=lambda: defaultdict(list))

    def find_common_exons(self, min_proportion=0.5, min_overlap=0.9):
        """
        min_proportion: minimum proportion of samples (that have a probe match) that have the exon.
            If a scaffold is too short, it is counted in the proportion anyway.
        min_overlap: min overlap between a target exon and the valid exons.
        """
        interval_tree_all_exons = IntervalST()
        for contig in self.cds_dict:
            cds = self.cds_dict[contig]
            if cds.is_empty():
                continue

            for exon_pair in cds.exon_correspondences:
                interval_tree_all_exons.put(Interval(exon_pair.exon_probe.start, exon_pair.exon_probe.end), contig)

        nodes = interval_tree_all_exons.tree_traversal_bsf_with_values()
        print("nodes", nodes)
        valid_nodes = []
        matching_samples = [contig for contig in self.cds_dict if not self.cds_dict[contig].is_empty()]
        # Check for all samples that are too short if they could have the exon
        for node in nodes:
            samples = node.value
            interval = node.interval
            possible_samples = len(samples)  # samples with the exon plus samples that are not full length
            missing_samples = set(matching_samples) - set(samples)
            # to the putative exons, add samples that are too short to cover it
            for sample in missing_samples:
                sample_cds = self.cds_dict[sample]
                sample_start = sample_cds.probe_start
                sample_end = sample_cds.probe_end
                if interval.hi <= sample_start or sample_end <= interval.lo:
                    possible_samples += 1
                    node.value.append(sample)
            if possible_samples / len(matching_samples) > min_proportion:
                valid_nodes.append(node)
        print("valid nodes", valid_nodes)
        self.common_exons = valid_nodes

    def eliminate_paralogs(self, overlap: int):
        if not self.common_exons:
            print("no valid exon")
            # add all sequences to orthologs
            self.paralogs.extend(self.cds_dict)
            return
        for scaffold, cds in self.cds_dict.items():
            presence = []
            valid_exons = self.common_exons
            for exon in cds.exon_correspondences:
                print(f"{exon=}")
                found_start = False
                found_end = False
                exon_start, exon_end = exon.exon_probe.start, exon.exon_probe.end
                while valid_exons:
                    target_exon = valid_exons[0]
                    print(f"{target_exon=}")
                    if scaffold in target_exon.value:
                        presence.append(True)
                        break
                    probe_start = target_exon.interval.lo
                    probe_end = target_exon.interval.hi
                    if max(0, exon_start - overlap) <= probe_start < exon_start + overlap:
                        found_start = True
                    if max(0, exon_end - overlap) <= probe_end < exon_end + overlap:
                        found_end = True
                    valid_exons = valid_exons[1:]
                    if not valid_exons and not all([found_start, found_end]):
                        presence.append(False)
                        break
                    if found_start and found_end:
                        presence.append(True)
                        break
                    elif not found_start and not found_end:
                        presence.append(False)
                        break
            print(f"{presence}")
            if all(presence):
                self.orthologs.append(scaffold)
            else:
                self.paralogs.append(scaffold)
            # print(f"{scaffold=}\t{presence=}")
        print(f"{self.orthologs=}")
        print(f"{self.paralogs=}")

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

    def chop_sequences(self):
        """
        Cut the different sequences into introns/exons/flanking
        Returns
        -------

        """
        #  Deal with the first exon and flanking left region
        print("chop chop")
        probe_exons = deque(self.common_exons)
        idx = 0
        previous_exon = None
        fully_covered_scaffolds = [] # scaffolds that have been fully chopped
        while probe_exons:
            current_exon = probe_exons.popleft()
            print(f"{current_exon=}")  # Node object
            exon_start, exon_end = current_exon.interval.lo, current_exon.interval.hi
            probe_exon = Exon(exon_start, exon_end)
            current_region = Region(exon_start, exon_end)
            print(f"{current_region=}")
            for scaffold in self.orthologs:
                # if scaffold != "Corymbium_cymosum_SRR6072274-7021_0_226":
                #     continue
                print(f"{scaffold=}")

                paracds = self.cds_dict[scaffold]
                print(paracds)
                print(f"{idx=}")
                current_accepted_exons = paracds.accepted_exons.copy()
                for correspondence in paracds.exon_correspondences:
                    if (has_exon:= correspondence.check_identity(probe_exon, 10)) is True:
                        print(f"correspondence {correspondence=}")
                        break
                else:
                    print(f"no correspondence found for probe exon {probe_exon=} in scaffold {scaffold=}")
                #  Case the exon is in the scaffold
                if not has_exon:
                    # Check that we have reached the end of the sequence
                    if not previous_exon:
                        continue
                    if scaffold not in fully_covered_scaffolds and paracds.accepted_exons:
                        exon_left = paracds.accepted_exons[-1].end
                        self.introns[previous_exon.end +1].append(Boundary(scaffold, exon_left, -1))
                        fully_covered_scaffolds.append(scaffold)
                        print(f"for scaffold {scaffold=} reaching the end for exon {previous_exon.end}, {probe_exon.end}")
                    continue

                scaffold_exon = paracds.find_scaffold_exon(probe_exon)
                print(f"found one matching exon: {scaffold_exon=}")
                if not scaffold_exon:
                    sys.exit(f"no scaffold exon found for probe exon {probe_exon=} in scaffold {scaffold=}")
                scaffold_start, scaffold_end = scaffold_exon.start, scaffold_exon.end+3
                print(f"adding exon {scaffold_start=}, {scaffold_end=}")
                self.exons[current_region].append(Boundary(scaffold, scaffold_start, scaffold_end))
                paracds.accepted_exons.append(Exon(scaffold_start, scaffold_end))
                print(
                    f"Case the exon is in the scaffold, scaffold boundaries are {scaffold_start=}, {scaffold_end=}")
                if idx ==0:
                    print("adding to left flank")
                    if scaffold_start != 0:
                        self.flank_left.append(Boundary(scaffold, 0, scaffold_start))

                if idx == len(self.common_exons)-1:
                    self.flank_right.append(Boundary(scaffold, scaffold_end, -1))
                    print("adding to right flank")

                #  Add intron information
                if current_accepted_exons:
                    print(f"{current_accepted_exons=}")
                    exon_left= current_accepted_exons[-1].end
                    #  Case no intron
                    if exon_left== scaffold_start:
                        continue
                    #  Case regular intron
                    print(f"assigning intron {exon_left =}, {scaffold_start=}")
                    self.introns[exon_start].append(Boundary(scaffold, exon_left , scaffold_start))
                    continue

                #  Case left truncated scaffold
                if previous_exon and not current_accepted_exons:
                    self.introns[exon_start].append(Boundary(scaffold, 0, scaffold_start))
                    continue
            idx += 1
            previous_exon = probe_exon
        print("introns",list(self.introns.items()))
        print("exons", list(self.exons.items()))
        print("flanking right,",list(self.flank_right))
        print("flanking left,",list(self.flank_left))

    def __repr__(self):
        return (f"probe_start={self.probe_start}, probe_end={self.probe_end},"
                f"{sorted(self.seq_names)}")

class TrimmedSeqs(NamedTuple):
    flank_left: list[SeqRecord]
    probe_match: list[SeqRecord]
    flank_right: list[SeqRecord]


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
    -min_global_identity: identity over all several fragments. Should be slightly lower than min_probe_contig_sim.

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
    min_global_identity: float = field(default=0.00)
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

    def _run_miniprot(self, record, probe_path_temp, contig_path_temp, paralogy=False) -> Cds:
        SeqIO.write(record, contig_path_temp, "fasta")
        miniprot_out = run_miniprot(
            probe_path_temp,
            contig_path_temp,
            self.treads,
            self.min_probe_contig_sim,
            self.min_fragment_cov,
        )
        if paralogy:
            return miniprot_out
        return Cds(miniprot_out)


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
    user_probe: str = field(default=None)
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
        print(f"{self.user_probe=}")
        with tempfile.TemporaryDirectory() as tmpdirname:
            tmp_probe_paths = {}
            # write all probe records in separate file
            if self.user_probe is not None and self.user_probe in self.probes_dict.keys():
                print(f"User selected probe is {self.user_probe}")
                fasta_name = f"{self.user_probe}.fas"
                contig_fasta = Path(tmpdirname) / fasta_name
                SeqIO.write(self.probes_dict[self.user_probe], contig_fasta, "fasta")
                tmp_probe_paths[self.user_probe] = contig_fasta

            else:
                print("creating probe list")
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
                    if not cds.is_empty() and cds.global_identity>=self.min_global_identity:
                        self.cds_dict[contig].append(cds)
                print("cds", self.cds_dict.get(contig))
            print("computing rank score")
            self._compute_rank_score()
            print("initial sequences:", len(self.cds_dict))
            print("filtered sequences:", len(self.filtered_cds_dict))
            print("sorting overlapping")
            self._sorting_overlapping()
            print("done sorting overlapping")
            #print("non_overlapping", self.non_overlapping)

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
            probe_ranks.items(), key=lambda x: x[1].get_mean_score(len(self.cds_dict))
        )[0]

        self.best_probe = best_combination[0]
        print(f"Best overall probe: {self.best_probe}")
        # Keep only contigs that are covered by the best probe
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
        interval_list = []
        for contig, cds in self.filtered_cds_dict.items():
            print(f"{contig=}, {cds=}")
            if cds.is_empty():
                continue
            start = cds.probe_start
            end = cds.probe_end
            interval_list.append(OverlappingSeqs(start, end, [contig]))
        self.non_overlapping = merge_intervals(interval_list)
        return self

    def _boundary_scorer_parser(self, boundary_scorer_path: Path) -> list[Exon]:
        """
        Parser of the miniprot_boundary_scorer output.
        """
        lines = boundary_scorer_path.read_text().split("\n")
        exons = []
        for line in lines:
            if not line:
                continue
            splt = line.split("\t")
            if splt[2] == 'CDS' and splt[6] == "+":
                exons.append(Exon(int(splt[3]) - 1, int(splt[4]) - 3))
        return exons

    def paralogy_search(self, substitution_matrix: Path):
        """
        On each OverlappingSeqs object perform a matching exon search.
        The output of miniprot is converted into a format that run_miniprot_boundary_scorer can parse.

        """
        if not self.non_overlapping:
            print("No overlapping group")
            return
        with tempfile.TemporaryDirectory() as tmpdirname:
            probe_record = self.probes_dict[self.best_probe]
            probe_name = f"{self.best_probe}.fas"
            probe_path = Path(tmpdirname) / probe_name
            SeqIO.write(probe_record, probe_path, "fasta")

            for overlapping_gp in self.non_overlapping:
                #  Run miniprot on each contig on each group of overlapping sequences
                for contig in overlapping_gp.seq_names:
                    record = self.contigs_dict[contig]
                    #print(f"working on paralogy contig {contig}")
                    contig_path = Path(tmpdirname) / f"{contig}.fas"
                    miniprot_out = self._run_miniprot(record, probe_path, contig_path, paralogy=True)
                    cds = ParalogyCds(miniprot_out)
                    cds.length = len(SeqIO.read(contig_path, "fasta").seq)
                    #print(contig, cds)
                    miniprot_out.seek(0)
                    boundary_scorer_out = Path(tmpdirname) / f"{contig_path.stem}_boundary.gfa"
                    # run miniprot_boundary scorer
                    run_miniprot_boundary_scorer(remove_gff(miniprot_out), boundary_scorer_out, substitution_matrix)
                    #print(self._boundary_scorer_parser(boundary_scorer_out))
                    cds.find_probe_exons(self._boundary_scorer_parser(boundary_scorer_out))
                    overlapping_gp.cds_dict[contig] = cds
                print(list(overlapping_gp.cds_dict.items()))
                overlapping_gp.find_common_exons()
                print("eliminate paralogs")
                #overlapping_gp.eliminate_paralogs(15)
                #overlapping_gp.chop_sequences()

    def save_records(self):
        import os
        os.chdir("/home/yjkbertrand/Documents/projects/nextpiper/debug/chopped_seqs")
        stem_name = self.probes_path.stem
        print(f"saving {stem_name }")
        
        for overlapping_gp in self.non_overlapping:
            prefix = f"{stem_name}_{overlapping_gp.probe_start}_{overlapping_gp.probe_end:}"
            if overlapping_gp.exons:
                for exon, boundaries in overlapping_gp.exons.items():
                    name = f"{prefix}_exon_{exon.start}_{exon.end}.fasta"
                    records = []
                    for boundary in boundaries:
                        print(boundary)
                        scaffold_name = boundary.scaffold_name
                        start = boundary.scaffold_start
                        end = boundary.scaffold_end
                        print(f"{self.contigs_dict[scaffold_name]=}")
                        chopped_seq = str(self.contigs_dict[scaffold_name].seq)[start:end]
                        new_record = SeqRecord(Seq(chopped_seq), name="", description="", id=scaffold_name)
                        print(f"new record: {new_record}")
                        records.append(new_record)
                    SeqIO.write(records, name, "fasta")
            if overlapping_gp.introns:
                for intron, boundaries in overlapping_gp.introns.items():
                    name = f"{prefix}_intron_{intron}.fasta"
                    records = []
                    for boundary in boundaries:
                        print(boundary)
                        scaffold_name = boundary.scaffold_name
                        start = boundary.scaffold_start
                        end = boundary.scaffold_end
                        print(f"{self.contigs_dict[scaffold_name]=}")
                        chopped_seq = str(self.contigs_dict[scaffold_name].seq)[start:end]
                        new_record = SeqRecord(Seq(chopped_seq), name="", description="", id=scaffold_name)
                        print(f"new record: {new_record}")
                        records.append(new_record)
                    SeqIO.write(records, name, "fasta")
                    
            if overlapping_gp.flank_left:
                left_flank_records = []
                for boundary in overlapping_gp.flank_left:
                    scaffold_name = boundary.scaffold_name
                    start = boundary.scaffold_start
                    end = boundary.scaffold_end
                    chopped_seq = str(self.contigs_dict[scaffold_name].seq)[start:end]
                    new_record = SeqRecord(Seq(chopped_seq), name="", description="", id=scaffold_name)
                    left_flank_records.append(new_record)
                SeqIO.write(left_flank_records,  f"{prefix}_left_flank.fasta", "fasta")
            if overlapping_gp.flank_right:
                right_flank_records = []
                for boundary in overlapping_gp.flank_right:
                    scaffold_name = boundary.scaffold_name
                    start = boundary.scaffold_start
                    seq = str(self.contigs_dict[scaffold_name].seq)
                    if start == len(seq):
                        continue
                    chopped_seq = seq[start:]
                    new_record = SeqRecord(Seq(chopped_seq), name="", description="", id=scaffold_name)
                    right_flank_records.append(new_record)
                SeqIO.write(right_flank_records, f"{prefix}_right_flank.fasta", "fasta")
            
            if overlapping_gp.paralogs:
                paralog_records = [self.contigs_dict[para] for para in overlapping_gp.paralogs]
                SeqIO.write(paralog_records, f"{prefix}_paralogs.fasta", "fasta")


    def save_best_probe(self, fasta_file: str) -> None:
        fasta_path = Path(fasta_file)
        fasta_path.parent.mkdir(parents=True, exist_ok=True)
        print(f"Saving best probe to fasta file {fasta_file}")
        SeqIO.write([self.probes_dict[self.best_probe]], fasta_path, "fasta")



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


def main():
    # probe_dir = Path(
    #     '/home/yjkbertrand/Documents/projects/nextpiper/debug/NextPyper_hieracium/Merged_run/First_rna_targeted/results_full2/aster_kew_rna/homolog_prospection/region_separation/input_probes')
    # scaffolds_dir = Path(
    #     '/home/yjkbertrand/Documents/projects/nextpiper/debug/NextPyper_hieracium/Merged_run/First_rna_targeted/results_full2/aster_kew_rna/homolog_prospection/region_separation/input_scfs')
    # out_dir = Path(
    #     '/home/yjkbertrand/Documents/projects/nextpiper/debug/NextPyper_hieracium/Merged_run/First_rna_targeted/results_full2/aster_kew_rna/homolog_prospection/region_separation/output_scfs')
    probe_dir = Path('/home/yjkbertrand/Documents/projects/nextpiper/debug/centroids_noHMM2/homolog_prospection/region_separation/input_probes')
    scaffolds_dir = Path('/home/yjkbertrand/Documents/projects/nextpiper/debug/centroids_noHMM2/homolog_prospection/region_separation/input_scfs')
    out_dir = Path('/home/yjkbertrand/Documents/projects/nextpiper/debug/centroids_noHMM2/homolog_prospection/test')
    parameters = [8, 0.85, 0.1, 10, 0.7,'TIUZ_probe4471']
    #parameters = [8, 0.85, 0.1, 10, 0.7,]
    for scfs in scaffolds_dir.glob("*.fasta"):
        if not scfs.name == "4471.fasta":
            continue
        probes = probe_dir / scfs.name
        print(f"{scfs=},  {probes=}")
        olc = OverlappingCds(str(probes), str(scfs), *parameters)
        print("overlapping:", olc)
        matrix = "/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_paralogy_2/blosum62.csv"
        olc.paralogy_search(matrix)
        print("saving")
        #olc.save_records()
        # olc.save_scaffolds(out_dir)
        break

    # scfs = "/home/yjkbertrand/Documents/projects/nextpiper/debug/5899_test.fasta"
    # probes = probe_dir / "5899.fasta"
    # olc = OverlappingCds(str(probes), str(scfs), *parameters)
    # matrix = "/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_paralogy_2/blosum62.csv"
    # olc.paralogy_search(matrix)

if __name__ == "__main__":
    if "snakemake" in globals():
        snakemake_call(snakemake)
    else:
        main()
