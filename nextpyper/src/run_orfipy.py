#    Copyright (C) 2024
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.

"""
Functions and classes used to run orfipy and assemble the putative exons into the longest translated peptide.
For each probe, the longest path obtained by combining the possible orfs is calculated.
Stop codons ('*') are removed from the sequences and the orfs are combined.
AA sequences are saved in fasta format.
    Three parameters are use
    -sensitivity: between 0 and 1, is the proportion of the nucleotide sequence that gets translated to be deemed sufficient.
    -stop_per_1Kbp: proportion of stop codons within exons that are most likely sequencing errors
    -min_exon_length: minimum number of AA in an exon
#  Usage example:
    find_cds(Path("../data/probe.fasta"), Path("../data/longest_cds.fasta"), sensitivity, stop_per_1Kbp, min_exon_length)
"""

__version__ = "0.1"

# =======================================================================================
#               IMPORTS
# =======================================================================================
from collections import defaultdict, namedtuple
from dataclasses import dataclass, field, fields
import re
import sys
from typing import Final, Optional, Self, TypedDict, Literal, Any
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import orfipy_core

pattern = re.compile(
    r"""^ID=Seq_ORF\.(?P<name>\d+?);ORF_type=(?P<type>.*?);ORF_len=(?P<length>\d+?);ORF_frame=(?P<frame>.*?);Start.*$""",
    re.VERBOSE,
)

# Named tuple that keeps the nodes in a path either as their suffixes or their cds and the sum of their length.
GraphPath = namedtuple("GraphPath", ["path", "length"])

# =======================================================================================
#               EXCEPTIOMS
# =======================================================================================

# =======================================================================================
#               FUNCTIONS
# =======================================================================================


def select_best_cds(
    orf_list: list["Orf"], seq_length: int, max_stops: int, sensitivity: float
) -> GraphPath:
    """
    Create an interval graph for '+' and '-' list of orf. Explore all possible paths and determine the best path.
    Return the longest path that can be translated.
    Parameters
    ----------
    orf_list: list of Orf objects.
    seq_length: length in nucleotides of the original record.
    max_stops: maximum number of stops allowed on path.
    sensitivity: proportion of the nucleotide sequence that gets translated to be deemed sufficient.
    """
    partial_plus = []
    partial_minus = []
    for orf in orf_list:
        if orf.strand == "+":
            partial_plus.append(orf)
        else:
            partial_minus.append(orf)
    optima = None
    best_length = 0
    # perform the path search on each strand separately
    if partial_plus:
        sorted_partial_plus = sorted(partial_plus, key=lambda o: o.start)
        iter_graph_plus = ItervalGraph(
            sorted_partial_plus, seq_length, max_stops, sensitivity
        )
        optima = iter_graph_plus.get_best_path()
    if optima:
        best_length = optima.length
    if partial_minus:
        sorted_partial_minus = sorted(partial_minus, key=lambda o: o.start)
        # use the best length from the '+' strand
        iter_graph_minus = ItervalGraph(
            sorted_partial_minus, seq_length, max_stops, sensitivity, best_length
        )
        if (optimum := iter_graph_minus.get_best_path()) is not None:
            optima = optimum
    return optima


def find_cds(
    input_fasta: Path,
    output_fasta: Path,
    sensitivity: float = 0.95,
    stop_per_1Kbp: float = 2.0,
    min_exon_length: int = 20,
) -> None:
    """
    Controle function to run orfipan and use the identified orf to find the longest translated region.
    Parameters
    ----------
    input_fasta
    output_fasta
    sensitivity: proportion of the nucleotide sequence that gets translated to be deemed sufficient, between 0 and 1.
    stop_per_1Kbp: stops found within exons that are most likely sequencing errors.
    min_exon_length: discard the smallest exons.

    """
    assert input_fasta.exists(), f"{input_fasta} does not exist"
    Path(output_fasta.parent).mkdir(parents=True, exist_ok=True)
    new_records = []
    records = SeqIO.parse(input_fasta, "fasta")
    for record in records:
        print(f"working on record {record.id}")
        if (
            translated_record := AllOrfs(
                record, sensitivity, stop_per_1Kbp, min_exon_length
            ).get_best_pep()
        ) is not None:
            new_records.append(translated_record)
            trans_fraction = 100 * (len(translated_record.seq) * 3) / len(record.seq)
            print(
                f"orfipy translated {int(trans_fraction)} percent of the original sequence"
            )
        else:
            print("no translation occurred.")
        # break
    if new_records:
        SeqIO.write(new_records, output_fasta, "fasta")


# =======================================================================================
#               CLASSES
# =======================================================================================


@dataclass
class Orf:
    """
    Data structure used to store relevant data from a orfipy record.
    Attributes
    ----------
    -name: orfipy ID.
    -strand:
    -type: Currently not relevant.
    -start: start index on the nucleotide sequence (it would be the end on a rev-translated sequence).
    -end: end index on the nucleotide sequence (it would be the start on a rev-translated sequence).
    """

    name: str
    strand: Literal["+", "-"]
    type: Literal["complete", "3-prime-partial"]
    start: int
    end: int
    length: int
    frame: Literal[-1, -2, -3, 1, 2, 3]

    def invert_coordinates(self):
        """For reversed complemented sequences the start and end coordinates are reversed for better sorting"""
        left = self.start
        right = self.end
        self.start = right
        self.end = left

    def get_strand(self):
        return self.strand

    def get_length(self):
        return self.length

    def get_start(self):
        return self.start

    def get_end(self):
        return self.end

    def overlap(self, other: "Orf") -> bool:
        if self.get_end() > other.get_start():
            return True
        return False


@dataclass
class Node:
    name: str
    orf: Orf = field(repr=False)
    children: list[str] = field(init=False, default_factory=list)

    def overlap(self, other: "Node") -> bool:
        if self.orf.overlap(other.orf):
            return True
        return False

    def add_child(self, child: str) -> Self:
        self.children.append(child)
        return self

    def not_root(self) -> Self:
        self.root = False
        return self

    def has_children(self) -> bool:
        if self.children:
            return True
        return False

    def get_length(self) -> int:
        return self.orf.get_length()

    def get_children(self) -> list[str]:
        return self.children

    def get_orf(self) -> Orf:
        return self.orf


@dataclass
class ItervalGraph:
    """
    Data structure that create a graph of non overlapping intervals (nodes).
    Intervals are ordered by starting position. Starting with each target interval, we iterate over the rest
    of the intervals that start at greater position. For each interval, if it doesn't overlap with the target
    we draw an edge, otherwise we skip the interval. Once the graph is constructed, we iterarate over all nodes starting
    from the most 5' nodes. Finding the longest path of non-overlapping intervals uses a DSF from each starting node.
    Attributes
    ----------
    -orf_list: a list of Orf objects.
    -seq_length: length in nucleotides of the whole sequence.
    -sensitivity: proportion of the nucleotide sequence that gets translated to be deemed sufficient.
    -best_path_length: length of the best path found so far.

        Post Init
    -nodes: list of Node objects ordered by starting position.
    -node_dict: orf name as key and node as value.
    -best_graph_path: best path found so far.
    """

    orf_list: list["Orf"]
    seq_length: int
    max_stops: int
    sensitivity: float = field(default=0.95)
    best_path_length: int = field(default=0)
    nodes: list[Node] = field(init=False, default_factory=list)
    node_dict: dict[str, Node] = field(init=False, default_factory=dict)
    best_graph_path: GraphPath = field(init=False, default=None)

    def __post_init__(self):
        self._create_nodes()
        self._create_edges()
        self._explore_graph()

    def _create_nodes(self) -> Self:
        """
        Populate the nodes attribute.
        """
        sorted_orfs = sorted(self.orf_list, key=lambda x: x.start)
        for orf in sorted_orfs:
            node = Node(orf.name, orf)
            self.nodes.append(node)
        return self

    def _create_edges(self) -> Self:
        """Connect the node to all other nodes whose orf is left of the current node orf and whose sequences
        do not overlap"""
        self.node_dict = {node.name: node for node in self.nodes}
        for idx, node in enumerate(self.nodes):
            for other_node in self.nodes[idx + 1 :]:
                if node.overlap(other_node):
                    continue
                node.add_child(other_node.name)
        return self

    def _explore_graph(self) -> Self:
        """
        Find all paths between root nodes (starting node in the path) and leaf nodes (nodes without children).
        Populate best_graph_path, modify the best_path_length attribute.
        """

        def dfs_util(
            current_name: str,
            current_path=None,
            current_length: int = 0,
            number_stops=0,
        ) -> GraphPath:
            if current_path is None:
                current_path = []
            current_node = self.node_dict[current_name]
            current_end = current_node.orf.get_end()
            current_path_copy = current_path.copy()
            current_number_stop = number_stops
            current_path_copy.append(current_name)
            current_length += current_node.get_length()

            # case there is not enough distance from the current index to the end of the sequence
            # to discover a better path than the best current one
            if self.seq_length + current_length - current_end < self.best_path_length:
                return

            # case node is a leaf
            if not current_node.has_children():
                if current_length > self.best_path_length:
                    self.best_path_length = current_length
                    self.best_graph_path = GraphPath(
                        current_path_copy,
                        current_length,
                    )
                return
            else:
                for child in current_node.get_children():
                    child_node = self.node_dict[child]
                    child_start = child_node.orf.get_start()
                    child_frame = child_node.orf.frame
                    #  counting number of stops that separate adjacent exons.
                    if child_start - current_end <= 3:
                        current_number_stop += 1
                    if current_number_stop == self.max_stops:
                        if current_length > self.best_path_length:
                            self.best_path_length = current_length
                            self.best_graph_path = GraphPath(
                                current_path_copy,
                                current_length,
                            )
                        return
                    dfs_util(
                        child, current_path_copy, current_length, current_number_stop
                    )

        # check all nodes in the list, starting from 5' end
        for root_name in (x.name for x in self.nodes):
            root_node = self.node_dict[root_name]
            root_start = root_node.orf.get_start()
            # case there is not enough distance from the current index to the end of the sequence
            # to discover a better path than the best current one
            if self.seq_length - root_start < self.best_path_length:
                continue
            # case the path is good enough
            if self.best_path_length >= self.seq_length * self.sensitivity:
                break
            dfs_util(root_name)
        return self

    def get_best_path(self) -> Optional[GraphPath]:
        if not self.best_graph_path:
            return
        return self.best_graph_path


@dataclass
class AllOrfs:
    """
    Infer all possible orfs in both strands with orfipy.
    Select the longest orf or a path made out of orfs that provides the best translation.
    Attributes
    ----------
    -seq: Seq object from biopython.
    -sensitivity: proportion of the nucleotide sequence that gets translated to be deemed sufficient.
    -stop_per_1Kbp: stops found within exons that are most likely sequencing errors
    -min_exon_length: number of amino acids
        Post Init
    -candidates: keys are Orf numbers, value are Orf recovered with orfipy
    -best_path: best path found overall.
    """

    record: SeqRecord
    sensitivity: float = 0.99
    stop_per_1Kbp: float = 2.0
    min_exon_length: int = 20  # number of amino acids
    candidates: dict[str, Orf] = field(
        init=False,
        default_factory=dict,
    )
    best_path: GraphPath = field(init=False, default=None)  # name of orf(s)

    def __post_init__(self):
        if 0 < self.sensitivity <= 1:
            self.sensitivity = 0.95
        self._find_orfs()
        self._select_best()

    def _find_orfs(self):
        for start, stop, strand, description in orfipy_core.orfs(
            str(self.record.seq).upper(),
            strand="b",
            minlen=3,
            maxlen=20000,
            include_stop=False,
            partial3=True,
            partial5=True,
            between_stops=True,
            starts=["ATG"],
        ):
            match = pattern.match(description)
            name = match.group("name")
            length = int(match.group("length"))
            if length < self.min_exon_length * 3:
                continue
            frame_match = match.group("frame")
            frame = int(frame_match) if strand == "+" else -int(frame_match[1:])
            type_match = match.group("type")
            self.candidates[name] = Orf(
                name=name,
                strand=strand,
                type=type_match,
                start=start,
                end=stop,
                length=length,
                frame=frame,
            )

    def _select_best(self):
        length = len(self.record.seq)
        if not self.candidates:
            return
        sorted_orfs = sorted(
            self.candidates.values(), key=lambda o: o.length, reverse=True
        )
        if sorted_orfs[0].length >= length * self.sensitivity:
            best_orf_name = sorted_orfs[0].name
            best_orf_length = sorted_orfs[0].length
            self.best_path = GraphPath([best_orf_name], best_orf_length)
            return
        max_stop = int((length / 1000) * self.stop_per_1Kbp)
        self.best_path = select_best_cds(
            list(self.candidates.values()), length, max_stop, self.sensitivity
        )

    def get_best_pep(self) -> Optional[SeqRecord]:
        if not self.best_path:
            return
        new_sequence = []
        best_path = self.best_path.path
        strand = self.candidates[best_path[0]].strand
        if strand == "+":
            for orf_name in best_path:
                orf = self.candidates[orf_name]
                new_sequence.append(
                    str(self.record.seq[orf.start : orf.end].translate())
                )
        else:
            best_path.reverse()
            for orf_name in best_path:
                orf = self.candidates[orf_name]
                new_sequence.append(
                    str(
                        self.record.seq[orf.start : orf.end]
                        .reverse_complement()
                        .translate()
                    )
                )
        pep_seq = "".join(new_sequence).replace("X", "")
        return SeqRecord(
            Seq(pep_seq),
            id=self.record.id,
            description=self.record.description,
            name=self.record.name,
        )


def snakemake_call(snakemake) -> None:
    with open(snakemake.log[0], "w") as outlog:
        sys.stdout = sys.stderr = outlog

        probes = Path(snakemake.input[0])
        translated_probes = Path(snakemake.output[0])

        translated_prop = snakemake.params.translated_prop
        stop_per_1Kbp = snakemake.params.stop_per_1Kbp
        min_exon_length = snakemake.params.min_exon_length

        find_cds(
            probes, translated_probes, translated_prop, stop_per_1Kbp, min_exon_length
        )


if __name__ == "__main__":
    if "snakemake" in globals():
        snakemake_call(snakemake)
    else:
        # probes_fasta = Path("/home/yjkbertrand/Documents/projects/nextpiper/debug/batrachium_prones/targets.fasta")
        probes_fasta = Path(
            "/home/yjkbertrand/Documents/projects/nextpiper/debug/brassica_probes/kasper_mega_probes.fasta"
        )
        # output_fasta = Path("/home/yjkbertrand/Documents/projects/nextpiper/debug/batrachium_prones/test_translation.fasta")
        output_fasta = Path(
            "/home/yjkbertrand/Documents/projects/nextpiper/debug/brassica_probes/test_translation.fasta"
        )
        find_cds(probes_fasta, output_fasta)
