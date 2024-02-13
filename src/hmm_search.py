#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2024
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.

"""
Classes used to perform HMM search of hmm profiles against nucleotide sequences.
The sequences are derived from the nodes of a blunted and compacted assembly graph.
For some reason the regular pyhmmer pipeline does not search both strands, so I am using the LongTargetsPipeline.
RC sequences are indicated with the target_to and target_from in reversed order.

Example usage:

probe_fasta = Path("probe_3_aa.fasta")
profile_hmm = Path("profile.hmm")
hmm = Hmmer_result(probe_fasta, profile_hmm)
for k, v in hmm.get_node_hits.items():
    print(k, v)
Node_8 {'probe_3_aa_0': Profile_hits(hmm_name='probe_3_aa_0', domain_hits=[Domain(profile_start=60, profile_end=452, node_start=343, node_end=1, e_value=7.85e-102)]),
    'probe_3_aa_10': Profile_hits(hmm_name='probe_3_aa_10', domain_hits=[Domain(profile_start=120, profile_end=471, node_start=317, node_end=2, e_value=1.416e-67)])}
Node_9 {'probe_3_aa_0': Profile_hits(hmm_name='probe_3_aa_0', domain_hits=[Domain(profile_start=453, profile_end=579, node_start=1, node_end=127, e_value=6.11e-46)]),
    'probe_3_aa_10': Profile_hits(hmm_name='probe_3_aa_10', domain_hits=[Domain(profile_start=473, profile_end=598, node_start=1, node_end=126, e_value=1.01e-31)])}


"""
from collections import defaultdict, namedtuple
from dataclasses import dataclass, field
from pathlib import Path
from typing import Final, Optional, Self, TypedDict, Literal, Any

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pyhmmer
from pyhmmer.easel import SequenceFile
from pyhmmer.plan7 import HMMFile, HMM, Pipeline, LongTargetsPipeline

ALPHABET = pyhmmer.easel.Alphabet.dna()
BACKGROUND = pyhmmer.plan7.Background(ALPHABET)

# object that contains the result of a hmm search of a single profile against a single sequence.
Domain = namedtuple(
    "Domain",
    ["profile_start", "profile_end", "node_start", "node_end", "score", "e_value"],
)


@dataclass
class Node_hit:
    """
    A dataclass to store the results of a single HMMer's search.
    Unlike the Profile_hit dataclass, this class is used to store the results of hmmer's search
    per node and per probe
    Attributes
    ----------
        -node_name: name of the contig that has a hmm hit.
        -profile_name: name of the probe that hits.
        -domain_hits: list of regions (Domains) on the node that have a hmm hit.
        -domain: the regions in domain_hits are concatenated to keep only the beginning and end coordinates.
    """

    node_name: str
    profile_name: str
    domain_hits: list[Domain] = field(default_factory=list, repr=False, init=False)
    domain: Domain = field(init=False)

    def concatenate_domains(self) -> Self:
        """
        When a single node has several hmm hits, concatenate all Domains and create a new Domain object
        with min start and max end coordinates. The score and e_value attribute are averaged over all domains.
        :return: populates the domain object.
        """
        if len(self.domain_hits) == 1:
            self.domain = self.domain_hits[0]
        else:
            new_profile_start = self.domain_hits[0].profile_start
            new_profile_end = self.domain_hits[-1].profile_end
            new_node_start = self.domain_hits[0].node_start
            new_node_end = self.domain_hits[-1].node_end
            new_score = sum([x.score for x in self.domain_hits]) / len(self.domain_hits)
            new_e_value = sum([x.e_value for x in self.domain_hits]) / len(
                self.domain_hits
            )
            self.domain = Domain(
                new_profile_start,
                new_profile_end,
                new_node_start,
                new_node_end,
                new_score,
                new_e_value,
            )
        return self

    def invert_RC(self) -> Self:
        self.domain = Domain(
            self.domain.profile_start,
            self.domain.profile_end,
            self.domain.node_end,
            self.domain.node_start,
            self.domain.score,
            self.domain.e_value,
        )
        return self


@dataclass
class Path_nodes:
    node_names: tuple[str]
    node_hits: list[Node_hit]
    score: int


@dataclass
class Profile_hits:
    """
    A dataclass to store the results of a single HMMer's search.
    Attributes
    ----------
        -hmm_name: the name of the profile matched by the sequence.
        -domain_hits: a list of regions where the profile and the sequence match.
    """

    hmm_name: str
    domain_hits: list[Domain] = field(default_factory=list, repr=False)
    # domain: Domain = field(init=False)


#
# @dataclass
# class Best_path:
#     path_nodes: list[str] = field(default_factory=list)
#     path_profiles: list[Profile_hits] = field(default_factory=list)
#     mean_score: float = field(init=False)
#
#     def __post_init__(self):
#         self.get_mean_score()
#
#     def get_mean_score(self) -> Self:
#         self.mean_score = sum([x.domain.score for x in self.path_profiles]) / len(
#             self.path_profiles
#         )
#         return self


@dataclass
class Hmmer_result:
    """
    A dataclass to store the results of a HMMer's seach of multiple profiles against multiple sequences.
    Attributes
    ----------
        -sequences_fasta: path to the node sequences.
        -hmm_file: path to the precomputed hmmer profiles.
        -max_evalue: maximum expected e-value for the hmmer match on a single domain.
        -best_path: dict with edge names as keys Path_nodes objects as values.
        -node_hits: dictionary with node names as keys and a dictionary with evidence as the values,
            where each key is a profile name and each value a concatenated domain.
    """

    sequences_fasta: Path
    hmm_file: Path
    max_evalue: float = field(default=1.0e-20)
    best_paths: dict[str, Path_nodes] = field(default_factory=dict, repr=False)
    node_hits: defaultdict[dict] = field(default_factory=lambda: defaultdict(dict))

    def __post_init__(self):
        self._search_hmm()
        self._create_node_hits()

    def get_node_hits(self) -> defaultdict[dict]:
        return self.node_hits

    def _create_node_hits(self) -> Self:
        """
        Convert the elements of best_paths that are profile centered into elements
        of node_hits that are node centered.
        :return: populate the node_hits class attribute.
        """
        for path in self.best_paths.values():
            for hit in path.node_hits:
                self.node_hits[hit.node_name][hit.profile_name] = hit
        return self

    def _search_hmm(self) -> Self:
        """
        Searching each profile against the node sequences.
        This search is performed on both strands of the sequence which makes it more efficient than
        the scan_seq approach that requires both strands (direct and reverse complemented) to be generated.
        The nodes that are hit by a profile are stored in a pseudo-path.
        When two pseudo-paths are identical, the one with the best score is stored.
        :return: populate the best_paths class attribute.
        """
        with SequenceFile(
            self.sequences_fasta, digital=True, alphabet=ALPHABET, format="fasta"
        ) as seq_file:  # format 'afa' is an aligned fasta!
            sequences = seq_file.read_block()
            pipeline = LongTargetsPipeline(ALPHABET, background=BACKGROUND)

            with HMMFile(self.hmm_file) as hmm_record:
                while True:
                    try:
                        # iterate over all hmms in profile
                        hmm = next(hmm_record)
                        hits = pipeline.search_hmm(query=hmm, sequences=sequences)
                        profile_name = hits.query_name.decode()
                        all_hits = []
                        if hits:
                            for hit in hits:
                                node_name = hit.name.decode()
                                node_hit = Node_hit(node_name, profile_name)
                                for domain in hit.domains:
                                    if domain.i_evalue > self.max_evalue:
                                        continue
                                    domain = Domain(
                                        domain.alignment.hmm_from,
                                        domain.alignment.hmm_to,
                                        domain.alignment.target_from,
                                        domain.alignment.target_to,
                                        domain.score,
                                        domain.i_evalue,
                                    )
                                    node_hit.domain_hits.append(domain)
                                if node_hit.domain_hits:
                                    node_hit.concatenate_domains()
                                    all_hits.append(node_hit)
                            if all_hits:
                                total_score = sum(
                                    [node.domain.score for node in all_hits]
                                )
                                node_names = tuple(
                                    sorted([node.node_name for node in all_hits])
                                )
                                path = Path_nodes(node_names, all_hits, total_score)
                                # Check that the pseudo-path has not already been found
                                if node_names in self.best_paths.keys():
                                    if self.best_paths[node_names].score > total_score:
                                        pass
                                    else:
                                        self.best_paths[node_names] = path
                                else:
                                    self.best_paths[node_names] = path

                    except StopIteration:
                        break
        return self


def extract_seq_from_gfa(gfa_path: Path, fasta: str):
    """
    Take a GFA graph and extract the node sequences into a fasta file.
    :param gfa: input GFA file.
    :param fasta: output fasta file save on the graph's path.
    :return: creates a fasta file.
    """
    records = []
    with open(gfa_path, "r") as fin:
        for line in fin.readlines():
            if line.startswith("S"):
                lst = line.strip().split("\t")[1:]
                node_id, seq = lst[0], lst[1]
                new_record = SeqRecord(Seq(seq), id=node_id, description="", name="")
                records.append(new_record)
    fasta_name = gfa_path.parent / fasta
    SeqIO.write(records, fasta_name, "fasta")


if __name__ == "__main__":
    # test star conformation
    sequences_file = r"/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_hmm/probe_3/H1_C8_nodes.fasta"
    hmm_file = r"/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_hmm/probe_3/probe_3.hmm"
    pkl_hmm_file = "/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_hmm/probe_3/hmm_H1_C8.pkl"

    # test four separated contigs
    sequences_file = r"/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_hmm/probe_3/H1_B6_nodes.fasta"
    hmm_file = r"/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_hmm/probe_3/probe_3.hmm"
    pkl_hmm_file = "/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_hmm/probe_3/hmm_H1_B6.pkl"

    # test a really complex graph
    sequences_file = r"/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_hmm/probe_3/H1_A8_nodes.fasta"
    hmm_file = r"/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_hmm/probe_3/probe_3.hmm"
    pkl_hmm_file = "/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_hmm/probe_3/hmm_H1_A8.pkl"

    hmm = Hmmer_result(Path(sequences_file), Path(hmm_file))
    for k, v in hmm.best_paths.items():
        print(k, v)

    for k, v in hmm.node_hits.items():
        print(k, v)
    import pickle

    pickle.dump(hmm, open(pkl_hmm_file, "ab"))

    # test convert gfa to node sequences
    # gfa = Path(
    #     "/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_hmm/probe_3/H1_A8_blunted_compacted.gfa"
    # )
    # fasta = "H1_A8_nodes.fasta"
    #
    # extract_seq_from_gfa(gfa, fasta)
