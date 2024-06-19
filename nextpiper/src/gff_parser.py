#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2024
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.

__version__ = "0.1"

import sys

# =======================================================================================
#               IMPORTS
# =======================================================================================
from dataclasses import dataclass, field, fields
from pathlib import Path
from typing import Final, Optional, Self, TypedDict, Literal, Any
import os

AAs = (
    "A",
    "R",
    "N",
    "D",
    "C",
    "E",
    "Q",
    "G",
    "H",
    "I",
    "L",
    "K",
    "M",
    "F",
    "P",
    "S",
    "T",
    "W",
    "Y",
    "V",
    "X",
)


@dataclass
class Fragment:
    """
    Data structure for the fragment object produced by the parsing of the miniprot output.
    Each fragment is one exon.
    Attributes
    ----------
    -contig_name: the name of the target of miniprot.
    -probe_name: the name of the query of miniprot.
    -target_start: in nucleotide space (start index at 0), start precedes end.
    -target_end: in nucleotide space (start index at 0).
    -query_start: in protein space (start index at 1), start precedes end.
    -query_end: in protein space (start index at 1).
    -score: alignment score.
    -strand of the contig.
    -frame of the contig.
    -identity: aa identity from miniprot.
    Post Init
    -correspondence: dictionary with keys as position on the query probe protein and values
        nucleotide positions on the contig. Difference in value between two consecutive keys should be at least 3.
        The index starts index at 0.

    """

    contig_name: str
    probe_name: str
    target_start: int
    target_end: int
    query_start: int
    query_end: int  # in protein space
    score: int
    strand: Literal[-1, 1]
    frame: Literal[0, 1, 2]
    identity: float
    correspondence: str = field(default_factory=dict, init=False)

    def get_correspondence(self) -> dict[int, int]:
        """
        :return:
        """
        return self.correspondence

    def get_contig_name(self) -> str:
        return self.contig_name

    def get_strand(self) -> Literal[-1, 1]:
        return self.strand

    def __repr__(self):
        try:
            description = (
                f"Fragment({self.contig_name}, {self.probe_name}, {self.target_start=}"
            )
            f" {self.target_end=}, {self.query_start=}, {self.query_end=},"
            f" {self.score=}, {self.strand=}, {self.frame=}, {self.identity=},"
            f" [{list(self.correspondence.items())[0]}-{list(self.correspondence.items())[-1]}]"
            return description
        except IndexError:
            return (
                f"Fragment({self.contig_name}, {self.probe_name}, {self.target_start=}"
                f" {self.target_end=}, {self.query_start=}, {self.query_end=},"
                f" {self.score=}, {self.strand=}, {self.frame=}, {self.identity=},"
            )


@dataclass
class Cds:
    """
    Data structure in charge of the parsing of the miniprot output.
    Attributes
    ----------
    -miniprot_output: lines of the miniprot output.
    Post Init
    -data: intermediate object that holds the output of the parsing.
    -fragments: list of Fragment objects. In theory miniprot should produce a single Fragment object.
    -mRNA_start: index of the start of the mRNA sequence (nucleotide space), start precedes end.
    -mRNA_end: index of the end of the mRNA sequence (nucleotide space).
    -probe_start: index of the start of the query (protein space), start precedes end.
    -probe_end: index of the end of the query (protein space).
    -target_nucleotides: sequence of the matching target (i.e. contig) in nucleotides.
    -target_AAs: sequence of the matching target (i.e. contig) in AAs.
    -query_AAs: sequence of the matching query (i.e. probe) in AAs.
    -global_identity: aa identity from miniprot over all exons (fragments).
    """

    miniprot_output: str = field(repr=False)
    data: int = field(default_factory=dict, init=False, repr=False)
    fragments: list[Fragment] = field(default_factory=list, init=False, repr=False)
    mRNA_start: int = field(init=False)
    mRNA_end: int = field(init=False)
    probe_start: int = field(init=False)
    probe_end: int = field(init=False)
    target_nucleotides: str = field(
        init=False, repr=False
    )  # target sequence in nucleotides
    target_AAs: str = field(init=False, repr=False)  # target sequence in aa
    query_AAs: str = field(init=False, repr=False)  # query sequence in aa
    global_identity: float = field(init=False)

    def __repr__(self):
        if self.fragments:
            return f"Cds({self.mRNA_start=},{self.mRNA_end=},{self.probe_start=},{self.probe_end=},{self.global_identity=})"
        else:
            return f"Cds()"

    def __post_init__(self):
        self._parse_gff()
        # If the data is not empty
        if list(self.data.keys()):
            self._find_correspondences()
            # Ensure that all variables have been set
            for field in fields(self):
                if field.name not in vars(self):
                    raise AttributeError(f"Field {field.name} has no attribute")
            assert (
                len(self.target_nucleotides)
                == len(self.target_AAs)
                == len(self.query_AAs)
            ), f"length of target_nucleotides, target_AAs and query_AAs do not match"

    def _parse_gff(self):
        for line in self.miniprot_output:
            if (
                line.startswith("##gff")
                or line.startswith("##PAF")
                or line.startswith("##AAS")
            ):
                continue
            # target (contig) nucleotides
            if line.startswith("##ATN"):
                self.target_nucleotides = line.removeprefix("##ATN\t").replace("\n", "")
                continue
            # target (contig) amino acids
            if line.startswith("##ATA"):
                self.target_AAs = line.removeprefix("##ATA\t").replace("\n", "")
                continue
            # query (probe) amino acids
            if line.startswith("##AQA"):
                self.query_AAs = line.removeprefix("##AQA\t").replace("\n", "")
                continue
            # mRNA & CDS lines
            record = line.strip().split("\t")
            sequence_name = record[0]
            source = record[1]
            feature = record[2]
            start = int(record[3])
            end = int(record[4])
            if record[5] != ".":
                score = record[5]
            else:
                score = None
            if record[6] == "+":
                strand = 1
            else:
                strand = -1
            if record[7] != ".":
                frame = record[7]
            else:
                frame = None
            attributes = record[8].split(";")
            attributes = {
                x.split("=")[0]: x.split("=")[1] for x in attributes if "=" in x
            }
            if not (sequence_name in self.data):
                self.data[sequence_name] = []
            alpha = {
                "source": source,
                "feature": feature,
                "start": int(start),
                "end": int(end),
                "score": int(score),
                "strand": strand,
                "frame": int(frame) if frame is not None else None,
            }
            for k, v in attributes.items():
                alpha[k] = v
            self.data[sequence_name].append(alpha)
        contigs = list(self.data.keys())
        # If the file is empty stop here
        if not contigs:
            return
        # Miniprot should map to a single location
        assert len(contigs) == 1, "more than one contig is present"
        contig = contigs[0]
        for item in self.data[contig]:
            if not item.get("Target"):
                continue

            target_splt = item["Target"].split()
            probe_name = target_splt[0]
            query_start = int(target_splt[1])
            query_end = int(target_splt[2])
            if item["feature"] == "CDS":
                fragment = Fragment(
                    contig,
                    probe_name,
                    item["start"] - 1,
                    item["end"] - 1,
                    # query_start,
                    query_start,
                    # query_end - 1,
                    query_end,
                    item["score"],
                    item["strand"],
                    item["frame"],
                    float(item["Identity"]),
                )
                self.fragments.append(fragment)
            else:
                self.mRNA_start = item["start"] - 1
                self.mRNA_end = item["end"] - 1
                self.probe_start = int(target_splt[1]) - 1
                self.probe_end = int(target_splt[2]) - 1
                self.global_identity = float(item["Identity"])

    def _find_correspondences(self):
        frg_idx = 0
        fragment = self.fragments[frg_idx]
        frg_limit = fragment.query_end
        idx = 0
        nuc_idx = self.mRNA_start if fragment.strand == 1 else self.mRNA_end
        query_AA_idx = self.probe_start
        probe_nucl_cor = {}
        while idx < len(self.target_nucleotides) + 1:
            try:
                target_nuc = self.target_nucleotides[idx]
            except:
                fragment.correspondence = probe_nucl_cor
                break
            target_AA = self.target_AAs[idx]
            query_AA = self.query_AAs[idx]
            # if query_AA in AAs:
            #     last_AA = query_AA
            if target_nuc == "-":
                idx += 1
                continue
            if target_AA in AAs and query_AA in AAs:
                probe_nucl_cor[query_AA_idx] = nuc_idx
            if query_AA in AAs:
                query_AA_idx += 1
            idx += 1
            nuc_idx += 1 * fragment.strand
            if query_AA_idx > frg_limit:
                fragment.correspondence = probe_nucl_cor
                probe_nucl_cor = {}
                frg_idx += 1
                if frg_idx == len(self.fragments):
                    break
                fragment = self.fragments[frg_idx]
                frg_limit = fragment.query_end

    def get_fragments(self) -> Optional[list[Fragment]]:
        if not self.data.keys():
            return None
        return self.fragments

    def get_global_sim(self):
        return self.global_identity

    def is_empty(self):
        return not bool(self.data)


if __name__ == "__main__":
    # os.chdir(
    #     "/home/yjkbertrand/Documents/projects/nextpiper/test_data/batrachium/exonerate"
    # )
    os.chdir(
        "/home/yjkbertrand/Documents/projects/nextpiper/test_data/test_clustering/gff_files"
    )

    # cds_list = CdsParser(open("test_8631_node3_rc.gff", "r"))
    # cds_list = CdsParser(open("test_8631_node3.gff", "r"))
    # cds_list = CdsParser(
    #     open("test_8631_probe_start_9_intron_1_micorinton1_strand_1.gff", "r")
    # )
    cds_list = Cds(open("gene_3_A3_1.gff", "r"))

    print(cds_list)
    fragments = cds_list.get_fragments()
    print(len(fragments))
    for frg in fragments:
        print(frg)
    for frg in fragments:
        print(list(frg.correspondence.items()))
