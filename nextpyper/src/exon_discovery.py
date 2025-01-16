#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2025
#    Simón Villanueva CORRALES: simon.corrales@ibot.cas.cz
#    Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#
#       All rights reserved.
"""
Function and classes for defining the best
"""
import sys
from collections import deque
from dataclasses import dataclass, field, fields
from itertools import chain, product
from typing import Optional, Self, NamedTuple
from intervaltree import Interval, IntervalTree

from exon_intron import ItervalGraph, Exon

# =======================================================================================
#               FUNCTIONS
# =======================================================================================


def fuse_intervals(intervals: list[Interval]) -> Interval:
    """Fuse several Interval objects into a single Interval.
    data attributes are fused into a single list, lower and higher bounds
    are set to min and max values respectively."""
    min_value = sorted(intervals, key=lambda i: i.begin)[0].begin
    max_value = sorted(intervals, key=lambda i: i.end)[-1].end
    data = list(chain(*[x.data for x in intervals]))
    return Interval(min_value, max_value, data)


def correctly_ordered(
    interval_0: tuple[Interval, str], interval_1: tuple[Interval, str]
) -> bool:
    """
    Interval_0 and interval_1 are tuples of interval and 'start'/'stop'.
    Check that two consecutive intervals are correctly ordered, which mean
    checking that if they contain the same name, the second idx cannot larger than the first one.
    """
    # two consecutive starts or stops
    if interval_0[1] == interval_1[1]:
        return True
    sequence_0_dict = {endpt.name: endpt.idx for endpt in interval_0[0].data}
    sequence_1_dict = {endpt.name: endpt.idx for endpt in interval_1[0].data}
    common_names = set(sequence_0_dict) & set(sequence_1_dict)
    if not common_names:
        return True
    if interval_0[1] == "stop":
        return False
    return True


def process_it(intervals: list[Interval]) -> IntervalTree:
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


def find_longest_exon_stretch(
    putative_exons: list["PutativeExon"], expansion_interval: int = 0
) -> list["Exon"]:
    """
    Arrange the putative exons to find the longest exon stretch.
    The search if done with and ItervalGraph structure.
    -expansion_interval: number of AAs to substract from the start of the exon, to decide if it overlaps with the end.
    """
    possible_exons = []
    for putative_exon in putative_exons:
        start_points = set([endpt.idx for endpt in putative_exon.start_sequences])
        end_points = set([endpt.idx for endpt in putative_exon.end_sequences])
        valid_combinations = [
            x
            for x in product(start_points, end_points)
            if x[0] - expansion_interval < x[1]
        ]
        for combination in valid_combinations:
            possible_exons.append(Exon(*combination))
    IG = ItervalGraph(possible_exons)
    return IG.get_best_path().path


@dataclass
class PutativeExon:
    """
    Keep the information about the actual sequence starting and ending points.
    Unlike the Exon object, here the start and end are a list of all empirical coordinates that are deemed
    to correspond to the same biological exon.
    """

    __slots__ = ["start_sequences", "end_sequences"]
    start_sequences: list[str]
    end_sequences: list[str]


@dataclass
class EndPoint:
    """
    Container for start or end of an exon.
    """

    __slots__ = ["idx", "name"]
    idx: int
    name: str


@dataclass
class DiscoverExons:
    """
    Data structure for selecting the most common exon boundaries
    Attributes
    ----------
    -exon_intervals: list of exon intervals found by miniprot boundary scorer for the scaffolds.
    -expansion_threshold: length proportion of the exon that are expanded around the exons endpoints in oder
        to form an interval. Start/stop intervals are then compared between putative exons.
    Post Init
    -starts: for each start point, the expansion produces an interval.
    -ends: for each end point, the expansion produces an interval.
    -exons: list of Exon objects.
    """

    exon_intervals: list[Interval]
    expansion_threshold: float = field(default=0.1)
    starts: list[Interval] = field(default_factory=list, init=False, repr=True)
    ends: list[Interval] = field(default_factory=list, init=False, repr=True)
    exons: list[Exon] = field(default_factory=list, init=False, repr=True)

    def __post_init__(self):
        self._get_starts_ends()
        self._find_exons()

    def get_exons(self) -> list[Exon]:
        return self.exons

    def _get_starts_ends(self) -> Self:
        """
        Populate the starts and ends attributes.
        MAX_EXPANSION_INTERVAL is the range in AA over an index that can encompass a second index.
        """
        MAX_EXPANSION_INTERVAL = 10
        for elt in self.exon_intervals:
            # Remove the shortest exons
            if elt.end - elt.begin <= MAX_EXPANSION_INTERVAL:
                continue
            expansion = max(
                MAX_EXPANSION_INTERVAL,
                int((elt.end - elt.begin) * self.expansion_threshold),
            )
            start_min, start_max = max(elt.begin - expansion, 0), elt.begin + expansion
            end_min, end_max = elt.end - expansion, elt.end + expansion
            self.starts.append(
                Interval(start_min, start_max, [EndPoint(elt.begin, elt.data)])
            )
            self.ends.append(Interval(end_min, end_max, [EndPoint(elt.end, elt.data)]))
        return self

    def _find_exons(self) -> Self:
        """
        Populate the exon attribute.
        """
        # Fuse overlapping intervals for start and end points of exon intervals using interval tree.
        intervals_start = sorted(process_it(self.starts), key=lambda i: i.begin)
        intervals_end = sorted(process_it(self.ends), key=lambda i: i.end)
        # Process the sorted list of start and end intervals.
        all_intervals = []
        all_intervals.extend((interval, "start") for interval in intervals_start)
        all_intervals.extend((interval, "stop") for interval in intervals_end)
        tmp_sorted_intervals = deque(sorted(all_intervals, key=lambda x: x[0].begin))
        first_interval = tmp_sorted_intervals.popleft()
        sorted_intervals = [first_interval]
        # Reorder the intervals
        while tmp_sorted_intervals:
            second_interval = tmp_sorted_intervals.popleft()
            if not correctly_ordered(first_interval, second_interval):
                sorted_intervals.pop()
                sorted_intervals.extend([second_interval, first_interval])
            else:
                sorted_intervals.append(second_interval)
                first_interval = second_interval

        print("sorted intervals", sorted_intervals)
        putative_exons = []
        start = sorted_intervals[0]
        start_interval = start[0]
        current_color = start[1]
        stop_interval = None
        queue = deque(sorted_intervals[1:])
        while queue:
            print(f"{start_interval=}, {stop_interval=}, {current_color=}")
            elt = queue.popleft()
            # print("length of queue:", len(queue))
            elt_interval = elt[0]
            elt_color = elt[1]
            # print(f"new elm {elt}")
            if current_color == "start":
                if elt_color == "start":
                    if len(elt_interval.data) > len(start_interval.data):
                        start_interval = elt_interval
                    if len(elt) == 0:
                        # Something went wrong and the last colour is not "stop"
                        break
                elif elt_color == "stop":
                    current_color = "stop"
                    stop_interval = elt_interval
                    # End of list
                    if len(elt) == 0:
                        # print("appending")
                        putative_exons.append(
                            PutativeExon(start_interval.data, stop_interval.data)
                        )

            elif current_color == "stop":
                if elt_color == "start":
                    # print("appending")
                    putative_exons.append(
                        PutativeExon(start_interval.data, stop_interval.data)
                    )
                    start_interval = elt_interval
                    current_color = "start"
                    stop_interval = None
                elif elt_color == "stop":
                    # switch to this new stop if there are more sequences that support this stop
                    # or if there are the same number but contain at least one common sequence withthe start.

                    if len(elt_interval.data) > len(stop_interval.data):
                        stop_interval = elt_interval
                    elif len(elt_interval.data) == len(stop_interval.data):
                        if {d.name for d in elt_interval.data} & {
                            d.name for d in stop_interval.data
                        }:
                            stop_interval = elt_interval
            if len(queue) == 0:
                # print("appending")
                putative_exons.append(
                    PutativeExon(start_interval.data, stop_interval.data)
                )
        longest_stretch_exons = find_longest_exon_stretch(putative_exons)
        if not longest_stretch_exons:
            return self
        # The boundaries need to be adjusted
        queue = deque(longest_stretch_exons)
        query = queue.popleft()
        self.exons.append(query)
        while queue:
            target = queue.popleft()
            if query.end != target.start:
                self.exons.append(target)
            else:
                self.exons.append(Exon(target.start + 1, target.end))
            query = target
        return self


if __name__ == "__main__":

    inputs = [
        (40, 70, "x"),
        (40, 70, "a"),
        (50, 70, "b"),
        (70, 140, "c"),
        (71, 110, "a"),
        (107, 142, "b"),
        (112, 142, "a"),
        (143, 178, "a"),
        (143, 170, "b"),
        (112, 191, "d"),
    ]

    inputs = [
        (86, 180, "a"),
        (89, 142, "b"),
        (91, 191, "c"),
        (93, 191, "d"),
    ]

    input_intervals = [Interval(*x) for x in inputs]
    disco = DiscoverExons(input_intervals)
    print("putative exons", disco.exons)
    # inputs_start = [
    #     (38, 42, [("x", 40)]),
    #     (38, 42, [("a", 40)]),
    #     (48, 55, [("b", 51)]),
    #     (58, 62, [("c", 60)]),
    #     (69, 73, [("a", 71)]),
    #     (105, 110, [("b", 107)]),
    #     (109, 115, [("c", 111)]),
    #     (140, 150, [("a", 143)]),
    #     (148, 152, [("b", 150)]),
    #     (188, 192, [("d", 190)]),
    #     (40, 44, [("X", 42)]),
    # ]
    # inputs_end = [
    #     (68, 72, [("x", 70)]),
    #     (68, 72, [("a", 70)]),
    #     (68, 72, [("b", 70)]),
    #     (78, 82, [("c", 80)]),
    #     (108, 112, [("a", 110)]),
    #     (138, 144, [("b", 142)]),
    #     (138, 144, [("c", 142)]),
    #     (176, 180, [("a", 178)]),
    #     (168, 172, [("b", 170)]),
    #     (190, 200, [("d", 195)]),
    #     (188, 192, [("X", 190)]),
    # ]
    # t_start = IntervalTree.from_tuples(inputs_start)
    # used_centers = []
    # for item in inputs_start:
    #     center = item[2][0][1]
    #     if center in used_centers:
    #         continue
    #     centered_intervals = sorted(t_start[center])
    #     new_interval = fuse_intervals(centered_intervals)
    #     del t_start[center]
    #     t_start.add(new_interval)
    #     used_centers.extend(x[1] for x in new_interval.data)
    # t_end = IntervalTree.from_tuples(inputs_end)
    # used_centers = []
    # for item in inputs_end:
    #     center = item[2][0][1]
    #     if center in used_centers:
    #         continue
    #     centered_intervals = sorted(t_end[center])
    #     new_interval = fuse_intervals(centered_intervals)
    #     del t_end[center]
    #     t_end.add(new_interval)
    #     used_centers.extend(x[1] for x in new_interval.data)
    # print(t_start)
    # print(t_end)
    # intervals_start = sorted(t_start, key=lambda i: i.begin)
    # intervals_end = sorted(t_end, key=lambda i: i.end)
    # all_intervals = []
    # all_intervals.extend((interval, "start") for interval in intervals_start)
    # all_intervals.extend((interval, "stop") for interval in intervals_end)
    # sorted_intervals = sorted(all_intervals, key=lambda x: x[0].begin)
    # print("sorted intervals", sorted_intervals)
    # exon_intervals = []
    # start = sorted_intervals[0]
    # start_interval = start[0]
    # current_color = start[1]
    # stop_interval = None
    # queue = deque(sorted_intervals[1:])
    # while queue:
    #     print(f"{start_interval=}, {stop_interval=}, {current_color=}")
    #     elt = queue.popleft()
    #     elt_interval = elt[0]
    #     elt_color = elt[1]
    #     print(f"new elm {elt}")
    #     if current_color == "start":
    #         if elt_color == "start":
    #             if len(elt_interval.data) > len(start_interval.data):
    #                 start_interval = elt_interval
    #         elif elt_color == "stop":
    #             current_color = "stop"
    #             stop_interval = elt_interval
    #             # End of list
    #             if len(elt) == 0:
    #                 start = min([x[1] for x in start_interval.data])
    #                 stop = max([x[1] for x in stop_interval.data])
    #                 new_data = [x[0] for x in stop_interval.data]
    #                 print("appending")
    #                 exon_intervals.append(PutativeExon(start, stop, new_data))
    #
    #     elif current_color == "stop":
    #         if elt_color == "start":
    #             start = min([x[1] for x in start_interval.data])
    #             stop = max([x[1] for x in stop_interval.data])
    #             new_data = [x[0] for x in stop_interval.data]
    #             print("appending")
    #             exon_intervals.append(PutativeExon(start, stop, new_data))
    #             start_interval = elt_interval
    #             current_color = "start"
    #             stop_interval = None
    #         elif elt_color == "stop":
    #             # switch to this new stop if there are more sequences that support this stop
    #             # or if there are the same number but contain at least one common sequence withthe start.
    #             if len(elt_interval.data) > len(stop_interval.data):
    #                 end_interval = elt_interval
    #             elif len(elt_interval.data) == len(stop_interval.data):
    #                 if {d[0] for d in elt_interval.data} & {
    #                     d[0] for d in stop_interval.data
    #                 }:
    #                     end_interval = elt_interval
    # print("final", exon_intervals)
