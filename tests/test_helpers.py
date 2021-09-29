from itertools import combinations, product
from typing import Dict, List, Tuple

from biosurfer.core.helpers import BisectDict
from biosurfer.plots.plotting import generate_subtracks
from hypothesis import given
from hypothesis.strategies import characters, data, integers, lists, dictionaries, sampled_from
from more_itertools import chunked

def force_list_length_to_be_even(lst: List) -> List:
    N = len(lst)
    return lst[:-1] if N % 2 else lst
    
non_overlapping_intervals = (
    lists(integers(min_value=0), unique=True, min_size=2).
    map(force_list_length_to_be_even).
    map(sorted).
    map(lambda lst: chunked(lst, 2, strict=True)).
    map(list)
)

label_intervals_mappings = dictionaries(
    characters(whitelist_categories=['Lu']),
    non_overlapping_intervals,
    min_size = 1
)

def extract_intervals_labels(label_to_intervals):
    labeled_intervals = ((label, interval) for label, intervals in label_to_intervals.items() for interval in intervals)
    labels, intervals = zip(*labeled_intervals)
    return intervals, labels, label_to_intervals

def intervals_overlap(intervals1, intervals2):
    return any(interval1[0] < interval2[1] and interval2[0] < interval1[1]
        for interval1, interval2 in product(intervals1, intervals2))

### TESTS BEGIN HERE ###

@given(
    dictionary = dictionaries(
        integers(min_value=1),
        characters(whitelist_categories=['Lu']),
        min_size = 1
    ),
    data = data()
)
def test_bisect_dict(dictionary, data):
    bdict = BisectDict(dictionary.items())
    breakpoints = sorted(dictionary.keys())
    i = data.draw(integers(min_value=0, max_value=len(breakpoints)-1))
    stop = breakpoints[i]
    start = breakpoints[i-1] if i > 0 else 0
    key = data.draw(integers(min_value=start, max_value=stop-1))
    assert bdict[key] == dictionary[stop]

@given(label_intervals_mappings.map(extract_intervals_labels))
def test_generate_subtracks(arg):
    intervals, labels, label_to_intervals = arg
    label_to_subtrack, _ = generate_subtracks(intervals, labels)
    subtrack_to_labels = dict()
    for label, subtrack in label_to_subtrack.items():
        subtrack_to_labels.setdefault(subtrack, []).append(label)
    for labels in subtrack_to_labels.values():
        for label1, label2 in combinations(labels, 2):
            assert not intervals_overlap(
                label_to_intervals[label1],
                label_to_intervals[label2]
                ), label_to_subtrack
