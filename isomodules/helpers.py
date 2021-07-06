from portion.interval import Interval, Atomic, closed
from portion.const import Bound
from itertools import accumulate
from typing import Collection, Tuple
import matplotlib.colors as mc
import colorsys
import numpy as np


class IntegerInterval(Interval):
    """Represents an interval over the integers."""
    def __init__(self, *intervals):
        super().__init__(*intervals)
        # merge any atomic intervals that are separated by a gap of length 1
        # ex: [1, 3] | [4, 5] is the same as [1, 5]
        i = 0
        while i+1 < len(self._intervals):
            current = self._intervals[i]
            successor = self._intervals[i+1]
            if current.upper + 1 == successor.lower and Bound.CLOSED in (current.right, successor.left):
                union = Atomic(current.left, current.lower, successor.upper, successor.right)
                self._intervals.pop(i)  # pop current
                self._intervals.pop(i)  # pop successor
                self._intervals.insert(i, union)
            else:
                i += 1
    
    def to_tuples(self):
        return tuple((a.lower, a.upper) for a in self)

    @staticmethod
    def from_tuples(*intervals: Tuple[int, int]):
        return IntegerInterval(*[closed(a, b) for a, b in intervals])


def get_start_end_from_lengths(lengths: Collection[int]) -> Tuple[Tuple[int, int]]:
    """Convert a series of interval lengths to a series of interval (start, end) coordinates.
    Coordinates are 1-indexed and inclusive."""
    starts = accumulate(lengths[:-1], initial=1)
    ends = accumulate(lengths)
    return tuple(zip(starts, ends))


# originally from https://gist.github.com/ihincks/6a420b599f43fcd7dbd79d56798c4e5a
def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.
    
    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = np.array(colorsys.rgb_to_hls(*mc.to_rgb(c)))
    return colorsys.hls_to_rgb(c[0],1-amount * (1-c[1]),c[2])