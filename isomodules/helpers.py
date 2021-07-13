from portion.interval import Interval, Atomic, inf, closed
from portion.const import Bound
from itertools import accumulate
from typing import Collection, Tuple


class IntegerInterval(Interval):
    """Represents an interval over the integers."""

    # convert all atomic intervals to open intervals, merge, then convert back to closed
    # from https://github.com/AlexandreDecan/portion/issues/24#issuecomment-604456362
    @staticmethod
    def expand(s):
        lower, upper = s.lower, s.upper
        if s.left is Bound.CLOSED:
            lower -= 1
        if s.right is Bound.CLOSED:
            upper += 1
        return Atomic(Bound.OPEN, lower, upper, Bound.OPEN)

    @staticmethod
    def reduce(s):
        if s.lower == inf or s.upper == -inf:
            return Atomic(Bound.OPEN, inf, -inf, Bound.OPEN)
        lower, upper = s.lower, s.upper
        if s.left is Bound.OPEN:
            lower += 1  
        if s.right is Bound.OPEN:
            upper -= 1
        return Atomic(Bound.CLOSED, lower, upper, Bound.CLOSED)

    def __init__(self, *intervals):
        super().__init__(*intervals)
        self._intervals = self.apply(IntegerInterval.expand).apply(IntegerInterval.reduce)._intervals
    
    def to_tuples(self):
        if self.empty:
            return ()
        return tuple((a.lower, a.upper) for a in self)

    @staticmethod
    def from_tuples(*intervals):
        return IntegerInterval(*[closed(a, b) for a, b in intervals])
    
    # Return number of integers contained in interval.
    @property
    def size(self):
        if self.empty:
            return 0
        return sum(s.upper - s.lower + 1 for s in self._intervals)


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