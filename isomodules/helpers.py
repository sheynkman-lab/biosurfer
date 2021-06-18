from portion.interval import Interval, Atomic, closed
from portion.const import Bound
from typing import Sequence, Tuple

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