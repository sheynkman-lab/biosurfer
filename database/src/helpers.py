# utility functions that don't fit in other modules
from bisect import bisect
from collections.abc import Mapping
from operator import itemgetter
from typing import Generic, Iterable, Tuple, TypeVar, Iterator
T = TypeVar('T')

class BisectDict(Mapping, Generic[T]):
    def __init__(self, items: Iterable[Tuple[int, T]]):
        self.breakpoints, self._values = zip(*sorted(items, key=itemgetter(0)))
    
    def __getitem__(self, key: int) -> T:
        if key < 0:
            raise KeyError('Key must be non-negative')
        i = bisect(self.breakpoints, key)
        try:
            return self._values[i]
        except IndexError as e:
            raise KeyError(key) from e
    
    def __iter__(self) -> Iterator[int]:
        yield from (0,) + self.breakpoints[:-1]
    
    def __len__(self):
        return len(self.breakpoints)