from biosurfer.isomodules.helpers import IntegerInterval as ii
from random import randint

def _helper(atoms, union):
    for a, b in atoms:
        for k in range(a, b+1):
            assert k in union, f"{atoms} did not match {union}!"

def test_interval_union():
    atoms = [(0, 3), (1, 2), (5, 7), (8, 9)]
    _helper(atoms, ii.from_tuples(*atoms))

def test_interval_union_random():
    for _ in range(20):
        atoms = [((a := randint(0, 9)), a + randint(1, 4)) for _ in range(5)]
        _helper(atoms, ii.from_tuples(*atoms))

