# utility functions that don't fit in other modules
import sys
import traceback
from bisect import bisect
from collections.abc import Mapping
from contextlib import AbstractContextManager
from copy import copy
from dataclasses import dataclass, field, fields
from enum import Enum
from itertools import chain, count
from operator import itemgetter
from typing import TYPE_CHECKING, Callable, Generic, Iterable, Iterator, List, Optional, Tuple, TypeVar

from graph_tool import Graph
from intervaltree import Interval, IntervalTree
from sqlalchemy.dialects.sqlite.dml import insert

if TYPE_CHECKING:
    from io import TextIOBase

T = TypeVar('T')


class OrderedEnum(Enum):
    # https://docs.python.org/3/library/enum.html#orderedenum
    def __ge__(self, other):
        if self.__class__ is other.__class__:
            return self.value >= other.value
        return NotImplemented
    def __gt__(self, other):
        if self.__class__ is other.__class__:
            return self.value > other.value
        return NotImplemented
    def __le__(self, other):
        if self.__class__ is other.__class__:
            return self.value <= other.value
        return NotImplemented
    def __lt__(self, other):
        if self.__class__ is other.__class__:
            return self.value < other.value
        return NotImplemented


class StringEnum(Enum):
    def __str__(self):
        return self.value


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


def frozendataclass(cls):
    frozencls = dataclass(cls, frozen=True)
    field_names = {field.name for field in fields(frozencls)}
    def replace(self, **kwargs):
        """Return new instance of frozendataclass with updated values."""
        new_field_values = {name: kwargs.get(name, getattr(self, name)) for name in field_names}
        return frozencls(**new_field_values)
    frozencls.replace = replace
    return frozencls


class ExceptionLogger(AbstractContextManager):
    def __init__(self, info=None):
        self.info = info
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type is not None:
            sys.stderr.write('---------\n')
            if self.info:
                sys.stderr.write(str(self.info) + '\n')
            traceback.print_exc()
            sys.stderr.write('---------\n')
            return True


def run_length_encode(text: str) -> str:
    if not text:
        return ''
    encoding = []
    run_length = 1
    prev_char = text[0]
    for char in text[1:]:
        if char == prev_char:
            run_length += 1
        else:
            encoding.append(f'{run_length}{prev_char}')
            prev_char = char
            run_length = 1
    encoding.append(f'{run_length}{prev_char}')
    return ','.join(encoding)


def run_length_decode(encoding: str) -> str:
    return ''.join(int(token[:-1]) * token[-1] for token in encoding.split(',')) if encoding else ''


def get_interval_overlap_graph(intervals: Iterable[Tuple[int, int]], labels: Optional[Iterable] = None, label_type: str = 'string') -> 'Graph':
    # inspired by https://stackoverflow.com/a/19088519
    # build graph of labels where labels are adjacent if their intervals overlap
    if not labels:
        labels = count()
    g = Graph(directed=False)
    g.vp.label = g.new_vertex_property(label_type)
    label_to_vertex = dict()
    active_labels = set()
    boundaries = sorted(
        chain.from_iterable(
            [(a, True, label), (b, False, label)]
            for (a, b), label in zip(intervals, labels)
        ),
        key = itemgetter(0, 1)
    )
    for _, start_of_interval, label in boundaries:
        if start_of_interval:
            if label not in label_to_vertex:
                v = g.add_vertex()
                g.vp.label[v] = label
                label_to_vertex[label] = v
            for other_label in active_labels:
                i = label_to_vertex[label]
                j = label_to_vertex[other_label]
                g.add_edge(i, j)
            active_labels.add(label)
        else:
            active_labels.discard(label)
    return g, label_to_vertex


# Helper functions/classes for loading into database
@dataclass
class FastaHeaderFields:
    transcript_id: str = None
    protein_id: str = None


def count_lines(file_handle: 'TextIOBase', only: Optional[Callable[..., bool]] = None):
    lines = sum(1 for _ in filter(only, file_handle))
    file_handle.seek(0)
    return lines


def bulk_upsert(session, table, records, primary_keys=('accession',)):
    if records:
        fields = [field for field in records[0] if field not in primary_keys]
        with session.begin():
            stmt = insert(table)
            session.execute(
                stmt.on_conflict_do_update(
                    index_elements = primary_keys,
                    set_ = {field: stmt.excluded[field] for field in fields}
                ),
                records
            )
        records[:] = []


def read_gtf_line(line: str) -> list:
    """Read and parse a single gtf line

    Args:
        line (str): unbroken line of a gtf file

    Returns:
        list: gtf attributes
        chromosome : str
        source : str
        feature : str
        start : int
        stop : int
        score : str
        strand : str
        phase : str
        attributes: dict
        tags: list

    """
    chromosome, source, feature, start, stop, score, strand, phase, attributes = line.split('\t')
    start = int(start)
    stop = int(stop)
    attributes = attributes.split(';')[:-1]
    attributes = [att.strip(' ').split(' ') for att in attributes]
    tags = [att[1].strip('"') for att in attributes if att[0] == 'tag']
    attributes = {att[0]: att[1].strip('"') for att in attributes if att[0] != 'tag'}
    return chromosome, source, feature, start, stop, score, strand, phase, attributes, tags


def get_ids_from_gencode_fasta(header: str):
    fields = [field for field in header.split('|') if field and not field.startswith(('UTR', 'CDS'))]
    transcript_id = next((field for field in fields if field.startswith('ENST')))
    protein_id = next((field for field in fields if field.startswith('ENSP')), None)
    return FastaHeaderFields(transcript_id, protein_id)


def skip_par_y(header: str):  # these have duplicate ENSEMBL accessions and that makes SQLAlchemy very sad
    return 'PAR_Y' in header
