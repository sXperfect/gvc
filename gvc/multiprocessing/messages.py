"""Picklable message types for the GVC multiprocessing pipeline."""

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class WorkItem:
    block_id: int
    raw_block: object


@dataclass(frozen=True)
class EncodedBlock:
    block_id: int
    parameter_set: object
    block: object


@dataclass(frozen=True)
class StopWork:
    pass


@dataclass(frozen=True)
class WorkerDone:
    worker_id: int


@dataclass(frozen=True)
class ReaderDone:
    total_blocks: int


@dataclass(frozen=True)
class WriterDone:
    total_blocks: int


@dataclass(frozen=True)
class WorkerError:
    stage: str
    worker_id: Optional[int]
    block_id: Optional[int]
    error_type: str
    message: str
    traceback: str
