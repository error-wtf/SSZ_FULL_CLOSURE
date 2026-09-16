from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np

Array = np.ndarray


@dataclass(frozen=True)
class P5Background:
    r: Array
    u: Array
    phi: Array
    f: Array
    h: Array
    phi_r: Array
    X: Array
    A0prime: Array
    region: Array | None = None


@dataclass(frozen=True)
class ActionMember:
    release: str
    epsilon_y: float
    vector_branch: Literal["A0prime=0"]
    source_definition: Path


@dataclass(frozen=True)
class Coefficients41:
    background: P5Background
    slots: dict[str, Array]
    provenance_id: str


@dataclass(frozen=True)
class ConstraintPivots:
    Dh1: Array
    DeltaV: Array
    pivotA0: Array


@dataclass(frozen=True)
class ReducedOperator:
    L: int
    r: Array
    K: Array
    R: Array
    G: Array
    S: Array
    M: Array
    pivots: ConstraintPivots
    provenance_id: str


@dataclass(frozen=True)
class GateResult:
    name: str
    status: Literal["PASS", "FAIL", "WARN", "OPEN", "PASS_EXACT", "PASS_NUMERICAL"]
    value: object
    criterion: str
    evidence: str


@dataclass(frozen=True)
class DirectKRGMCertificate:
    pass_: bool
    release: str
    action_sha256: str
    coefficient_stream_sha256: str
    L_values: tuple[int, ...]
    gates: tuple[GateResult, ...]
