"""Non-executable target API skeleton for Codex.

This file is a design contract, not production code. Copy signatures into the package
and implement them in the modules specified by CODEX_IMPLEMENTATION_SPEC.md.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Literal
import numpy as np

Array = np.ndarray
SLOT_NAMES = tuple([f"a{i}" for i in range(1,10)] + [f"b{i}" for i in range(1,6)] +
                   [f"c{i}" for i in range(1,7)] + [f"d{i}" for i in range(1,5)] +
                   [f"e{i}" for i in range(1,5)] + [f"v{i}" for i in range(1,14)])
DEFAULT_L = (6,12,20,42,110,420,1000)

@dataclass(frozen=True)
class P5Background:
    r: Array; u: Array; phi: Array; f: Array; h: Array; phi_r: Array; X: Array; A0prime: Array
    region: Array | None = None

@dataclass(frozen=True)
class ActionMember:
    release: str; epsilon_y: float; vector_branch: Literal["A0prime=0"]; source_definition: Path

@dataclass(frozen=True)
class Coefficients41:
    background: P5Background; slots: dict[str, Array]; provenance_id: str

@dataclass(frozen=True)
class ConstraintPivots:
    Dh1: Array; DeltaV: Array; pivotA0: Array

@dataclass(frozen=True)
class ReducedOperator:
    L: int; r: Array; K: Array; R: Array; G: Array; S: Array; M: Array
    pivots: ConstraintPivots; provenance_id: str

@dataclass(frozen=True)
class GateResult:
    name: str; status: str; value: object; criterion: str; evidence: str

# geometry
def load_frozen_p5(source: Path) -> P5Background: raise NotImplementedError
def validate_background(bg: P5Background) -> list[GateResult]: raise NotImplementedError
def kappa(bg: P5Background) -> Array: return bg.h * bg.phi_r**2

def locate_light_rings(bg: P5Background) -> dict[str,float]: raise NotImplementedError

# jets
def derivative(x: Array, y: Array, order: int=1, *, window: int=9, degree: int=8) -> Array: raise NotImplementedError

# action/coefficient production
def load_action_member(path: Path) -> ActionMember: raise NotImplementedError
def emit_mh_patch(action_patch, bg: P5Background) -> Coefficients41: raise NotImplementedError
def validate_41_schema(c: Coefficients41) -> list[GateResult]: raise NotImplementedError

def apply_epsilon_y_deformation(c: Coefficients41, epsilon_y: float) -> Coefficients41:
    """Frozen production law: only v1/v10 change on A0prime=0."""
    raise NotImplementedError

# constraints/reduction
def constraint_pivots(c: Coefficients41, L: int) -> ConstraintPivots: raise NotImplementedError
def reduce_profile(c: Coefficients41, L: int) -> ReducedOperator: raise NotImplementedError
def validate_operator(op: ReducedOperator) -> list[GateResult]: raise NotImplementedError

# global export
def build_global_coefficients(config) -> Coefficients41: raise NotImplementedError
def build_global_krgm(config, L_values: tuple[int,...]=DEFAULT_L) -> dict[int,ReducedOperator]: raise NotImplementedError
def write_direct_certificate(*args, **kwargs) -> Path: raise NotImplementedError

# qnm guard
def require_direct_krgm_certificate(path: Path): raise NotImplementedError
