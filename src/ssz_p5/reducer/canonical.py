"""Public profile reduction API backed by the existing golden implementation."""

import numpy as np
import pandas as pd

from ..coefficients.identities import require_v12
from ..coefficients.schema41 import validate_41_schema
from ..config import SLOT_NAMES
from ..geometry.p5 import validate_background
from ..numerics import module
from ..policy import numerical_policy
from ..types import ConstraintPivots, ReducedOperator


def coefficient_frame(c):
    validate_background(c.background)
    if any(g.status != "PASS" for g in validate_41_schema(c)):
        raise ValueError("invalid 41-slot schema")
    require_v12(c.slots, c.background.h)
    b = c.background
    d = pd.DataFrame(
        dict(x=b.r, u=b.u, phi=b.phi, f=b.f, h=b.h, phiprime=b.phi_r, A0prime=b.A0prime, X=b.X)
    )
    for name in SLOT_NAMES:
        d[name] = c.slots[name]
    residual = c.slots["v7"] - c.slots["v2"] ** 2 / (4 * c.slots["v1"])
    if (
        not np.all(np.isfinite(residual))
        or np.max(np.abs(residual)) > numerical_policy()["branch_identity_abs"]
    ):
        raise ValueError("auxiliary identity violated")
    return d


def require_pivots(m):
    for name in ("Dh1", "DeltaV", "pivotA0"):
        x = np.asarray(m[name])
        if not np.all(np.isfinite(x)) or np.min(np.abs(x)) <= numerical_policy()["pivot_floor"]:
            raise ValueError(f"forbidden constraint pivot: {name}")


def validate_operator(op):
    policy = numerical_policy()
    for name in ("K", "R", "G", "S", "M"):
        a = np.asarray(getattr(op, name))
        if a.shape != (len(op.r), 3, 3) or not np.all(np.isfinite(a)):
            raise ValueError(f"invalid matrix: {name}")
    for name, sign in [("K", -1), ("G", -1), ("M", -1), ("S", 1)]:
        a = getattr(op, name)
        err = np.max(np.abs(a + sign * a.swapaxes(1, 2))) / max(1.0, np.max(np.abs(a)))
        if err > policy["matrix_symmetry_scaled"]:
            raise ValueError(f"matrix symmetry: {name}")
    if (
        np.max(np.abs(op.R)) / max(1.0, np.max(np.abs(op.K)), np.max(np.abs(op.G)))
        > policy["R_zero_scaled"]
    ):
        raise ValueError("nonzero R")
    require_pivots(vars(op.pivots))
    return op


def reduce_profile(c, L):
    if L < 6 or int(L) != L:
        raise ValueError("L must be an integer >=6")
    d = coefficient_frame(c)
    red = module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    maps = red.fm.maps(d, L)
    require_pivots(maps)
    result = red.canonical_audit(d, L)
    pivots = ConstraintPivots(maps["Dh1"], maps["DeltaV"], maps["pivotA0"])
    return validate_operator(
        ReducedOperator(
            L,
            c.background.r,
            *(result[k] for k in ("K", "R", "G", "S", "M")),
            pivots,
            c.provenance_id,
        )
    )
