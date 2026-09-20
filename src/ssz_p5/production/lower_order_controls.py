"""Background-null lower-order action-control map for the Inner handover.

The selected normal coordinates are
    q1 = f2_{phi F}, q2 = f2_{X phi}, q3 = f2_{phi phi}.
On the locked static P5 background their lower-order response is triangular.
For the representative used here the independent diagonal responses are

    dv5/dq1 = r^2 sqrt(h/f) A0'
    dc3/dq2 = 1/2 r^2 sqrt(f h) (phi')^2
    de3/dq3 = -1/2 r^2 sqrt(f/h)

Taken as formal jet coordinates, these slots do not enter the local principal
response map.  However, a *single smooth background-null deformation* must also
satisfy the Hessian chain rule H @ (phi',X',F') = 0.  Consequently these three
mixed/phi jets are not independent of f2XX/f2XF/f2FF.  This module certifies
only the local triangular response inverse; common-action holonomy is audited
separately in :mod:`ssz_p5.production.holonomic_hessian`.

At A0'=0 the v5 response vanishes.  The C-infinity Inner target also has
v5=0 there, so the endpoint is handled by the reduced (c3,e3) block and no
singular division is performed.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

CONTROL_NAMES = ("f2phiF_control", "f2Xphi_control", "f2phiphi_control")
TARGET_NAMES = ("v5", "c3", "e3")


def response_diagonal(background: pd.DataFrame) -> np.ndarray:
    """Return the three exact diagonal response factors on the radial grid."""
    r = background.x.to_numpy(float)
    f = background.f.to_numpy(float)
    h = background.h.to_numpy(float)
    A = background.A0prime.to_numpy(float)
    if "phiprime" in background:
        ph = background.phiprime.to_numpy(float)
    elif "X" in background:
        ph = -np.sqrt(np.maximum(0.0, -2.0 * background.X.to_numpy(float) / h))
    else:
        raise KeyError("phiprime or X required for lower-order control map")
    if not all(np.isfinite(v).all() for v in (r, f, h, A, ph)):
        raise ValueError("nonfinite lower-order control background")
    if np.any(f <= 0) or np.any(h <= 0) or np.any(np.abs(ph) < 1e-14):
        raise ValueError("lower-order control map requires f>0,h>0,phiprime!=0")
    return np.column_stack(
        [
            r**2 * np.sqrt(h / f) * A,
            0.5 * r**2 * np.sqrt(f * h) * ph**2,
            -0.5 * r**2 * np.sqrt(f / h),
        ]
    )


def invert_lower_order_targets(
    background: pd.DataFrame,
    targets: pd.DataFrame | np.ndarray,
    *,
    tolerance: float = 1e-10,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Invert the lower-order normal-coordinate map and replay its response.

    ``targets`` may be a frame with v5_target/c3_target/e3_target columns or an
    (N,3) array.  The returned action-control frame is an explicit radial
    section of the three normal jets.  ``emitted`` contains the lower-order
    coefficients obtained from that action-control section, not values copied
    independently into a 41-slot table.
    """
    if isinstance(targets, pd.DataFrame):
        target = targets[[f"{k}_target" for k in TARGET_NAMES]].to_numpy(float)
    else:
        target = np.asarray(targets, float)
    if target.shape != (len(background), 3) or not np.isfinite(target).all():
        raise ValueError("invalid lower-order targets")

    J = response_diagonal(background)
    q = np.zeros_like(target)
    for j, name in enumerate(TARGET_NAMES):
        factor = J[:, j]
        scale = np.maximum(1.0, np.abs(target[:, j]))
        zero = np.abs(factor) <= 1e-14
        impossible = zero & (np.abs(target[:, j]) > tolerance * scale)
        if impossible.any():
            row = int(np.flatnonzero(impossible)[0])
            raise ValueError(f"unreachable {name} target at row {row}: zero response")
        q[~zero, j] = target[~zero, j] / factor[~zero]
    if not np.isfinite(q).all():
        raise ValueError("nonfinite lower-order action controls")

    replay = J * q
    error = np.abs(replay - target) / np.maximum(1.0, np.abs(target))
    max_error = float(np.max(error))
    if max_error > tolerance:
        raise ValueError(f"lower-order control replay failed: {max_error:.3e}")

    action = background[[c for c in ("u", "x", "phi", "f", "h", "X", "A0prime", "S_SVT", "T_H") if c in background]].copy()
    for j, name in enumerate(CONTROL_NAMES):
        action[name] = q[:, j]
    emitted = action[[c for c in ("u", "x") if c in action]].copy()
    for j, name in enumerate(TARGET_NAMES):
        emitted[name] = replay[:, j]

    determinant = J[:, 0] * J[:, 1] * J[:, 2]
    vector_active = np.abs(background.A0prime.to_numpy(float)) > 1e-14
    reduced_det = J[:, 1] * J[:, 2]
    report = {
        "status": "PASS",
        "map": "triangular lower-order normal-jet response",
        "controls": list(CONTROL_NAMES),
        "targets": list(TARGET_NAMES),
        "max_scaled_replay_error": max_error,
        "all_controls_finite": bool(np.isfinite(q).all()),
        "vector_active_rows": int(vector_active.sum()),
        "vector_free_rows": int((~vector_active).sum()),
        "min_abs_full_det_vector_active": (
            float(np.min(np.abs(determinant[vector_active]))) if vector_active.any() else None
        ),
        "min_abs_reduced_c3_e3_det": float(np.min(np.abs(reduced_det))),
        "background_null": False,
        "background_null_reason": (
            "local second-jet response only; a smooth deformation is background-null only after "
            "the common Hessian chain rules H@(phi',X',F')=0 are satisfied"
        ),
        "principal_symbol_unchanged": False,
        "principal_symbol_reason": (
            "holonomic completion generally induces f2XX/f2XF/f2FF, so unchanged principal "
            "symbol cannot be certified from the split lower inverse alone"
        ),
        "endpoint_rule": "when A0prime=0 require v5_target=0 and use reduced c3/e3 block",
    }
    return action, emitted, report
