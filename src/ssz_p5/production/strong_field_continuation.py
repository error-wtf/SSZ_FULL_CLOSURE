"""Predictor/corrector helpers for the action-first P5 strong-field branch.

The corrector treats (f3,f3X,f4,f4X,f4XX,...) as free on-curve action data and
solves the three robust static electric background equations algebraically for
(f2F,f2,f2X):

  JA = 0, E00 = 0, E11 = 0.

Mixed phi jets are then completed holonomically and the 41-slot stream is
absolutely re-emitted.  This module never fits historical perturbation
coefficients.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
import pandas as pd

from ..numerics import module
from .full_action_lower import emit_lower_slots
from .svt_background_eom import evaluate_svt_background, scalar_ode_identity


@dataclass(frozen=True)
class CorrectedAction:
    action: pd.DataFrame
    direct41: pd.DataFrame
    residuals: dict


def _arr(d: pd.DataFrame, name: str, default: float = 0.0) -> np.ndarray:
    if name in d:
        return d[name].to_numpy(float)
    return np.full(len(d), float(default))


def algebraic_background_correct(action: pd.DataFrame) -> pd.DataFrame:
    """Reimpose JA=E00=E11=0 on one static electric SVT action table.

    The geometry and electric background (f,h,phi',A0') are held fixed.  The
    corrector changes only f2F, f2 and f2X, which are exactly the three first
    action jets entering the robust current/metric equations linearly.
    """
    d = action.copy().sort_values("x").reset_index(drop=True)
    zk = module("ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py")
    r = _arr(d, "x")
    f = _arr(d, "f")
    h = _arr(d, "h")
    ph = _arr(d, "phiprime")
    A = _arr(d, "A0prime")
    f3 = _arr(d, "f3")
    f3X = _arr(d, "f3X")
    f4 = _arr(d, "f4")
    f4X = _arr(d, "f4X")
    f4XX = _arr(d, "f4XX")
    tf4 = _arr(d, "tf4")

    fp = zk.dr(r, f, 1, 9, 8)
    hp = zk.dr(r, h, 1, 9, 8)

    # Zero-current electric branch.  This solves the constitutive bracket even
    # where A0' later becomes small, rather than dividing by A0'.
    f2F = -(
        4.0 * r * h * ph * f3
        + 8.0 * (1.0 - h) * f4
        + 2.0 * h**2 * ph**2 * (f4X + 2.0 * tf4)
    ) / r**2

    # E00 = 0 solved for f2.
    quartic00 = 4.0 * (h - 1.0) * f4 - h**2 * ph**2 * (f4X + 2.0 * tf4)
    f2 = (
        r * f * hp
        - f * (1.0 - h)
        + r**2 * h * A**2 * f2F
        + 2.0 * r * h**2 * ph * A**2 * f3
        - h * A**2 * quartic00
    ) / (r**2 * f)

    # E11 = 0 solved for f2X.
    quartic11 = (
        4.0 * (3.0 * h - 1.0) * f4
        - h * (9.0 * h - 4.0) * ph**2 * f4X
        + h**3 * ph**4 * f4XX
        - 10.0 * h**2 * ph**2 * tf4
    )
    numer = (
        r * h * fp
        - f * (1.0 - h)
        - r**2 * f * f2
        + r**2 * h * A**2 * f2F
        + 2.0 * r * h**2 * ph * A**2 * (3.0 * f3 - h * ph**2 * f3X)
        - h * A**2 * quartic11
    )
    denom = r**2 * f * h * ph**2
    if np.min(np.abs(denom)) < 1e-13:
        raise RuntimeError("metric corrector denominator is singular")
    f2X = numer / denom

    d["f2F"] = f2F
    d["f2"] = f2
    d["f2X"] = f2X
    return d


def absolute_reemit(action: pd.DataFrame) -> CorrectedAction:
    """Correct robust background equations and absolutely emit the SVT 41er.

    The accepted exact background-null G2XX lift is reattached through its
    action-derived c2/c6/e2 response when the necessary columns are present.
    """
    d = algebraic_background_correct(action)
    zk = module("ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py")
    r = d.x.to_numpy(float)
    f = d.f.to_numpy(float)
    h = d.h.to_numpy(float)
    X = d.X.to_numpy(float)
    A = d.A0prime.to_numpy(float)
    ph = d.phiprime.to_numpy(float)
    F = h * A * A / (2.0 * f)
    Y = 4.0 * X * F
    Xp = zk.dr(r, X, 1, 9, 8)
    Fp = zk.dr(r, F, 1, 9, 8)
    Yp = zk.dr(r, Y, 1, 9, 8)
    if "f2Y" not in d:
        d["f2Y"] = 0.0
    d["f2phi"] = (
        zk.dr(r, d.f2.to_numpy(float), 1, 9, 8)
        - d.f2X.to_numpy(float) * Xp
        - d.f2F.to_numpy(float) * Fp
        - d.f2Y.to_numpy(float) * Yp
    ) / ph
    lower, completed = emit_lower_slots(d)
    out = zk.emit(
        completed,
        selected_v5=lower.v5.to_numpy(float),
        selected_c3=lower.c3.to_numpy(float),
        selected_e3=lower.e3.to_numpy(float),
        v6_phi_selector="action",
    )
    out["a5"] = (
        zk.dr(r, out.a2.to_numpy(float), 1, 9, 8)
        - zk.dr(r, out.a1.to_numpy(float), 2, 9, 8)
        - zk.dr(r, out.A0prime.to_numpy(float) * out.v4.to_numpy(float) / 2.0, 1, 9, 8)
        + out.A0prime.to_numpy(float) * out.v5.to_numpy(float) / 2.0
    )
    if "delta_c2_G2_lift" in completed:
        dc2 = completed.delta_c2_G2_lift.to_numpy(float)
        out["c2"] = out.c2.to_numpy(float) + dc2
        out["c6"] = out.c6.to_numpy(float) - 0.25 * ph * dc2
        out["e2"] = out.e2.to_numpy(float) - dc2 / ph

    bg = evaluate_svt_background(completed)
    sid = scalar_ode_identity(completed)
    residuals = {
        "max_abs_E00": float(np.max(np.abs(bg.E00))),
        "max_abs_E11": float(np.max(np.abs(bg.E11))),
        "max_abs_JA": float(np.max(np.abs(bg.JA))),
        "max_abs_scalar_identity": float(np.max(np.abs(sid.residual))),
    }
    return CorrectedAction(action=completed, direct41=out, residuals=residuals)
