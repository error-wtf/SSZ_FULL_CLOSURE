"""Static electric U(1)-SVT background equations and numerical audit helpers.

The formulas implemented here are the static, spherically symmetric,
purely-electric background equations of Heisenberg & Tsujikawa (2018),
Eqs. (2.14)--(2.20), specialized to M_Pl^2=1 and to the current production
member's f2Y=tilde-f4=0 branch.

Important numerical point
-------------------------
The scalar equation is J_phi' = P_phi.  Directly differentiating J_phi is a
high-cancellation diagnostic on the present radial tables.  Therefore the
release audit keeps two logically separate statements:

* E00, E11 and the zero-current electric branch JA=0 are strict residual gates.
* The direct finite-difference scalar residual is a *resolution diagnostic*.
  It is compared against the archived integrable SVT reference with the exact
  same derivative service.  This comparison does not relax the global
  background_residual_abs release policy and does not by itself certify the
  scalar equation at machine precision.

This distinction prevents a known differentiation floor from being mistaken
for an off-shell physical failure while still keeping the strict scalar
certificate open until an algebraic/identity evaluator is reproduced.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict

import numpy as np
import pandas as pd

from ..jets.jet9d8 import profile_derivative


@dataclass(frozen=True)
class SVTBackgroundResiduals:
    E00: np.ndarray
    E11: np.ndarray
    JA: np.ndarray
    JA_prime: np.ndarray
    Jphi: np.ndarray
    Pphi: np.ndarray
    Ephi_direct: np.ndarray


def _col(d: pd.DataFrame, name: str, default: float = 0.0) -> np.ndarray:
    if name in d.columns:
        return d[name].to_numpy(float)
    return np.full(len(d), float(default))


def evaluate_svt_background(
    d: pd.DataFrame,
    *,
    window: int = 9,
    degree: int = 8,
    require_f2y_zero: bool = True,
) -> SVTBackgroundResiduals:
    """Evaluate the electric static SVT background equations on ``d``.

    ``d`` must carry the action jets used by the Central production member.
    The published compact background equations are written for f2(phi,X,F).
    The current Electric-Hybrid Central branch has f2Y=0, so that specialization
    is exact here.  Nonzero f2Y is rejected rather than silently mapped onto a
    different convention.
    """
    r = _col(d, "x")
    f = _col(d, "f")
    h = _col(d, "h")
    ph = _col(d, "phiprime")
    ap = _col(d, "A0prime")

    f2 = _col(d, "f2")
    f2x = _col(d, "f2X")
    f2f = _col(d, "f2F")
    f2y = _col(d, "f2Y")
    f2phi = _col(d, "f2phi")
    f3 = _col(d, "f3")
    f3x = _col(d, "f3X")
    f3phi = _col(d, "f3phi")
    f4 = _col(d, "f4")
    f4x = _col(d, "f4X")
    f4xx = _col(d, "f4XX")
    f4phi = _col(d, "f4phi")
    f4xphi = _col(d, "f4phiX")
    tf4 = _col(d, "tf4")
    tf4phi = _col(d, "tf4phi")

    if require_f2y_zero and np.max(np.abs(f2y)) > 1e-13:
        raise ValueError(
            "Heisenberg--Tsujikawa compact background evaluator is locked to "
            "the current f2Y=0 production branch"
        )

    def dr(y: np.ndarray, order: int = 1) -> np.ndarray:
        return profile_derivative(r, np.asarray(y, float), order, window, degree)

    fp = dr(f)
    hp = dr(h)

    # HT2018 Eqs. (2.14), (2.15), M_Pl^2=1.
    e00 = r * f * hp - (
        f * (1.0 - h)
        + r**2 * (f * f2 - h * ap**2 * f2f)
        - 2.0 * r * h**2 * ph * ap**2 * f3
        + h * ap**2 * (4.0 * (h - 1.0) * f4 - h**2 * ph**2 * (f4x + 2.0 * tf4))
    )
    e11 = r * h * fp - (
        f * (1.0 - h)
        + r**2 * (f * f2 + f * h * ph**2 * f2x - h * ap**2 * f2f)
        - 2.0 * r * h**2 * ph * ap**2 * (3.0 * f3 - h * ph**2 * f3x)
        + h
        * ap**2
        * (
            4.0 * (3.0 * h - 1.0) * f4
            - h * (9.0 * h - 4.0) * ph**2 * f4x
            + h**3 * ph**4 * f4xx
            - 10.0 * h**2 * ph**2 * tf4
        )
    )

    # HT2018 Eq. (2.20).  The physical electric branch used by P5 fixes the
    # integration constant to zero, so JA itself is the robust branch gate.
    ja = np.sqrt(h / f) * ap * (
        r**2 * f2f
        + 4.0 * r * h * ph * f3
        + 8.0 * (1.0 - h) * f4
        + 2.0 * h**2 * ph**2 * (f4x + 2.0 * tf4)
    )

    # HT2018 Eqs. (2.18), (2.19): scalar current and source.
    jphi = -np.sqrt(h / f) * (
        r**2 * f * f2x * ph
        - 2.0 * h * ap**2 * (2.0 * h * tf4 + (3.0 * h - 2.0) * f4x) * ph
        + 2.0 * r * h**2 * ap**2 * f3x * ph**2
        + h**3 * ap**2 * f4xx * ph**3
        - 2.0 * r * h * ap**2 * f3
    )
    pphi = 1.0 / np.sqrt(f * h) * (
        r**2 * f * f2phi
        + h
        * ap**2
        * (
            4.0 * f4phi
            + 2.0 * h * (r * ph * f3phi - 2.0 * f4phi)
            + h**2 * (f4xphi + 2.0 * tf4phi) * ph**2
        )
    )
    ephi = dr(jphi) - pphi
    return SVTBackgroundResiduals(
        E00=e00,
        E11=e11,
        JA=ja,
        JA_prime=dr(ja),
        Jphi=jphi,
        Pphi=pphi,
        Ephi_direct=ephi,
    )


def residual_metrics(values: SVTBackgroundResiduals, mask: np.ndarray) -> Dict[str, Dict[str, float]]:
    mask = np.asarray(mask, bool)
    out: Dict[str, Dict[str, float]] = {}
    for name in ("E00", "E11", "JA", "JA_prime", "Jphi", "Pphi", "Ephi_direct"):
        a = np.asarray(getattr(values, name), float)[mask]
        out[name] = {
            "max_abs": float(np.max(np.abs(a))),
            "median_abs": float(np.median(np.abs(a))),
        }
    scale = np.abs(values.Pphi[mask]) + np.abs(
        values.Ephi_direct[mask] + values.Pphi[mask]
    )
    out["Ephi_direct"]["max_scaled"] = float(
        np.max(np.abs(values.Ephi_direct[mask]) / np.maximum(1.0, scale))
    )
    return out
