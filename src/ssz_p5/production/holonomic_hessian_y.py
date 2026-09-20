"""4D background-null f2 Hessian controls for the electric SVT Inner band.

The relevant even-parity action dependence contains f2(phi, X, F, Y).  For an
*additive background-null deformation* whose value and first action jets vanish
on the selected radial curve, its symmetric Hessian H obeys

    H @ (phi', X', F', Y') = 0.

The six transverse entries (XX, XF, XY, FF, FY, YY) therefore determine the
four phi-containing entries.  This module implements that completion and the
Appendix-A f2-only response in (v5,c3,e3,v1,v4,c2).

Important limitation
--------------------
The six transverse action jets are not automatically six independent coefficient
controls on the static electric background.  In particular the five algebraic
channels (v5,c3,v1,v4,c2) can have local rank below five.  Callers must audit
reachability rather than infer it from Hessian dimension counting.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse

from ..jets.jet9d8 import _jet_weights, profile_derivative

TRANSVERSE = ("f2XX", "f2XF", "f2XY", "f2FF", "f2FY", "f2YY")
MIXED = ("f2phiphi", "f2phiX", "f2phiF", "f2phiY")
RESPONSES = ("v5", "c3", "e3", "v1", "v4", "c2")
ALGEBRAIC_RESPONSES = ("v5", "c3", "v1", "v4", "c2")


def _arrays(background: pd.DataFrame):
    required = ("x", "f", "h", "X", "A0prime")
    missing = [name for name in required if name not in background]
    if missing:
        raise KeyError(f"missing background columns: {missing}")
    x, f, h, X, A = (background[name].to_numpy(float) for name in required)
    if not all(np.isfinite(v).all() for v in (x, f, h, X, A)):
        raise ValueError("nonfinite background")
    if np.any(f <= 0) or np.any(h <= 0):
        raise ValueError("requires f>0 and h>0")
    if "phiprime" in background:
        ph = background.phiprime.to_numpy(float)
    else:
        radicand = -2.0 * X / h
        if np.any(radicand < -1e-13):
            raise ValueError("cannot infer real phiprime")
        ph = -np.sqrt(np.maximum(0.0, radicand))
    if not np.isfinite(ph).all() or np.any(np.abs(ph) <= 1e-14):
        raise ValueError("requires finite nonzero phiprime")
    F = h * A**2 / (2.0 * f)
    Y = 4.0 * X * F
    return x, f, h, X, A, ph, F, Y


def tangent4(background: pd.DataFrame, *, window: int = 9, degree: int = 8) -> np.ndarray:
    """Return (phi', X', F', Y') on the radial curve."""
    x, _f, _h, X, _A, ph, F, Y = _arrays(background)
    return np.column_stack(
        [
            ph,
            profile_derivative(x, X, 1, window, degree),
            profile_derivative(x, F, 1, window, degree),
            profile_derivative(x, Y, 1, window, degree),
        ]
    )


def completion_coefficients(
    background: pd.DataFrame, *, window: int = 9, degree: int = 8
) -> dict[str, np.ndarray]:
    """Linear maps from the six transverse jets to phi-containing Hessian jets.

    Each returned array has shape (N,6), in TRANSVERSE column order.
    """
    t = tangent4(background, window=window, degree=degree)
    ph, Xp, Fp, Yp = t.T
    n = len(background)
    phiX = np.zeros((n, 6))
    phiF = np.zeros((n, 6))
    phiY = np.zeros((n, 6))
    # q = (XX,XF,XY,FF,FY,YY)
    phiX[:, 0] = -Xp / ph
    phiX[:, 1] = -Fp / ph
    phiX[:, 2] = -Yp / ph

    phiF[:, 1] = -Xp / ph
    phiF[:, 3] = -Fp / ph
    phiF[:, 4] = -Yp / ph

    phiY[:, 2] = -Xp / ph
    phiY[:, 4] = -Fp / ph
    phiY[:, 5] = -Yp / ph

    phiphi = -(Xp[:, None] * phiX + Fp[:, None] * phiF + Yp[:, None] * phiY) / ph[:, None]
    return {"f2phiphi": phiphi, "f2phiX": phiX, "f2phiF": phiF, "f2phiY": phiY}


def complete_hessian(
    background: pd.DataFrame,
    transverse: pd.DataFrame | np.ndarray,
    *,
    window: int = 9,
    degree: int = 8,
) -> tuple[np.ndarray, pd.DataFrame]:
    """Complete one symmetric 4x4 background-null Hessian from transverse jets."""
    if isinstance(transverse, pd.DataFrame):
        q = transverse[list(TRANSVERSE)].to_numpy(float)
    else:
        q = np.asarray(transverse, float)
    if q.shape != (len(background), 6) or not np.isfinite(q).all():
        raise ValueError("transverse Hessian must have shape (N,6) and be finite")
    maps = completion_coefficients(background, window=window, degree=degree)
    pp = np.einsum("nij,nj->ni", np.stack([maps[k] for k in MIXED], axis=1), q)
    phiphi, phiX, phiF, phiY = pp.T

    H = np.zeros((len(background), 4, 4), dtype=float)
    H[:, 0, 0] = phiphi
    H[:, 0, 1] = H[:, 1, 0] = phiX
    H[:, 0, 2] = H[:, 2, 0] = phiF
    H[:, 0, 3] = H[:, 3, 0] = phiY
    H[:, 1, 1] = q[:, 0]
    H[:, 1, 2] = H[:, 2, 1] = q[:, 1]
    H[:, 1, 3] = H[:, 3, 1] = q[:, 2]
    H[:, 2, 2] = q[:, 3]
    H[:, 2, 3] = H[:, 3, 2] = q[:, 4]
    H[:, 3, 3] = q[:, 5]

    out = background[[c for c in ("u", "x") if c in background]].copy()
    for j, name in enumerate(TRANSVERSE):
        out[name] = q[:, j]
    for values, name in zip((phiphi, phiX, phiF, phiY), MIXED, strict=True):
        out[name] = values
    return H, out


def chain_residual(background: pd.DataFrame, H: np.ndarray, *, window: int = 9, degree: int = 8):
    if H.shape != (len(background), 4, 4):
        raise ValueError("H shape mismatch")
    t = tangent4(background, window=window, degree=degree)
    residual = np.einsum("nij,nj->ni", H, t)
    denom = np.sum(np.abs(H * t[:, None, :]), axis=2)
    normalized = np.divide(np.abs(residual), denom, out=np.zeros_like(residual), where=denom > 0)
    return residual, normalized


def response_from_hessian(
    background: pd.DataFrame,
    H: np.ndarray,
    *,
    window: int = 9,
    degree: int = 8,
) -> pd.DataFrame:
    """Return the f2-only Appendix-A delta in the six audited coefficient slots."""
    if H.shape != (len(background), 4, 4):
        raise ValueError("H shape mismatch")
    x, f, h, _X, A, ph, F, _Y = _arrays(background)
    pp = H[:, 0, 0]
    phiX, phiF, phiY = H[:, 0, 1], H[:, 0, 2], H[:, 0, 3]
    XX, XF, XY = H[:, 1, 1], H[:, 1, 2], H[:, 1, 3]
    FF, FY, YY = H[:, 2, 2], H[:, 2, 3], H[:, 3, 3]

    v5 = x**2 * np.sqrt(h / f) * A * (phiF - 2.0 * h * ph**2 * phiY)
    c3 = 0.5 * x**2 * np.sqrt(f * h) * ph**2 * phiX - 0.5 * x**2 * A**2 * np.sqrt(h / f) * (
        phiF - 4.0 * h * ph**2 * phiY
    )
    current = x**2 * np.sqrt(f * h) * ph * (phiX + 4.0 * F * phiY)
    e3 = -0.5 * profile_derivative(x, current, 1, window, degree) - 0.5 * x**2 * np.sqrt(f / h) * pp
    v1 = (
        x**2
        * h**1.5
        * A**2
        / (2.0 * f**1.5)
        * (FF - 4.0 * h * ph**2 * FY + 4.0 * h**2 * ph**4 * YY)
    )
    v4 = 2.0 * x**2 * h**2.5 * ph * A**3 / f**1.5 * (
        2.0 * h * ph**2 * YY - FY
    ) + h**1.5 * A / np.sqrt(f) * (2.0 * h * ph**3 * x**2 * XY - x**2 * ph * XF)
    c2 = (
        -(x**2) * h**2.5 * ph * A**4 / f**1.5 * (4.0 * h * ph**2 * YY - FY)
        - h**1.5 * A**2 / (2.0 * np.sqrt(f)) * (6.0 * h * ph**3 * x**2 * XY - ph * x**2 * XF)
        - 0.5 * x**2 * np.sqrt(f * h) * h * ph**3 * XX
    )
    out = background[[c for c in ("u", "x") if c in background]].copy()
    for name, values in zip(RESPONSES, (v5, c3, e3, v1, v4, c2), strict=True):
        out[name] = values
    return out


def algebraic_response_matrix(
    background: pd.DataFrame, *, window: int = 9, degree: int = 8
) -> np.ndarray:
    """Return N x 5 x 6 local map to (v5,c3,v1,v4,c2)."""
    x, f, h, _X, A, ph, _F, _Y = _arrays(background)
    maps = completion_coefficients(background, window=window, degree=degree)
    Bx, Bf, By = maps["f2phiX"], maps["f2phiF"], maps["f2phiY"]
    n = len(background)
    M = np.zeros((n, 5, 6))
    M[:, 0, :] = (x**2 * np.sqrt(h / f) * A)[:, None] * (Bf - 2.0 * (h * ph**2)[:, None] * By)
    M[:, 1, :] = (0.5 * x**2 * np.sqrt(f * h) * ph**2)[:, None] * Bx - (
        0.5 * x**2 * A**2 * np.sqrt(h / f)
    )[:, None] * (Bf - 4.0 * (h * ph**2)[:, None] * By)
    pref = x**2 * h**1.5 * A**2 / (2.0 * f**1.5)
    M[:, 2, 3] = pref
    M[:, 2, 4] = -4.0 * h * ph**2 * pref
    M[:, 2, 5] = 4.0 * h**2 * ph**4 * pref

    p1 = 2.0 * x**2 * h**2.5 * ph * A**3 / f**1.5
    p2 = h**1.5 * A / np.sqrt(f)
    M[:, 3, 5] += p1 * 2.0 * h * ph**2
    M[:, 3, 4] -= p1
    M[:, 3, 2] += p2 * 2.0 * h * ph**3 * x**2
    M[:, 3, 1] -= p2 * x**2 * ph

    p1 = -(x**2) * h**2.5 * ph * A**4 / f**1.5
    p2 = -(h**1.5) * A**2 / (2.0 * np.sqrt(f))
    M[:, 4, 5] += p1 * 4.0 * h * ph**2
    M[:, 4, 4] -= p1
    M[:, 4, 2] += p2 * 6.0 * h * ph**3 * x**2
    M[:, 4, 1] -= p2 * ph * x**2
    M[:, 4, 0] += -0.5 * x**2 * np.sqrt(f * h) * h * ph**3
    return M


def derivative_matrix(x: np.ndarray, *, window: int = 9, degree: int = 8) -> sparse.csr_matrix:
    """Sparse matrix exactly matching the accepted JET local derivative weights."""
    inds, weights = _jet_weights(np.asarray(x, float), 1, window, degree)
    n, width = weights.shape
    rows = np.repeat(np.arange(n), width)
    return sparse.csr_matrix((weights.ravel(), (rows, inds.ravel())), shape=(n, n))
