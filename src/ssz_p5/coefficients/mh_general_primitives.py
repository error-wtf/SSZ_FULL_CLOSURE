"""General Maxwell--Horndeski Appendix-A coefficient emitter from action-derived primitives.

This module regenerates the common 41-slot unreduced even-parity coefficient
stream from the primitive on-shell profiles
    a1, c2, c4, F, G, H
plus the fixed background geometry.  It implements the Appendix-A identities
of Kase & Tsujikawa, Phys. Rev. D 107, 104045 (2023), with the project JET9D8
radial derivative convention and the 2026-09-17 lower-order representative
c3=e3=v5=0.

For the punctured Horndeski core used by the locked regional production member
A0'=0 and the vector sector is the Maxwell baseline.  The function keeps A0'
as an argument for algebraic completeness but nonzero-A0' Horndeski+SVT
handover assembly is handled separately at the action-contribution level.
"""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import pandas as pd

SLOTS = [
    *(f"a{i}" for i in range(1, 10)),
    *(f"b{i}" for i in range(1, 6)),
    *(f"c{i}" for i in range(1, 7)),
    *(f"d{i}" for i in range(1, 5)),
    *(f"e{i}" for i in range(1, 5)),
    *(f"v{i}" for i in range(1, 14)),
]


def jet9d8(
    x: Iterable[float], y: Iterable[float], order: int = 1, window: int = 9, degree: int = 8
) -> np.ndarray:
    """Use the existing project derivative service; no second implementation."""
    from ..jets.jet9d8 import derivative

    return derivative(x, y, order, window=window, degree=degree)


def _arr(d: pd.DataFrame, name: str, default=None) -> np.ndarray:
    if name in d:
        return d[name].to_numpy(float)
    if default is None:
        raise ValueError(f"missing required column {name!r}")
    if np.isscalar(default):
        return np.full(len(d), float(default))
    return np.asarray(default, dtype=float)


def emit_from_primitives(
    df: pd.DataFrame, *, window: int = 9, degree: int = 8, regularize_photon_root: bool = False
) -> pd.DataFrame:
    """Emit 41 coefficients from general Horndeski primitive profiles.

    Required columns (aliases accepted for geometry):
      x or r_over_rs, f, h, phi_r or phiprime,
      a1, c2, c4, F_tensor, G_tensor, H_tensor.

    Optional: u, phi, A0prime.  The current punctured-core production use has
    A0prime=0.  ``regularize_photon_root`` is intentionally False by default:
    this emitter never silently fills removable 0/0 expressions.  A caller may
    request continuous local interpolation only after a convergence check.
    """
    d = df.copy().reset_index(drop=True)
    if "x" in d:
        r = _arr(d, "x")
    elif "r_over_rs" in d:
        r = _arr(d, "r_over_rs")
    else:
        raise ValueError("need x or r_over_rs")
    if np.any(np.diff(r) <= 0):
        raise ValueError("emitter requires strictly increasing radial coordinate x")

    f = _arr(d, "f")
    h = _arr(d, "h")
    if "phi_r" in d:
        ph = _arr(d, "phi_r")
    elif "phiprime" in d:
        ph = _arr(d, "phiprime")
    elif "X" in d:
        ph = -np.sqrt(np.maximum(0.0, -2.0 * _arr(d, "X") / h))
    else:
        raise ValueError("need phi_r/phiprime/X")
    Ap = _arr(d, "A0prime", 0.0)

    a1 = _arr(d, "a1")
    c2 = _arr(d, "c2")
    c4 = _arr(d, "c4")
    F = _arr(d, "F_tensor")
    G = _arr(d, "G_tensor")
    H = _arr(d, "H_tensor")

    fp = jet9d8(r, f, 1, window, degree)
    hp = jet9d8(r, h, 1, window, degree)
    phipp = jet9d8(r, ph, 1, window, degree)
    Hp = jet9d8(r, H, 1, window, degree)
    Gp = jet9d8(r, G, 1, window, degree)
    sqfh = np.sqrt(f * h)

    # Maxwell baseline in the Kase--Tsujikawa old-v notation.
    ov1 = 0.5 * r**2 * np.sqrt(h / f)
    ov2 = Ap * ov1
    ov3 = -Ap * ov1
    ov4 = np.zeros_like(r)
    ov5 = np.zeros_like(r)
    ov6 = 0.25 * Ap**2 * ov1
    ov8 = 1.0 / (2.0 * sqfh)
    ov7 = -2.0 * h * Ap * ov8
    ov9 = -f * h * ov8

    a4 = 0.5 * sqfh * H
    a1p = jet9d8(r, a1, 1, window, degree)
    a4p = jet9d8(r, a4, 1, window, degree)

    # A2 and the corrected holonomic A5 identity.  Because v4=v5=0 in the
    # pure Maxwell--Horndeski core, their terms vanish here.
    P = phipp / ph + 0.5 * hp / h
    R = r / ph * (fp / f - hp / h) * a4
    a2 = a1p - P * a1 + R + 0.5 * Ap * ov4
    a3 = -0.5 * ph * a1 - r * a4
    a3p = jet9d8(r, a3, 1, window, degree)
    a6 = -np.sqrt(f) / (2.0 * np.sqrt(h) * ph) * (Hp + H / r - F / r)

    # a5 = a2' - a1'' - (A0' v4/2)' + A0' v5/2, written without taking a
    # second numerical derivative of a1.  This is exactly the locked 17-Sep
    # holonomic form for the Maxwell core.
    Pp = jet9d8(r, P, 1, window, degree)
    Rp = jet9d8(r, R, 1, window, degree)
    a5 = -Pp * a1 - P * a1p + Rp + 0.5 * Ap * ov5

    a7 = a3p - 0.5 * Ap**2 * ov1 - 0.25 * ph * Ap * ov4
    a8 = -a4 / (2.0 * h)
    a9 = a4p + (1.0 / r - 0.5 * fp / f) * a4

    b1 = a4 / (2.0 * f)
    b2 = -2.0 * a1 / f
    b3 = -2.0 * (a2 - a1p) / f + Ap * ov4 / f
    b4 = -2.0 * a3 / f
    b5 = -2.0 * b1

    c1 = -a1 / (f * h)
    c3 = np.zeros_like(r)
    c5 = -h * ph * c4 - 0.5 * sqfh * G / r - 0.5 * fp * a4 / f
    c6 = (
        fp * ph * a1 / (8.0 * f)
        + 0.5 * fp * r * a4 / f
        - 0.25 * ph * c2
        + 0.5 * h * ph * r * c4
        + 0.25 * sqfh * G
        + 0.25 * Ap**2 * ov1
        + 0.125 * ph * Ap * ov4
    )

    d1 = a4 / (2.0 * f)
    d2 = 2.0 * h * c4

    photon = fp * r - 2.0 * f
    # Along the one-dimensional background curve, the holonomic derivative is
    # the locked direct-closure representative for partial_phi a4.
    da4_dphi = a4p / ph
    Q = (
        2.0 * phipp / (h * ph * r)
        + fp**2 / f**2
        - fp * hp / (f * h)
        - 2.0 * fp / (f * r)
        + 2.0 * hp / (h * r)
        + hp / (h**2 * r)
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        d3 = (
            -((2.0 * phipp / ph + hp / h) * a1) / r**2
            + 2.0 * f / (photon * ph) * Q * a4
            + photon / (f * r) * da4_dphi
            + np.sqrt(f) / (ph * np.sqrt(h) * r**2) * F
            - f**1.5
            / (np.sqrt(h) * photon * ph)
            * (fp / (f * r) + 2.0 * phipp / (ph * r) + hp / (h * r) - 2.0 / r**2)
            * G
        )
    d4 = 0.5 * sqfh * G / r**2

    e1 = ((fp / f + 0.5 * hp / h) * a1 - 2.0 * a1p + a2 - 2.0 * r * h * a6 - 0.5 * Ap * ov4) / (
        ph * f * h
    )
    e2 = -(fp * a1 / f + 2.0 * c2 + 4.0 * h * r * c4 + Ap * ov4) / (2.0 * ph)
    e3 = np.zeros_like(r)
    c4p = jet9d8(r, c4, 1, window, degree)
    with np.errstate(divide="ignore", invalid="ignore"):
        e4 = (
            c4p / ph
            - 0.5 * fp * a4p / (f * ph**2 * h)
            - 0.5 * np.sqrt(f) * Gp / (ph**2 * np.sqrt(h) * r)
            + (phipp / ph + 0.5 * hp / h) * a1 / (h * ph * r**2)
            + a4
            / (4.0 * h * ph**2)
            * (
                (photon - 4.0 * f) * fp / (f**2 * r)
                + hp * (photon + 6.0 * f) / (h * r * f)
                - 4.0 * f * (2.0 * phipp * h + hp * ph) / (ph * h**2 * r * photon)
            )
            + 0.5 * hp * c4 / (h * ph)
            - 0.5 * photon * da4_dphi / (f * h * r * ph)
            + 0.5 * (fp * h * r - f) * F / (r**2 * np.sqrt(f) * ph**2 * h**1.5)
            + 0.5
            * np.sqrt(f)
            * G
            / (r * ph**2 * h**1.5)
            * (
                f * (2.0 * phipp * h + hp * ph) / (h * ph * photon)
                + 0.5 * (2.0 * f - fp * h * r) / (f * r)
            )
        )

    if regularize_photon_root:
        # Only fill isolated non-finite removable-root samples from nearby
        # finite values.  This is deliberately opt-in and should be paired
        # with a convergence audit by the caller.
        for arr in (d3, e4):
            bad = ~np.isfinite(arr)
            if bad.any():
                good = ~bad
                if good.sum() < 4:
                    raise RuntimeError("insufficient finite neighbors for photon-root continuation")
                arr[bad] = np.interp(r[bad], r[good], arr[good])

    # Correct P1/K scalar oracle from the same primitives.
    mu = 2.0 * (ph * a1 + 2.0 * r * a4) / sqfh
    Y = f * r**4 * H**4 / (mu**2 * h)
    P1 = h * mu / (2.0 * f * r**2 * H**2) * jet9d8(r, Y, 1, window, degree)
    K = 2.0 * P1 - F

    out = pd.DataFrame(
        {
            "u": _arr(d, "u", 1.0 / r),
            "x": r,
            "phi": _arr(d, "phi", np.nan),
            "f": f,
            "h": h,
            "phiprime": ph,
            "A0prime": Ap,
            "a1": a1,
            "a2": a2,
            "a3": a3,
            "a4": a4,
            "a5": a5,
            "a6": a6,
            "a7": a7,
            "a8": a8,
            "a9": a9,
            "b1": b1,
            "b2": b2,
            "b3": b3,
            "b4": b4,
            "b5": b5,
            "c1": c1,
            "c2": c2,
            "c3": c3,
            "c4": c4,
            "c5": c5,
            "c6": c6,
            "d1": d1,
            "d2": d2,
            "d3": d3,
            "d4": d4,
            "e1": e1,
            "e2": e2,
            "e3": e3,
            "e4": e4,
            # Map old Maxwell-Horndeski v1..v9 to the common 13-slot schema.
            "v1": ov1,
            "v2": ov2,
            "v3": ov3,
            "v4": ov4,
            "v5": ov5,
            "v6": np.zeros_like(r),
            "v7": ov6,
            "v8": ov7,
            "v9": ov8,
            "v10": ov9,
            "v11": np.zeros_like(r),
            "v12": np.zeros_like(r),
            "v13": np.zeros_like(r),
            "F_tensor": F,
            "G_tensor": G,
            "H_tensor": H,
            "mu": mu,
            "P1": P1,
            "K_scalar": K,
            "photon_factor": photon,
        }
    )
    return out


def scaled_rel(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    return np.abs(a - b) / np.maximum(1.0, np.abs(b))
