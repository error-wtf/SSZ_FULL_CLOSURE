"""Physical contract and diagnostics for the P5 right strong-field transition.

The final construction treats 0.70 < u < 0.715 as one action problem.  The
historical seam u=0.71 is not a physical matching surface.  The left anchor is
the electric-hybrid bulk, the right anchor is the pure-Horndeski core, and the
stable inner light ring lies inside the solve.

This module deliberately contains *diagnostics and contract helpers*, not a
coefficient-target optimizer.  The 41 perturbation coefficients are outputs of
an action and are never accepted as independent transition targets here.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from ..reducer.kinetic_schur import kinetic_schur
from .action_handover import compact_bump01

U_LEFT = 0.70
U_RIGHT = 0.715
U_HISTORICAL_SEAM = 0.71
U_INNER_LIGHT_RING = 0.7061345733
REQUIRED_L = (6, 12, 20, 42, 110, 420, 1000)


@dataclass(frozen=True)
class ModePoint:
    L: int
    u: float
    eigenvalues: tuple[float, float, float]
    weakest_vector: tuple[float, float, float]


def transition_coordinate(u) -> np.ndarray:
    """Map the physical transition [0.70,0.715] to s in [0,1]."""
    return (np.asarray(u, float) - U_LEFT) / (U_RIGHT - U_LEFT)


def compact_transition_bump(u) -> np.ndarray:
    """C-infinity compact support on the physical transition interval."""
    return compact_bump01(transition_coordinate(u))


def transition_basis(u, order: int = 4) -> np.ndarray:
    """Small smooth basis for continuation controls.

    Every column is C-infinity-flat at both physical endpoints.  Polynomial
    factors resolve different radial shapes without introducing a fake seam at
    u=0.71.
    """
    if order < 1:
        raise ValueError("order must be >=1")
    s = transition_coordinate(u)
    b = compact_transition_bump(u)
    z = 2.0 * s - 1.0
    return np.column_stack([b * z**k for k in range(order)])


def electric_constitutive_pivot(action: pd.DataFrame) -> np.ndarray:
    """Return chi_E = partial J_A / partial A0' on the on-curve action.

    For the static electric branch

      J_A = sqrt(h/f) A0' C,

    with C containing r^2 f2_F plus F-independent f3/f4 terms.  Holding phi,X
    fixed while varying A0' gives dF/dA0'=h A0'/f and dY/dA0'=4X dF/dA0'.
    Hence f2_F contributes through f2_FF + 4 X f2_FY.  This is the relevant
    local constitutive-root nondegeneracy diagnostic.  It is not a substitute
    for the full perturbation kinetic gate.
    """
    d = action
    r = d.x.to_numpy(float)
    f = d.f.to_numpy(float)
    h = d.h.to_numpy(float)
    ph = d.phiprime.to_numpy(float)
    A = d.A0prime.to_numpy(float)
    X = d.X.to_numpy(float)
    f2F = d.f2F.to_numpy(float)
    f2FF = d.f2FF.to_numpy(float)
    f2FY = d.f2FY.to_numpy(float) if "f2FY" in d else np.zeros(len(d))
    f3 = d.f3.to_numpy(float)
    f4 = d.f4.to_numpy(float)
    f4X = d.f4X.to_numpy(float)
    tf4 = d.tf4.to_numpy(float) if "tf4" in d else np.zeros(len(d))
    C = (
        r**2 * f2F
        + 4.0 * r * h * ph * f3
        + 8.0 * (1.0 - h) * f4
        + 2.0 * h**2 * ph**2 * (f4X + 2.0 * tf4)
    )
    dC_dA = r**2 * (h * A / f) * (f2FF + 4.0 * X * f2FY)
    return np.sqrt(h / f) * (C + A * dC_dA)


def _sym(a: np.ndarray) -> np.ndarray:
    return (a + a.swapaxes(1, 2)) / 2.0


def mode_points(stream: pd.DataFrame, L_values=REQUIRED_L, *, u=U_INNER_LIGHT_RING):
    """Kinetic eigensystem at one physical radius for mode identification."""
    idx = int(np.argmin(np.abs(stream.u.to_numpy(float) - float(u))))
    rows = []
    for L in L_values:
        red = kinetic_schur(stream, int(L))
        w, v = np.linalg.eigh(_sym(red["K"])[idx])
        rows.append(
            ModePoint(
                int(L),
                float(stream.u.iloc[idx]),
                tuple(float(x) for x in w),
                tuple(float(x) for x in v[:, 0]),
            )
        )
    return rows


def track_modes(K: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Track symmetric-matrix eigenmodes by adjacent eigenvector overlap.

    Returns tracked eigenvalues and tracked eigenvectors.  At every step the
    permutation maximizing total absolute overlap with the previous frame is
    selected.  For the 3x3 even sector exhaustive six-permutation matching is
    deterministic and avoids an external assignment dependency.
    """
    import itertools

    Ks = _sym(np.asarray(K, float))
    n = len(Ks)
    vals = np.empty((n, 3))
    vecs = np.empty((n, 3, 3))
    w, v = np.linalg.eigh(Ks[0])
    vals[0], vecs[0] = w, v
    perms = list(itertools.permutations(range(3)))
    for i in range(1, n):
        w, v = np.linalg.eigh(Ks[i])
        overlap = np.abs(vecs[i - 1].T @ v)
        p = max(perms, key=lambda q: sum(overlap[j, q[j]] for j in range(3)))
        vals[i] = w[list(p)]
        vecs[i] = v[:, list(p)]
        # remove arbitrary sign flips to make vector diagnostics continuous
        for j in range(3):
            if np.dot(vecs[i - 1, :, j], vecs[i, :, j]) < 0:
                vecs[i, :, j] *= -1.0
    return vals, vecs


def first_zero_crossing(u, y) -> float | None:
    """Linear estimate of the first sign crossing when ordered in increasing u."""
    u = np.asarray(u, float)
    y = np.asarray(y, float)
    p = np.argsort(u)
    u, y = u[p], y[p]
    for i in range(len(u) - 1):
        if y[i] == 0:
            return float(u[i])
        if y[i] * y[i + 1] < 0:
            t = -y[i] / (y[i + 1] - y[i])
            return float(u[i] + t * (u[i + 1] - u[i]))
    return None


def kinetic_transition_scan(stream: pd.DataFrame, L_values=REQUIRED_L) -> tuple[pd.DataFrame, dict]:
    """Diagnose the unmodified member across 0.70<u<0.71.

    This is intentionally a diagnostic extrapolation, not a certification of
    the hybrid outside its locked 0.62<u<0.70 production bulk.
    """
    u = stream.u.to_numpy(float)
    mask = (u >= U_LEFT) & (u < U_HISTORICAL_SEAM)
    inds = np.flatnonzero(mask)
    records = []
    summary = {}
    for L in L_values:
        red = kinetic_schur(stream, int(L))
        K = _sym(red["K"])
        eig, vec = track_modes(K[inds])
        # true minimum, independently of tracking label
        raw = np.linalg.eigvalsh(K[inds])
        crossing = first_zero_crossing(u[inds], raw[:, 0])
        j = int(np.argmin(np.abs(u[inds] - U_INNER_LIGHT_RING)))
        summary[str(L)] = {
            "first_min_eigenvalue_zero_u": crossing,
            "min_eig_K_transition": float(np.min(raw[:, 0])),
            "light_ring_u_grid": float(u[inds][j]),
            "light_ring_min_eig_K": float(raw[j, 0]),
            "light_ring_weakest_vector": [float(x) for x in np.linalg.eigh(K[inds][j])[1][:, 0]],
            "min_abs_Dh1": float(np.min(np.abs(red["Dh1"][inds]))),
            "min_abs_auxiliary_determinant": float(np.min(np.abs(red["auxiliary_determinant"][inds]))),
        }
        for q, idx in enumerate(inds):
            records.append(
                {
                    "L": int(L),
                    "u": float(u[idx]),
                    "x": float(stream.x.iloc[idx]),
                    "eig_min": float(raw[q, 0]),
                    "eig_mid": float(raw[q, 1]),
                    "eig_max": float(raw[q, 2]),
                    "tracked_0": float(eig[q, 0]),
                    "tracked_1": float(eig[q, 1]),
                    "tracked_2": float(eig[q, 2]),
                    "tracked0_psi": float(vec[q, 0, 0]),
                    "tracked0_dphi": float(vec[q, 1, 0]),
                    "tracked0_V": float(vec[q, 2, 0]),
                }
            )
    return pd.DataFrame(records), summary


def g2xx_kinetic_null_test(stream: pd.DataFrame, action: pd.DataFrame, amplitude: float = 1.0e4) -> dict:
    """Falsification test: can the existing background-null G2XX lift repair K?

    The known G2XX action response changes c2,c6,e2.  We apply a compact
    transition perturbation and compare the kinetic matrix.  A tiny response is
    useful rank information, not a failure of the transition program.
    """
    d0 = stream.sort_values("x").reset_index(drop=True)
    a = action.sort_values("x").reset_index(drop=True)
    if len(d0) != len(a) or not np.allclose(d0.x, a.x, rtol=0, atol=1e-13):
        raise ValueError("action/direct41 grid mismatch")
    bump = compact_transition_bump(d0.u.to_numpy(float))
    ph = d0.phiprime.to_numpy(float)
    dc2 = a.dc2_dG2XX.to_numpy(float) * amplitude * bump
    changed = d0.copy()
    changed["c2"] = changed.c2.to_numpy(float) + dc2
    changed["c6"] = changed.c6.to_numpy(float) - 0.25 * ph * dc2
    changed["e2"] = changed.e2.to_numpy(float) - dc2 / ph
    mask = (d0.u.to_numpy(float) >= U_LEFT) & (d0.u.to_numpy(float) < U_HISTORICAL_SEAM)
    inds = np.flatnonzero(mask)
    result = {}
    for L in REQUIRED_L:
        k0 = _sym(kinetic_schur(d0, L)["K"])[inds]
        k1 = _sym(kinetic_schur(changed, L)["K"])[inds]
        result[str(L)] = {
            "max_abs_delta_K": float(np.max(np.abs(k1 - k0))),
            "max_abs_delta_min_eigenvalue": float(
                np.max(np.abs(np.linalg.eigvalsh(k1)[:, 0] - np.linalg.eigvalsh(k0)[:, 0]))
            ),
        }
    return {
        "amplitude": float(amplitude),
        "shape": "C-infinity compact bump on 0.70<u<0.715",
        "response": result,
        "interpretation": "G2XX is a null/near-null direction for the right-transition even kinetic obstruction if all reported deltas are negligible.",
    }
