"""Action-level handover diagnostics for the electric-hybrid P5 member.

This module deliberately separates three logically different statements:

1. ``flat_step01`` / ``compact_bump01`` provide an exact C-infinity partition
   with all endpoint derivatives flat.  They are safe glue *coordinates*.
2. ``one_control_a2_probe`` is a falsification test.  It asks whether a
   prescribed smooth v6 bridge can be supported by the central electric A2
   relation using f4 alone.  The accepted 19-Sep checkpoint says no; this
   function makes that failure reproducible instead of hiding it.
3. ``compact_horndeski_control_rank`` checks the six primitive
   Maxwell--Horndeski response directions (a1,c2,c4,F,G,H) with compact
   support.  It is a rank/conditioning diagnostic only.  It is *not* a
   background EOM solve and therefore cannot certify a handover by itself.

The remaining production gate after this module is the joint background-
constrained Horndeski+SVT solve in the two buffers 0.61<u<0.62 and
0.70<u<0.71.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

from ..coefficients.mh_general_primitives import SLOTS, emit_from_primitives
from ..numerics import module
from .electric_hybrid_onshell_central import CentralBuild, build_onshell_central

PRIMITIVE_CONTROLS = ("a1", "c2", "c4", "F_tensor", "G_tensor", "H_tensor")


@dataclass(frozen=True)
class HandoverBand:
    name: str
    lo: float
    hi: float
    orientation: str
    old_background: Path
    old_direct: Path


def handover_bands(root: Path) -> tuple[HandoverBand, HandoverBand]:
    root = Path(root)
    snap = root / "archive/full_working_snapshot"
    return (
        HandoverBand(
            "outer",
            0.61,
            0.62,
            "old_to_central",
            snap / "ssz_p5_F2_outer_same_action_RESOLVED_background_jets_2026-09-15.csv",
            root / "data/prestaged/direct41/outer_resolved_raw_41of41.csv",
        ),
        HandoverBand(
            "inner",
            0.70,
            0.71,
            "central_to_old",
            snap / "ssz_p5_F2_inner_same_action_RESOLVED_v2_background_jets_2026-09-15.csv",
            root / "data/prestaged/direct41/inner_selected_candidate_41of41.csv",
        ),
    )


def flat_step01(t) -> np.ndarray:
    """C-infinity step: exactly 0 for t<=0 and 1 for t>=1.

    On 0<t<1 this is b(t)/(b(t)+b(1-t)) with b(t)=exp(-1/t).
    The logistic form below is algebraically identical and numerically stable.
    """
    t = np.asarray(t, dtype=float)
    out = np.zeros_like(t)
    out[t >= 1.0] = 1.0
    m = (t > 0.0) & (t < 1.0)
    if np.any(m):
        z = 1.0 / t[m] - 1.0 / (1.0 - t[m])
        z = np.clip(z, -700.0, 700.0)
        out[m] = 1.0 / (1.0 + np.exp(z))
    return out


def flat_step_derivative_u(u, lo: float, hi: float) -> tuple[np.ndarray, np.ndarray]:
    u = np.asarray(u, dtype=float)
    t = (u - lo) / (hi - lo)
    w = flat_step01(t)
    wp = np.zeros_like(u)
    m = (t > 0.0) & (t < 1.0)
    if np.any(m):
        wp[m] = (
            w[m]
            * (1.0 - w[m])
            * (1.0 / t[m] ** 2 + 1.0 / (1.0 - t[m]) ** 2)
            / (hi - lo)
        )
    return w, wp


def compact_bump01(t) -> np.ndarray:
    """Unit-height C-infinity compact bump generated from the flat step."""
    s = flat_step01(t)
    return 4.0 * s * (1.0 - s)


def _pchip_u(df: pd.DataFrame, column: str) -> PchipInterpolator:
    q = df.sort_values("u")
    return PchipInterpolator(q.u.to_numpy(float), q[column].to_numpy(float), extrapolate=True)


def _old_augmented(background_path: Path) -> pd.DataFrame:
    d = pd.read_csv(background_path).sort_values("u").reset_index(drop=True)
    ph = -np.sqrt(np.maximum(0.0, -2.0 * d.X.to_numpy(float) / d.h.to_numpy(float)))
    out = d.copy()
    out["phiprime"] = ph
    out["v6_from_action"] = (
        2.0
        * out.h.to_numpy(float) ** 1.5
        * out.A0prime.to_numpy(float)
        / (out.x.to_numpy(float) * np.sqrt(out.f.to_numpy(float)))
        * (
            out.x.to_numpy(float) * ph * out.f3.to_numpy(float)
            - 4.0 * out.f4.to_numpy(float)
            + out.h.to_numpy(float) * ph**2 * out.N4.to_numpy(float)
        )
    )
    return out


def endpoint_action_diagnostics(root: Path, build: CentralBuild | None = None) -> dict:
    """Quantify old/new endpoint mismatches before any handover solve."""
    root = Path(root)
    build = build or build_onshell_central(root)
    ca = build.action.sort_values("u").reset_index(drop=True)
    rows = {}
    for band in handover_bands(root):
        old = _old_augmented(band.old_background)
        edge = band.lo if band.name == "outer" else band.hi
        old_edge = band.lo if band.name == "outer" else band.hi
        central_edge = band.lo if band.name == "outer" else min(band.hi, float(ca.u.max()))
        oldv = float(_pchip_u(old, "v6_from_action")(old_edge))
        cv = float(_pchip_u(ca, "v6_A2_resolved")(central_edge))
        rows[band.name] = {
            "old_edge_u": float(old_edge),
            "central_edge_u": float(central_edge),
            "old_v6": oldv,
            "central_v6": cv,
            "delta_v6": cv - oldv,
            "old_f4": float(_pchip_u(old, "f4")(old_edge)),
            "central_f4": float(_pchip_u(ca, "f4")(central_edge)),
            "old_f2": float(_pchip_u(old, "f2")(old_edge)),
            "central_f2": float(_pchip_u(ca, "f2")(central_edge)),
            "old_f2X": float(_pchip_u(old, "f2X")(old_edge)),
            "central_f2X": float(_pchip_u(ca, "f2X")(central_edge)),
            "old_f2phi": float(_pchip_u(old, "f2phi")(old_edge)),
            "central_f2phi": float(_pchip_u(ca, "f2phi")(central_edge)),
        }
    return rows


def _bridge_geometry(root: Path, band: HandoverBand, build: CentralBuild) -> pd.DataFrame:
    """Geometry on a monotone-x 501-point buffer grid.

    The P5 background geometry is shared by the old on-shell handovers and the
    new Central member.  We interpolate from their union to avoid extrapolating
    across the inner endpoint u=0.71.
    """
    ca = build.action.copy()
    old = _old_augmented(band.old_background)
    common = ["u", "x", "phi", "f", "h", "phiprime", "A0prime"]
    # Old tables store X rather than phiprime; _old_augmented supplied it.
    union = pd.concat([ca[common], old[common]], ignore_index=True)
    union = union.sort_values("u").drop_duplicates("u", keep="last")
    u = np.linspace(band.lo, band.hi, 501)
    vals = {c: _pchip_u(union, c)(u) for c in common if c != "u"}
    g = pd.DataFrame({"u": u, **vals})
    # Appendix emitter differentiates in increasing x, so return that order.
    return g.sort_values("x").reset_index(drop=True)


def compact_horndeski_control_rank(root: Path, band_name: str, build: CentralBuild | None = None) -> dict:
    """Rank of six compact-support primitive response directions.

    This checks the coefficient-level response bundle, not the background EOM.
    Every perturbation is C-infinity-flat at both buffer endpoints.
    """
    root = Path(root)
    build = build or build_onshell_central(root)
    band = next(b for b in handover_bands(root) if b.name == band_name)
    g = _bridge_geometry(root, band, build)
    u = g.u.to_numpy(float)
    t = (u - band.lo) / (band.hi - band.lo)
    bump = compact_bump01(t)

    ref = g.copy()
    ref["a1"] = 0.0
    ref["c2"] = 0.0
    ref["c4"] = 0.0
    ref["F_tensor"] = 0.0
    ref["G_tensor"] = 0.0
    ref["H_tensor"] = 1.0
    e0 = emit_from_primitives(ref, regularize_photon_root=True)

    columns = []
    norms = {}
    max_endpoint_response = {}
    for control in PRIMITIVE_CONTROLS:
        d = ref.copy()
        d[control] = d[control].to_numpy(float) + bump
        e1 = emit_from_primitives(d, regularize_photon_root=True)
        delta = np.column_stack(
            [e1[s].to_numpy(float) - e0[s].to_numpy(float) for s in SLOTS]
        )
        # Do not let interpolation around the removable photon root dominate
        # the algebraic rank.  Non-finite samples are excluded jointly below.
        vec = delta.reshape(-1)
        columns.append(vec)
        norms[control] = float(np.linalg.norm(vec[np.isfinite(vec)]))
        endpoints = np.r_[delta[0], delta[-1]]
        max_endpoint_response[control] = float(np.nanmax(np.abs(endpoints)))

    M = np.column_stack(columns)
    finite_rows = np.all(np.isfinite(M), axis=1)
    M = M[finite_rows]
    colnorm = np.linalg.norm(M, axis=0)
    if np.any(colnorm == 0):
        rank = int(np.linalg.matrix_rank(M))
        singular = []
        cond = float("inf")
    else:
        Mn = M / colnorm
        singular = np.linalg.svd(Mn, compute_uv=False)
        tol = max(Mn.shape) * np.finfo(float).eps * singular[0]
        rank = int(np.sum(singular > tol))
        cond = float(singular[0] / singular[-1]) if singular[-1] > 0 else float("inf")
        singular = [float(v) for v in singular]

    return {
        "band": band.name,
        "u_interval": [band.lo, band.hi],
        "grid_rows": int(len(g)),
        "finite_response_rows": int(M.shape[0]),
        "controls": list(PRIMITIVE_CONTROLS),
        "rank": rank,
        "expected_rank": len(PRIMITIVE_CONTROLS),
        "condition_normalized": cond,
        "singular_values_normalized": singular,
        "response_norms": norms,
        "max_abs_endpoint_response": max_endpoint_response,
        "endpoint_flat_by_construction": True,
        "background_eom_certified": False,
    }


def one_control_a2_probe(root: Path, band_name: str, build: CentralBuild | None = None) -> dict:
    """Rejected one-control bridge: prescribe v6 and solve the central A2 relation for f4.

    A small A2 residual is *not* success here.  The diagnostic rejects the
    branch if the required f4 leaves the endpoint action values by O(1).  This
    is exactly the failure that motivates the joint Horndeski+SVT solve.
    """
    root = Path(root)
    build = build or build_onshell_central(root)
    band = next(b for b in handover_bands(root) if b.name == band_name)
    ca = build.action.sort_values("u").reset_index(drop=True)
    old = _old_augmented(band.old_background)
    zk = module("ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py")

    # Use the Central grid in the open buffer; add exact endpoints for reporting
    # through interpolation.  This is the same derivative convention as the
    # accepted Central A2 construction.
    mask = (ca.u.to_numpy(float) >= band.lo) & (ca.u.to_numpy(float) <= band.hi)
    u = ca.loc[mask, "u"].to_numpy(float)
    if band.name == "inner" and (len(u) == 0 or u.max() < band.hi):
        u = np.r_[u, band.hi]
    u = np.unique(np.sort(u))

    def pc(col):
        return _pchip_u(ca, col)

    def po(col):
        return _pchip_u(old, col)

    r = pc("x")(u)
    f = pc("f")(u)
    h = pc("h")(u)
    ph = pc("phiprime")(u)
    ap = pc("A0prime")(u)

    # JET9D8 background derivatives on the accepted full Central grid, then
    # interpolate to the buffer.  Inner endpoint is a tiny smooth extrapolation.
    rr = ca.x.to_numpy(float)
    fp = _pchip_u(pd.DataFrame({"u": ca.u, "z": zk.dr(rr, ca.f.to_numpy(float), 1, 9, 8)}), "z")(u)
    fpp = _pchip_u(pd.DataFrame({"u": ca.u, "z": zk.dr(rr, ca.f.to_numpy(float), 2, 9, 8)}), "z")(u)
    hp = _pchip_u(pd.DataFrame({"u": ca.u, "z": zk.dr(rr, ca.h.to_numpy(float), 1, 9, 8)}), "z")(u)
    app = _pchip_u(pd.DataFrame({"u": ca.u, "z": zk.dr(rr, ca.A0prime.to_numpy(float), 1, 9, 8)}), "z")(u)

    w, wu = flat_step_derivative_u(u, band.lo, band.hi)
    cv6, ov6 = pc("v6_A2_resolved"), po("v6_from_action")
    cf4x, of4x = pc("f4X"), po("N4")
    cf3x, of3x = pc("f3X"), po("f3X")
    if band.orientation == "old_to_central":
        va, vb = ov6(u), cv6(u)
        dva, dvb = ov6.derivative()(u), cv6.derivative()(u)
        f4x = (1.0 - w) * of4x(u) + w * cf4x(u)
        f3x = (1.0 - w) * of3x(u) + w * cf3x(u)
    else:
        va, vb = cv6(u), ov6(u)
        dva, dvb = cv6.derivative()(u), ov6.derivative()(u)
        f4x = (1.0 - w) * cf4x(u) + w * of4x(u)
        f3x = (1.0 - w) * cf3x(u) + w * of3x(u)
    v6 = (1.0 - w) * va + w * vb
    dvdu = (1.0 - w) * dva + w * dvb + wu * (vb - va)
    v6p = -u**2 * dvdu

    def predict(F4):
        c = 2.0 * h**1.5 * ap / (r * np.sqrt(f))
        F3 = (v6 / c + 4.0 * F4 - h * ph**2 * f4x) / (r * ph)
        F2F = -(4.0 * r * h * ph * F3 + 8.0 * (1.0 - h) * F4 + 2.0 * h**2 * ph**2 * f4x) / r**2
        V10 = -np.sqrt(f * h) / (2.0 * r) * (
            r * F2F + 2.0 * h * ph * F3 + (h * fp / f) * (r * ph * F3 - 4.0 * F4 + h * ph**2 * f4x)
        )
        A7 = (1.0 - (4.0 * h * ap**2 / f) * F4) / (4.0 * r**2 * np.sqrt(f * h))
        A4 = np.sqrt(f * h) / 2.0
        num = A4 * (
            2.0 * f**2 * (r * hp + 2.0 * h)
            + h * r**2 * fp**2
            - f * r * (r * fp * hp + 2.0 * h * (r * fpp + fp))
        ) - f * h * r**2 * (ap * (4.0 * V10 * ap + v6 * fp) + f * v6 * app + 8.0 * A7 * f**2)
        den = f**2 * h * r**2 * ap
        return num / den, F3, F2F

    p0, _, _ = predict(np.zeros_like(u))
    p1, _, _ = predict(np.ones_like(u))
    slope = p1 - p0
    solved_f4 = (v6p - p0) / slope
    pred, F3, F2F = predict(solved_f4)

    expected_left = po("f4")(band.lo) if band.orientation == "old_to_central" else pc("f4")(band.lo)
    expected_right = pc("f4")(band.hi) if band.orientation == "old_to_central" else po("f4")(band.hi)
    endpoint_error = max(abs(float(solved_f4[0] - expected_left)), abs(float(solved_f4[-1] - expected_right)))
    accepted = bool(endpoint_error < 1e-3 and np.nanmax(np.abs(pred - v6p)) < 1e-8)
    return {
        "band": band.name,
        "u_interval": [band.lo, band.hi],
        "rows": int(len(u)),
        "max_abs_A2_residual": float(np.nanmax(np.abs(pred - v6p))),
        "min_abs_f4_sensitivity": float(np.nanmin(np.abs(slope))),
        "solved_f4_endpoints": [float(solved_f4[0]), float(solved_f4[-1])],
        "expected_f4_endpoints": [float(expected_left), float(expected_right)],
        "max_abs_f4_endpoint_error": float(endpoint_error),
        "solved_f4_range": [float(np.nanmin(solved_f4)), float(np.nanmax(solved_f4))],
        "solved_f3_range": [float(np.nanmin(F3)), float(np.nanmax(F3))],
        "solved_f2F_range": [float(np.nanmin(F2F)), float(np.nanmax(F2F))],
        "f3X_bridge_range": [float(np.nanmin(f3x)), float(np.nanmax(f3x))],
        "accepted": accepted,
        "interpretation": "REJECT_ONE_CONTROL" if not accepted else "PASS_ONE_CONTROL",
    }
