#!/usr/bin/env python3
"""Projected on-shell controllability audit for the right strong-field transition.

This audit implements the local test

    B dp = 0,     B = d(E00,E11,JA,Ephi_id)/dp,
    d lambda_min = g . dp,

using *action-level* local C-infinity perturbations and fresh absolute 41-slot
reemission.  Historical coefficient profiles are never targets.

Two nested tangent spaces are reported:
  * strict metric+vector tangent: ker d(E00,E11,JA), whose three equations are
    already strict background gates in the project;
  * extended scalar-identity tangent: ker d(E00,E11,JA,Ephi_id), where Ephi_id
    is the algebraic first-order scalar identity evaluator.  This second result
    remains diagnostic until the scalar identity itself is promoted globally.

The primary coordinate-invariant existence test is

    rank([B; g]) > rank(B),

which is equivalent to P_ker(B) g != 0.  A small linear program additionally
asks whether one native-coordinate tangent can increase the weakest kinetic
mode for all required finite L simultaneously.

No result from this file is a physical transition certificate.  It is the
pre-continuation rank/control gate requested by the strong-field contract.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from ssz_p5.config import SLOT_NAMES
from ssz_p5.numerics import module
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central
from ssz_p5.production.full_action_lower import emit_lower_slots
from ssz_p5.production.strong_field_transition import (
    REQUIRED_L,
    U_INNER_LIGHT_RING,
    U_LEFT,
)
from ssz_p5.production.svt_background_eom import evaluate_svt_background, scalar_ode_identity
from ssz_p5.reducer.kinetic_schur import kinetic_schur

OUTDIR = ROOT / "data/generated/strong_field_transition_2026-09-20"
REPORT = OUTDIR / "TRANSITION_PROJECTED_CONTROLLABILITY.json"
TABLE = OUTDIR / "TRANSITION_PROJECTED_CONTROLLABILITY.csv"

# These are genuine on-curve action profiles, not 41-slot coefficients.
SVT_DIRECTIONS = (
    "f2",
    "f2X",
    "f2F",
    "f3",
    "f3X",
    "f4",
    "f4X",
    "f4XX",
)
# Exact Horndeski normal sections whose first background variation vanishes.
EXACT_NULL_DIRECTIONS = ("H_G2XX", "H_G3XX_TUBULAR")
DIRECTIONS = SVT_DIRECTIONS + EXACT_NULL_DIRECTIONS

SAMPLE_POINTS = (0.7013, 0.7030, U_INNER_LIGHT_RING, 0.7080)
REQUESTED_RIGHT_DIAGNOSTIC = 0.7120
FD_EPS = 1.0e-5
LS = (6, 20, 42, 1000)


def _sym(a: np.ndarray) -> np.ndarray:
    return (a + a.swapaxes(1, 2)) / 2.0


def _local_bump(u: np.ndarray, u0: float, width: float) -> np.ndarray:
    z = (np.asarray(u, float) - float(u0)) / float(width)
    out = np.zeros_like(z)
    m = np.abs(z) < 1.0
    out[m] = np.exp(-1.0 / (1.0 - z[m] ** 2)) / np.exp(-1.0)
    return out


def _width(u0: float) -> float:
    # Compact support stays inside the available absolute Hybrid action table
    # and away from the left physical anchor.  The final point uses a smaller
    # support because the current absolute Hybrid emitter stops just below .71.
    if u0 >= 0.7088:
        return 5.0e-4
    if u0 <= 0.7015:
        return 7.0e-4
    return 8.0e-4


def _recomplete_and_emit(action: pd.DataFrame):
    """Holonomically complete and absolutely re-emit one SVT action table.

    The existing background-null Horndeski G2XX lift is then reattached through
    its exact action-derived c2/c6/e2 response.  This is the same decomposition
    used by the accepted Electric-Hybrid checkpoint, not a historical slot fit.
    """
    d = action.copy().sort_values("x").reset_index(drop=True)
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
    # f2phi is a first action jet and must follow the on-curve chain rule when
    # f2/f2X/f2F are varied.  Never keep the historical value frozen.
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
    # Locked holonomic a5 identity.
    out["a5"] = (
        zk.dr(r, out.a2.to_numpy(float), 1, 9, 8)
        - zk.dr(r, out.a1.to_numpy(float), 2, 9, 8)
        - zk.dr(r, out.A0prime.to_numpy(float) * out.v4.to_numpy(float) / 2.0, 1, 9, 8)
        + out.A0prime.to_numpy(float) * out.v5.to_numpy(float) / 2.0
    )
    # Existing exact background-null Horndeski G2XX lift.
    dc2 = d.delta_c2_G2_lift.to_numpy(float)
    out["c2"] = out.c2.to_numpy(float) + dc2
    out["c6"] = out.c6.to_numpy(float) - 0.25 * ph * dc2
    out["e2"] = out.e2.to_numpy(float) - dc2 / ph
    return completed, out


def _eom_vector(completed: pd.DataFrame, i: int) -> np.ndarray:
    e = evaluate_svt_background(completed)
    sid = scalar_ode_identity(completed)
    return np.array([e.E00[i], e.E11[i], e.JA[i], sid.residual[i]], dtype=float)


def _ranks_and_projection(B: np.ndarray, g: np.ndarray) -> dict:
    """Coordinate-robust rank test plus a native-coordinate projection."""
    B = np.asarray(B, float)
    g = np.asarray(g, float)
    # Row-normalize; row scaling does not alter the nullspace.
    rn = np.linalg.norm(B, axis=1)
    keep = rn > max(1e-14, 1e-12 * (np.max(rn) if rn.size else 0.0))
    Br = B[keep].copy()
    if len(Br):
        Br /= np.linalg.norm(Br, axis=1)[:, None]
    # Column equilibration improves numerical rank only; it is an invertible
    # reparameterization on nonzero response columns and therefore does not
    # alter the exact structural controllability statement.
    col = np.sqrt(np.sum(Br * Br, axis=0) + g * g)
    col = np.where(col > 1e-14, col, 1.0)
    Be = Br / col[None, :]
    ge = g / col

    def rank(A):
        if A.size == 0:
            return 0, []
        s = np.linalg.svd(A, compute_uv=False)
        tol = max(A.shape) * np.finfo(float).eps * (s[0] if len(s) else 1.0) * 50.0
        return int(np.sum(s > tol)), [float(x) for x in s]

    rB, sB = rank(Be)
    grow = ge / max(np.linalg.norm(ge), 1e-300)
    rA, sA = rank(np.vstack([Be, grow]))
    # Native-coordinate orthogonal projection is secondary/diagnostic; the
    # invariant existence statement is the rank increase above.
    if len(Br):
        U, S, Vt = np.linalg.svd(Br, full_matrices=True)
        tol = max(Br.shape) * np.finfo(float).eps * (S[0] if len(S) else 1.0) * 50.0
        rr = int(np.sum(S > tol))
        N = Vt[rr:].T
    else:
        N = np.eye(B.shape[1])
    proj = N @ (N.T @ g)
    return {
        "rank_B": rB,
        "rank_B_plus_g": rA,
        "structurally_controllable": bool(rA > rB),
        "singular_values_B_equilibrated": sB,
        "singular_values_augmented_equilibrated": sA,
        "native_gradient_norm": float(np.linalg.norm(g)),
        "native_projected_gradient_norm": float(np.linalg.norm(proj)),
        "native_projected_fraction": float(np.linalg.norm(proj) / max(np.linalg.norm(g), 1e-300)),
        "native_projected_gradient": proj,
    }


def _common_gain_lp(B: np.ndarray, gradients: dict[int, np.ndarray]) -> dict:
    """Maximize the minimum first-order K gain for all required L.

    Native action amplitudes are bounded by |dp_j|<=1.  The magnitude is not a
    physical norm; only positivity of the optimum is used as an existence
    diagnostic.  A later continuation solver must supply its own trust metric.
    """
    B = np.asarray(B, float)
    rn = np.linalg.norm(B, axis=1)
    keep = rn > max(1e-14, 1e-12 * (np.max(rn) if rn.size else 0.0))
    Beq = B[keep].copy()
    if len(Beq):
        Beq /= np.linalg.norm(Beq, axis=1)[:, None]
    n = B.shape[1]
    # x=(dp_1...dp_n,t), maximize t => minimize -t.
    c = np.zeros(n + 1)
    c[-1] = -1.0
    Aub = []
    bub = []
    for L in LS:
        row = np.zeros(n + 1)
        row[:n] = -gradients[L]
        row[-1] = 1.0
        Aub.append(row)
        bub.append(0.0)
    Aeq = None
    beq = None
    if len(Beq):
        Aeq = np.column_stack([Beq, np.zeros(len(Beq))])
        beq = np.zeros(len(Beq))
    bounds = [(-1.0, 1.0)] * n + [(None, None)]
    sol = linprog(c, A_ub=np.asarray(Aub), b_ub=np.asarray(bub), A_eq=Aeq, b_eq=beq,
                  bounds=bounds, method="highs")
    if not sol.success:
        return {"success": False, "message": sol.message}
    dp = sol.x[:n]
    gains = {str(L): float(gradients[L] @ dp) for L in LS}
    return {
        "success": True,
        "max_common_first_order_gain_native_box": float(sol.x[-1]),
        "all_required_L_gain_positive": bool(min(gains.values()) > 0.0),
        "gains": gains,
        "direction": dp,
        "max_abs_native_component": float(np.max(np.abs(dp))),
    }


def _point_audit(full_action: pd.DataFrame, full_stream: pd.DataFrame, u0: float) -> tuple[dict, list[dict]]:
    width = _width(u0)
    guard = max(0.0024, 3.0 * width)
    hi = min(float(full_action.u.max()), u0 + guard)
    lo = max(float(full_action.u.min()), u0 - guard)
    loc = full_action[(full_action.u >= lo) & (full_action.u <= hi)].copy().sort_values("x").reset_index(drop=True)
    if len(loc) < 35:
        raise RuntimeError(f"insufficient guard rows around u={u0}")
    completed0, stream0 = _recomplete_and_emit(loc)
    u = loc.u.to_numpy(float)
    i = int(np.argmin(np.abs(u - u0)))
    shape = _local_bump(u, u0, width)
    if shape[i] < 0.99:
        raise RuntimeError("local bump did not peak at audit point")

    # Verify local absolute re-emission against the accepted global build.
    ug = full_stream.u.to_numpy(float)
    jg = int(np.argmin(np.abs(ug - u[i])))
    replay_err = max(
        abs(float(stream0[s].iloc[i]) - float(full_stream[s].iloc[jg]))
        / max(1.0, abs(float(full_stream[s].iloc[jg])))
        for s in SLOT_NAMES
    )

    E0 = _eom_vector(completed0, i)
    base_K = {}
    base_lam = {}
    base_vec = {}
    for L in LS:
        K = _sym(kinetic_schur(stream0, L)["K"])
        w, v = np.linalg.eigh(K[i])
        base_K[L] = K
        base_lam[L] = float(w[0])
        base_vec[L] = [float(x) for x in v[:, 0]]

    nd = len(DIRECTIONS)
    B = np.zeros((4, nd))
    gradients = {L: np.zeros(nd) for L in LS}
    a1_grad = np.zeros(nd)

    for k, name in enumerate(DIRECTIONS):
        if name in SVT_DIRECTIONS:
            # One-sided tangent derivative.  The action/background maps are
            # algebraic-linear in these local profile amplitudes; a dedicated
            # finite-difference convergence check is stored separately at the
            # light ring before this audit may be promoted.
            d = loc.copy()
            d[name] = d[name].to_numpy(float) + FD_EPS * shape
            comp, st = _recomplete_and_emit(d)
            B[:, k] = (_eom_vector(comp, i) - E0) / FD_EPS
            a1_grad[k] = (float(st.a1.iloc[i]) - float(stream0.a1.iloc[i])) / FD_EPS
            for L in LS:
                K = _sym(kinetic_schur(st, L)["K"])
                gradients[L][k] = (float(np.linalg.eigvalsh(K[i])[0]) - base_lam[L]) / FD_EPS
        elif name == "H_G2XX":
            # Exact background-null Horndeski G2XX response changes only
            # c2,c6,e2, none of which enter the locked kinetic Schur map.
            B[:, k] = 0.0
            a1_grad[k] = 0.0
            for L in LS:
                gradients[L][k] = 0.0
        elif name == "H_G3XX_TUBULAR":
            # This exact tubular response changes only c2,c3,e3.  The locked
            # kinetic Schur map contains none of those slots, so its K response
            # is identically zero; no expensive reducer call is needed.
            B[:, k] = 0.0
            a1_grad[k] = 0.0
            for L in LS:
                gradients[L][k] = 0.0
        else:
            raise AssertionError(name)

    metric_B = B[:3]
    full_B = B
    perL = {}
    for L in LS:
        p3 = _ranks_and_projection(metric_B, gradients[L])
        p4 = _ranks_and_projection(full_B, gradients[L])
        # JSON cannot carry numpy arrays.
        for p in (p3, p4):
            arr = p.pop("native_projected_gradient")
            inds = np.argsort(-np.abs(arr))[:6]
            p["dominant_native_projected_components"] = [
                {"direction": DIRECTIONS[int(j)], "component": float(arr[j])} for j in inds
            ]
        perL[str(L)] = {
            "lambda_min": base_lam[L],
            "weakest_vector": base_vec[L],
            "metric_vector_tangent": p3,
            "with_scalar_identity_tangent": p4,
        }

    lp3 = _common_gain_lp(metric_B, gradients)
    lp4 = _common_gain_lp(full_B, gradients)
    for lp in (lp3, lp4):
        if lp.get("success"):
            v = lp.pop("direction")
            inds = np.argsort(-np.abs(v))[:8]
            lp["dominant_native_components"] = [
                {"direction": DIRECTIONS[int(j)], "component": float(v[j])} for j in inds
            ]
            lp["induced_a1_first_order"] = float(a1_grad @ v)

    rows = []
    for k, name in enumerate(DIRECTIONS):
        row = {
            "u_requested": float(u0),
            "u_grid": float(u[i]),
            "direction": name,
            "dE00": float(B[0, k]),
            "dE11": float(B[1, k]),
            "dJA": float(B[2, k]),
            "dEphi_identity": float(B[3, k]),
            "da1": float(a1_grad[k]),
        }
        for L in LS:
            row[f"dlambda_L{L}"] = float(gradients[L][k])
        rows.append(row)

    return {
        "u_requested": float(u0),
        "u_grid": float(u[i]),
        "bump_width_u": float(width),
        "rows_local_grid": int(len(loc)),
        "absolute_reemit_max_scaled_vs_global": float(replay_err),
        "baseline_background": {
            "E00": float(E0[0]),
            "E11": float(E0[1]),
            "JA": float(E0[2]),
            "Ephi_identity": float(E0[3]),
        },
        "B_rank_metric_vector": _ranks_and_projection(metric_B, gradients[LS[0]])["rank_B"],
        "B_rank_with_scalar_identity": _ranks_and_projection(full_B, gradients[LS[0]])["rank_B"],
        "per_L": perL,
        "simultaneous_all_L": {
            "metric_vector_tangent": lp3,
            "with_scalar_identity_tangent": lp4,
        },
    }, rows


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    build = build_onshell_central(ROOT)
    action = build.action.sort_values("x").reset_index(drop=True)
    stream = build.direct41.sort_values("x").reset_index(drop=True)

    points = []
    table_rows = []
    for u0 in SAMPLE_POINTS:
        p, rows = _point_audit(action, stream, float(u0))
        points.append(p)
        table_rows.extend(rows)
        print(
            f"u={p['u_grid']:.9f} rank3={p['B_rank_metric_vector']} "
            f"rank4={p['B_rank_with_scalar_identity']} "
            f"common3={p['simultaneous_all_L']['metric_vector_tangent'].get('max_common_first_order_gain_native_box')} "
            f"common4={p['simultaneous_all_L']['with_scalar_identity_tangent'].get('max_common_first_order_gain_native_box')}",
            flush=True,
        )

    # The current absolute Hybrid action ends just below .71.  Do not fabricate
    # an absolute H+SVT member at .712 from the rejected historical Inner stream.
    right = {
        "u_requested": REQUESTED_RIGHT_DIAGNOSTIC,
        "status": "NOT_EVALUATED_NO_CERTIFIABLE_ABSOLUTE_ACTION_SEED",
        "reason": (
            "The accepted Electric-Hybrid absolute action table ends below u=0.71. "
            "The old Inner background/action data remain useful seeds but the full absolute "
            "H+SVT 41-emitter right of the historical seam is not yet certified.  Using the "
            "historical Inner 41 stream here would violate the action-first contract."
        ),
        "next_requirement": "extend the joint action continuation across u=0.71, then rerun this identical tangent audit",
    }

    # Structural conclusion: require simultaneous positive first-order gain at
    # every sampled point under the strict metric+vector tangent.  The scalar-
    # identity variant is reported separately because its baseline gate is not
    # yet globally promoted.
    controllable3 = all(
        p["simultaneous_all_L"]["metric_vector_tangent"].get("all_required_L_gain_positive", False)
        for p in points
    )
    controllable4 = all(
        p["simultaneous_all_L"]["with_scalar_identity_tangent"].get("all_required_L_gain_positive", False)
        for p in points
    )
    all_structural4 = all(
        all(v["with_scalar_identity_tangent"]["structurally_controllable"] for v in p["per_L"].values())
        for p in points
    )

    report = {
        "status": "PASS_LOCAL_TANGENT_CONTROLLABILITY" if controllable3 else "OPEN_LOCAL_TANGENT_CONTROLLABILITY",
        "scope": "pre-continuation local action tangent/rank audit; not a finite transition member",
        "action_member": "ELECTRIC_HYBRID_ONSHELL_CENTRAL diagnostic extension on u<0.71",
        "sample_points": points,
        "requested_right_seed_diagnostic": right,
        "directions": list(DIRECTIONS),
        "finite_difference_epsilon": FD_EPS,
        "required_L": list(LS),
        "global_summary": {
            "metric_vector_tangent_has_one_direction_improving_all_required_L_at_all_sampled_points": controllable3,
            "extended_scalar_identity_tangent_has_one_direction_improving_all_required_L_at_all_sampled_points": controllable4,
            "all_per_L_gradients_have_nonzero_projection_in_extended_tangent": all_structural4,
            "interpretation": (
                "A positive result is a local first-order controllability statement only. "
                "It licenses a predictor-corrector continuation attempt; it does not prove a finite healthy branch."
            ),
        },
        "promotion": {
            "finite_strong_field_transition": False,
            "absolute_direct41_transition": False,
            "global_KRGSM": False,
            "same_operator_QNM": False,
        },
    }
    REPORT.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    pd.DataFrame(table_rows).to_csv(TABLE, index=False)
    print(json.dumps(report["global_summary"], indent=2))
    print(f"WROTE {REPORT}")
    print(f"WROTE {TABLE}")
    return 0 if controllable3 else 2


if __name__ == "__main__":
    raise SystemExit(main())
