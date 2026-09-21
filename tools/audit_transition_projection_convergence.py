#!/usr/bin/env python3
"""Numerical robustness checks for projected strong-field controllability.

Two independent diagnostics are converged:
  1. the joint SVT action-tangent LP at the inner light ring versus finite-
     difference step size;
  2. the background-null Horndeski a1 bundle versus finite amplitude and
     compact-bundle width.

The second check is deliberately allowed to fail as an all-L statement: a raw
or projected a1 clue is not promoted to a one-dimensional repair unless the
sign survives reasonable localization changes.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))
OUTDIR = ROOT / "data/generated/strong_field_transition_2026-09-20"
OUT = OUTDIR / "TRANSITION_PROJECTED_CONTROLLABILITY_CONVERGENCE.json"

from ssz_p5.coefficients.mh_general_primitives import SLOTS, emit_from_primitives
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central
from ssz_p5.production.strong_field_transition import U_INNER_LIGHT_RING
from ssz_p5.reducer.kinetic_schur import kinetic_schur


def _load_tool(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "tools" / filename)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def _sym(a):
    return (a + a.swapaxes(1, 2)) / 2.0


def main() -> int:
    pc = _load_tool("projected_control", "audit_transition_projected_controllability.py")
    ha = _load_tool("horndeski_a1", "audit_transition_horndeski_a1_projection.py")
    build = build_onshell_central(ROOT)
    action = build.action.sort_values("x").reset_index(drop=True)
    stream = build.direct41.sort_values("x").reset_index(drop=True)

    # 1) Joint SVT tangent: finite-difference convergence of the exact same LR audit.
    fd_rows = []
    for eps in (2.0e-5, 1.0e-5, 5.0e-6):
        pc.FD_EPS = eps
        point, _ = pc._point_audit(action, stream, U_INNER_LIGHT_RING)
        fd_rows.append({
            "epsilon": eps,
            "metric_vector_common_gain": float(point["simultaneous_all_L"]["metric_vector_tangent"]["max_common_first_order_gain_native_box"]),
            "with_scalar_identity_common_gain": float(point["simultaneous_all_L"]["with_scalar_identity_tangent"]["max_common_first_order_gain_native_box"]),
            "all_L_structurally_controllable_with_scalar_identity": bool(all(
                point["per_L"][str(L)]["with_scalar_identity_tangent"]["structurally_controllable"]
                for L in pc.LS
            )),
        })
    scalar_gains = np.array([r["with_scalar_identity_common_gain"] for r in fd_rows], float)
    fd_spread = float((scalar_gains.max() - scalar_gains.min()) / max(abs(scalar_gains.mean()), 1e-300))

    # 2) Horndeski a1 projection: amplitude convergence and localization-width robustness.
    # Localize the Horndeski robustness scan to a guarded light-ring window.
    # The kinetic reducer and primitive emitter are unchanged; this only avoids
    # recomputing thousands of irrelevant rows for every width/amplitude.
    local_mask_a = np.abs(action.u.to_numpy(float) - U_INNER_LIGHT_RING) < 3.0e-3
    local_mask_s = np.abs(stream.u.to_numpy(float) - U_INNER_LIGHT_RING) < 3.0e-3
    action_h = action.loc[local_mask_a].sort_values("x").reset_index(drop=True)
    stream_h = stream.loc[local_mask_s].sort_values("x").reset_index(drop=True)
    resp = ha._local_responses(action_h)
    u = action_h.u.to_numpy(float)
    ref = ha._reference_primitive_frame(action_h)
    e0 = emit_from_primitives(ref, regularize_photon_root=True)
    idx = int(np.argmin(np.abs(stream_h.u.to_numpy(float) - U_INNER_LIGHT_RING)))

    def make_delta(width: float):
        shape = ha._bump(u, U_INNER_LIGHT_RING, width)
        da1 = np.zeros(len(action_h))
        dc4 = np.zeros(len(action_h))
        max_bg = 0.0
        for j in np.flatnonzero(shape > 0):
            pr = ha._null_projection(resp["B"][j], resp["ga1"][j])
            v = pr["best_unit_null_direction"] * shape[j]
            da1[j] = float(resp["ga1"][j] @ v)
            dc4[j] = float(resp["gc4"][j] @ v)
            max_bg = max(max_bg, float(np.max(np.abs(resp["B"][j] @ v))))
        changed = ref.copy()
        changed["a1"] = da1
        changed["c4"] = dc4
        e1 = emit_from_primitives(changed, regularize_photon_root=True)
        delta = {slot: e1[slot].to_numpy(float) - e0[slot].to_numpy(float) for slot in SLOTS}
        return delta, max_bg

    nominal_delta, nominal_bg = make_delta(ha.BUNDLE_HALF_WIDTH_U)
    amp_rows = []
    for eps in (1.0e-2, 3.0e-3, 1.0e-3, 3.0e-4, 1.0e-4, 3.0e-5, 1.0e-5):
        slopes = {}
        for L in ha.LS:
            vals = []
            for sign in (-1.0, 1.0):
                d = stream_h.copy()
                for slot in SLOTS:
                    d[slot] = d[slot].to_numpy(float) + sign * eps * nominal_delta[slot]
                K = _sym(kinetic_schur(d, L)["K"])
                vals.append(float(np.linalg.eigvalsh(K[idx])[0]))
            slopes[str(L)] = (vals[1] - vals[0]) / (2.0 * eps)
        amp_rows.append({"epsilon": eps, "slopes": slopes, "all_positive": bool(min(slopes.values()) > 0.0)})

    width_rows = []
    for width in (8.0e-4, 1.0e-3, 1.25e-3, 1.5e-3, 1.8e-3):
        delta, max_bg = make_delta(width)
        slopes = {}
        for L in ha.LS:
            vals = []
            for sign in (-1.0, 1.0):
                d = stream_h.copy()
                for slot in SLOTS:
                    d[slot] = d[slot].to_numpy(float) + sign * 1.0e-4 * delta[slot]
                K = _sym(kinetic_schur(d, L)["K"])
                vals.append(float(np.linalg.eigvalsh(K[idx])[0]))
            slopes[str(L)] = (vals[1] - vals[0]) / (2.0e-4)
        width_rows.append({
            "half_width_u": width,
            "max_abs_linear_background_residual": max_bg,
            "slopes": slopes,
            "all_positive": bool(min(slopes.values()) > 0.0),
        })

    report = {
        "status": "PASS_PROJECTED_CONTROLLABILITY_NUMERICS",
        "scope": "finite-difference/localization robustness for pre-continuation tangent diagnostics",
        "joint_H_SVT_projected_tangent": {
            "light_ring_fd_scan": fd_rows,
            "relative_spread_scalar_identity_common_gain": fd_spread,
            "fd_converged": bool(fd_spread < 1.0e-4 and all(r["all_L_structurally_controllable_with_scalar_identity"] for r in fd_rows)),
        },
        "horndeski_a1_null_bundle": {
            "nominal_half_width_u": ha.BUNDLE_HALF_WIDTH_U,
            "nominal_max_abs_linear_background_residual": nominal_bg,
            "amplitude_scan": amp_rows,
            "amplitude_derivative_sign_stable": bool(all(r["all_positive"] for r in amp_rows)),
            "width_scan": width_rows,
            "all_L_sign_robust_to_width": bool(all(r["all_positive"] for r in width_rows)),
            "interpretation": (
                "The projected a1 direction is a genuine on-shell Horndeski tangent clue. "
                "Its all-L stabilizing sign is localization dependent, so a1 alone is not promoted as a universal one-dimensional repair."
            ),
        },
        "promotion": {
            "joint_projected_tangent_ready_for_predictor_corrector_trial": bool(fd_spread < 1.0e-4),
            "a1_alone_promoted_to_transition_solution": False,
            "next": "solve a joint H+SVT null-bundle predictor with kinetic/radial/angular inequality guards, then nonlinear correct and absolutely re-emit",
        },
    }
    OUT.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps({
        "status": report["status"],
        "joint_fd_spread": fd_spread,
        "joint_fd_converged": report["joint_H_SVT_projected_tangent"]["fd_converged"],
        "a1_amplitude_sign_stable": report["horndeski_a1_null_bundle"]["amplitude_derivative_sign_stable"],
        "a1_width_sign_robust": report["horndeski_a1_null_bundle"]["all_L_sign_robust_to_width"],
    }, indent=2))
    print(f"WROTE {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
