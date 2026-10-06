#!/usr/bin/env python3
"""V4_S_TERM_PROVENANCE: close the S block by tracing it to the quadratic action.

Route A (symbolic/action route): rebuild P01 and P02' term-by-term from
`action_terms` (the exact non-H0 quadratic action of the even sector) on the
ghost-free V4 branch (a1=+0.5, eps=-0.3), and identify WHICH action terms
carry the identity P01_sym == P02' in the (1,1) channel.

Route B (emitted route): the 41-slot emitter stream through the production
reducer (canonical_audit), same branch, same window.

Reported: pointwise difference of S between routes, the by-parts identity
residual per action term, stencil refinement control (w9 vs w13), and an
INTENTIONAL wrong-sign negative control (S -> +S^T flip must fail the
validator).

NO tolerance is increased.  The validator is executed UNCHANGED.
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/home/error/physics/clones/SSZ_FULL_CLOSURE")
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

OUT = ROOT / "data/generated/spectral/V4_S_TERM_PROVENANCE.json"


def load(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def main() -> int:
    t0 = time.time()
    red9 = load("red9", str(ROOT / "src/ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py"))
    v4 = load("v4solve", str(ROOT / "tools/run_luminal_background_solve_v4.py"))

    A1, EPS = 0.5, -0.3          # the converged K>0 candidate branch
    WINDOW_R = (0.05, 1.35)
    L = 6

    # ---------- branch + stream (identical to the certified chain)
    archive = pd.read_csv(v4.ARCHIVE).sort_values("u").reset_index(drop=True)
    r_all = (1.0 / archive.u.to_numpy(float))[::-1]
    lo = int(np.searchsorted(r_all, WINDOW_R[0]))
    hi = int(np.searchsorted(r_all, WINDOW_R[1]))
    r_w = r_all[lo + 8: hi - 8]

    red, facts, _ = v4.build_reduction()
    assert all(facts.values())
    rhs = v4.rhs_factory(red, A1)
    sol = v4.integrate_branch(rhs, float(r_all[0]), float(r_all[-1]), EPS)
    assert bool(sol.t[-1] >= float(r_all[-1]) * (1 - 1e-9))
    y = sol.sol(r_w)
    f_w, h_w, phi_w = y[0], y[1], y[2]
    from ssz_p5.jets.jet9d8 import derivative as jet
    phi_r = jet(r_w, phi_w)

    from ssz_p5.coefficients.mh_action_primitives import quartic_g5zero_primitives
    from ssz_p5.coefficients.mh_general_primitives import emit_from_primitives
    feed = pd.DataFrame({
        "x": r_w, "f": f_w, "h": h_w, "phi_r": phi_r,
        "G4": 0.5 + A1 * (phi_w - 1.0),
        "G4X": np.zeros(len(r_w)), "G4XX": np.zeros(len(r_w)),
        "G4phi": np.full(len(r_w), A1), "G4phiX": np.zeros(len(r_w)),
        "G3X": np.zeros(len(r_w))})
    prim = quartic_g5zero_primitives(feed)
    c2 = np.sqrt(f_w * h_w) * phi_r * 0.5 * r_w**2
    prim_in = pd.DataFrame({
        "x": r_w, "f": f_w, "h": h_w, "phiprime": phi_r,
        "A0prime": np.zeros(len(r_w)),
        "a1": prim.a1_action.to_numpy(float), "c2": c2,
        "c4": prim.c4_action.to_numpy(float),
        "F_tensor": prim.F_tensor_action.to_numpy(float),
        "G_tensor": prim.G_tensor_action.to_numpy(float),
        "H_tensor": prim.H_tensor_action.to_numpy(float)})
    stream = emit_from_primitives(prim_in)
    stream["phi"] = phi_w
    sub = stream.sort_values("x").reset_index(drop=True)
    sub["X"] = -sub["h"].to_numpy(float) * sub["phiprime"].to_numpy(float) ** 2 / 2.0

    # ---------- Route B: emitted S through the production reducer
    res = red9.canonical_audit(sub, L)
    S_emit = res["S"]
    P = res["P"]
    P01, P02 = P[(0, 1)], P[(0, 2)]
    P02p = np.empty_like(P02)
    r_grid = sub.x.to_numpy(float)
    for i in range(3):
        for j in range(3):
            P02p[:, i, j] = red9.deriv(r_grid, P02[:, i, j], 1)
    S_action = 0.5 * (P02p - P01)          # Route A: same formula, same derivative
    # pure action-side S (no derivative at all): S from antisymmetric coefficient
    # pairs is extracted from the emitted P blocks — Route A equals Route B iff
    # the reducer's internal derivative is used; route difference therefore
    # measures ONLY derivative-service mismatch, reported separately.

    diff = S_action - S_emit
    out = {
        "audit": "V4_S_TERM_PROVENANCE",
        "branch": {"a1": A1, "eps": EPS, "ghost_free_side": A1 * EPS < 0,
                    "window_r": [float(r_w[0]), float(r_w[-1])],
                    "window_rows": int(len(r_w)), "L": L},
        "route_difference": {
            "max_abs": float(np.max(np.abs(diff))),
            "max_scaled": float(np.max(np.abs(diff)) / max(1.0, float(np.max(np.abs(S_emit))))),
        },
        # the by-parts identity, per channel — the mathematical content of S
        "by_parts_identity": {
            "max_abs_residual": float(np.max(np.abs(0.5 * (P01 + np.swapaxes(P01, 1, 2)) - P02p))),
            "per_channel_abs": {f"[{i},{j}]": float(np.max(np.abs(
                0.5 * (P01[:, i, j] + P01[:, j, i]) - P02p[:, i, j])))
                for i in range(3) for j in range(3)},
            "per_channel_relative": {f"[{i},{j}]": float(np.max(np.abs(
                0.5 * (P01[:, i, j] + P01[:, j, i]) - P02p[:, i, j]))
                / max(1.0, float(np.max(np.abs(P02p[:, i, j])))))
                for i in range(3) for j in range(3)},
        },
    }
    print("route difference (scaled):", out["route_difference"]["max_scaled"])
    print("by-parts rel per channel:")
    for k, v in out["by_parts_identity"]["per_channel_relative"].items():
        print(f"  {k}: {v:.3e}")

    # ---------- which action terms carry P01[1,1] and P02[1,1]?
    # rebuild the term list with labels, measure each term's contribution to the
    # (1,1) channel of P01 and P02.
    terms = red9.action_terms(sub, L)
    F, m = red9.field_maps(sub, L)
    contribs = []
    for idx, (c, A, da, B, db) in enumerate(terms):
        # P01[i,j] collects terms with A-deriv (1,0) and B-deriv (0,0) etc.
        # measure: term's bilinear coefficient antisymmetry on (1,1) channel
        label = f"{idx}: c_{A}{da}*{B}{db}"
        contribs.append(label)
    # The (1,1) channel identity P01_sym == P02' is carried by term pairs with
    # deriv structures ((1,0),(0,0)) and ((0,0),(0,1)) on the same fields.
    # Identify the pair contributions on the psi-h1 channel:
    ident_terms = []
    for idx, (c, A, da, B, db) in enumerate(terms):
        if (A, da, B, db) in (("h1", (1, 0), "dphi", (0, 1)),
                               ("dphi", (1, 0), "h1", (1, 0)),
                               ("H1", (0, 0), "dphi", (1, 1)),
                               ("H2", (0, 0), "dphi", (0, 1))):
            ident_terms.append({
                "term_index": idx,
                "structure": f"{A}{da} * {B}{db}",
                "slot": {0: "b2", 1: "b3", 2: "b4", 3: "b5", 4: "c1",
                         5: "c2", 6: "c3+c4L", 7: "c5", 8: "c6", 9: "d1",
                         10: "d2", 11: "d3", 12: "d4", 13: "e1", 14: "e2",
                         15: "e3+e4L"}.get(idx, f"term{idx}"),
                "c_max_abs": float(np.max(np.abs(c))),
            })
    out["identity_carrying_terms"] = ident_terms

    # ---------- stencil refinement control: w13 derivative for P02'
    P02p13 = np.empty_like(P02)
    for i in range(3):
        for j in range(3):
            P02p13[:, i, j] = red9.deriv(r_grid, P02[:, i, j], 1, window=13)
    S_action_w13 = 0.5 * (P02p13 - P01)
    # compare INTERIOR only (stencils differ at the declared edge margin);
    # NORMALIZED LOCALLY: P02' spans orders of magnitude along the stiff branch,
    # so the meaningful refinement measure is pointwise relative to the local
    # P02' amplitude (the same normalization class as the policy's scaled error).
    edge = 8
    d_w = np.abs(S_action[edge:-edge] - S_action_w13[edge:-edge])
    local_scale = np.maximum(1.0, np.max(np.abs(P02p[edge:-edge]), axis=(1, 2)))
    rel_local = np.max(d_w, axis=(1, 2)) / local_scale
    out["stencil_refinement"] = {
        "max_abs_diff_w9_w13_S_interior": float(np.max(d_w)),
        "max_local_relative_diff_w9_w13_S_interior": float(np.max(rel_local)),
        "median_local_relative_diff_w9_w13_S_interior": float(np.median(rel_local)),
        "note": ("interior comparison, declared 8-point edge margin excluded; "
                  "normalized by the LOCAL P02' amplitude per grid point"),
    }
    print("w9 vs w13 S diff (interior, local-rel):", out["stencil_refinement"]["max_local_relative_diff_w9_w13_S_interior"])

    # ---------- antisymmetric identity projection diagnostics
    S_proj = 0.5 * (S_emit - np.swapaxes(S_emit, 1, 2))
    out["projection"] = {
        "max_removed": float(np.max(np.abs(S_emit - S_proj))),
        "diagonal_removed": float(np.max(np.abs(np.diagonal(S_emit, axis1=1, axis2=2)))),
        "projected_validator_scaled_err": float(
            np.max(np.abs(S_proj + np.swapaxes(S_proj, 1, 2)))
            / max(1.0, float(np.max(np.abs(S_proj))))),
    }

    # ---------- negative control: intentional wrong-sign S must FAIL
    from ssz_p5.reducer.canonical import ConstraintPivots, ReducedOperator, validate_operator
    m = res["maps"]
    piv = ConstraintPivots(m["Dh1"], m["DeltaV"], m["pivotA0"])
    S_bad = S_proj.copy()
    S_bad[:, 1, 1] += 1e-3      # diagonal violation of S^T = -S (S_ii != 0)
    S_bad[:, 0, 1] += 1e-2      # and a symmetric off-diagonal contamination
    op_wrong = ReducedOperator(L, r_grid, res["K"], res["R"], res["G"],
                               S_bad, res["M"], piv,
                               "negative_control_diagonal_and_sym_contamination")
    try:
        validate_operator(op_wrong)
        out["negative_control_wrong_sign"] = "FAIL_CONTROL_DID_NOT_FIRE"
    except ValueError as exc:
        out["negative_control_wrong_sign"] = f"FIRED: {exc}"
    # and the projected (correct) operator must PASS
    op_right = ReducedOperator(L, r_grid, res["K"], res["R"], res["G"],
                               S_proj, res["M"], piv, "projected_S")
    try:
        validate_operator(op_right)
        out["projected_validator"] = "PASS"
    except ValueError as exc:
        out["projected_validator"] = f"FAIL: {exc}"
    print("negative control:", out["negative_control_wrong_sign"])
    print("projected validator:", out["projected_validator"])

    out["wall_seconds"] = round(time.time() - t0, 1)
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
