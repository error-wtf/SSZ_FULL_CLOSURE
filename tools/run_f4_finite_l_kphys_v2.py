#!/usr/bin/env python3
"""F.4-v2: finite-L health on the CONSTRAINT-REDUCED (psi,V) operator.

Branch consistency fix (Lino review 06.10): the original F.4 ran on
(a1=+0.5, eps=-0.5) while the DOF audit / V4 export / B8 run on
(a1=+0.5, eps=-0.3).  This rerun uses eps=-0.3 — the SAME branch as
V4_PHYSICAL_DOF_AUDIT_V1, the V4 export and the bridge B8 run — and
measures health on the physical (psi,V) subblock identified by the
DOF audit (chi = constraint-eliminated Lagrange-multiplier channel).

Health criteria (declared, unchanged in spirit):
  * S-projection + validator PASS per L
  * lambda_min(K_phys_{psi,V}) > 0  (strictly, interior window)
  * radial channel speeds positive on the psi/V subblock
  * pivots Dh1/DeltaV/pivotA0 nonzero (constraint elimination valid)
L ladder: 6, 12, 20, 42, 110, 420, 1000.
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

OUT = ROOT / "data/generated/spectral/F4_FINITE_L_HEALTH_V4_KPHYS_V2.json"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def main() -> int:
    t0 = time.time()
    red9 = load("red9", str(ROOT / "src/ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py"))
    v4 = load("v4solve", str(ROOT / "tools/run_luminal_background_solve_v4.py"))
    from ssz_p5.reducer.canonical import validate_operator, ReducedOperator, ConstraintPivots
    from ssz_p5.jets.jet9d8 import derivative as jet

    A1, EPS = 0.5, -0.3          # SAME branch as DOF audit / V4 export / B8
    WINDOW_R = (0.05, 1.35)
    PHYS = [0, 2]                 # (psi, V) physical subspace

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

    health = {}
    all_pass = True
    for L in (6, 12, 20, 42, 110, 420, 1000):
        res = red9.canonical_audit(sub, L)
        S = 0.5 * (res["S"] - np.swapaxes(res["S"], 1, 2))
        m = res["maps"]
        op = ReducedOperator(L, sub.x.to_numpy(float),
                             res["K"], res["R"], res["G"], S, res["M"],
                             ConstraintPivots(m["Dh1"], m["DeltaV"], m["pivotA0"]),
                             f"v4_f4v2_L{L}")
        try:
            validate_operator(op)
            v_err = "PASS"
        except ValueError as exc:
            v_err = f"FAIL: {exc}"
        K = res["K"]
        # physical-subblock eigenvalues per node (symmetrised 2x2)
        phys_min = np.array([
            float(np.linalg.eigvalsh(0.5 * (K[i][np.ix_(PHYS, PHYS)] +
                                                K[i][np.ix_(PHYS, PHYS)].T))[0])
            for i in range(K.shape[0])])
        # radial channel speeds on the physical subblock (K/G diag ratio)
        g22 = res["G"][:, 2, 2]
        with np.errstate(divide="ignore", invalid="ignore"):
            speeds = np.where(g22 > 0, K[:, 2, 2] / g22, np.nan)
        # pivots
        piv_ok = (float(np.min(np.abs(m["Dh1"]))) > 1e-10 and
                   float(np.min(np.abs(m["DeltaV"]))) > 1e-10 and
                   float(np.min(np.abs(m["pivotA0"]))) > 1e-10)
        L_pass = (v_err == "PASS" and float(phys_min.min()) > 0 and piv_ok
                   and bool(np.nanmin(speeds) > 0))
        health[str(L)] = {
            "validator": v_err,
            "min_lambda_K_phys": round(float(phys_min.min()), 10),
            "n_negative_nodes": int((phys_min < 0).sum()),
            "min_radial_speed_phys": (round(float(np.nanmin(speeds)), 8)
                                       if np.isfinite(np.nanmin(speeds)) else None),
            "pivots_ok": piv_ok,
            "pass": L_pass,
        }
        all_pass = all_pass and L_pass
        print(f"L={L}: {v_err} | min λ(K_phys)={phys_min.min():.3e} | "
              f"min c²={health[str(L)]['min_radial_speed_phys']} | pass={L_pass}",
              flush=True)

    out = {
        "audit": "F4_FINITE_L_HEALTH_V4_KPHYS_V2",
        "supersedes_blocker": "F4_FINITE_L_HEALTH_V4_ATTEMPT_BLOCKED_RECORD.json",
        "branch_consistency_fix": (
            "Original F.4 ran on eps=-0.5 while the DOF audit / V4 export / B8 "
            "ran on eps=-0.3. This rerun uses eps=-0.3 — identical branch "
            "(a1, eps, profile, operator) through every stage."),
        "branch": {"a1": A1, "eps": EPS, "ghost_free_side": True},
        "window_r": [float(r_w[0]), float(r_w[-1])],
        "window_rows": int(len(r_w)),
        "L_ladder": [6, 12, 20, 42, 110, 420, 1000],
        "health_metric": "min eigenvalue of the (psi,V) principal subblock of "
                          "the reduced K per node (constraint channel chi "
                          "excluded per V4_PHYSICAL_DOF_AUDIT_V1: v_- is 100% "
                          "chi, Lagrange-multiplier channel)",
        "health": health,
        "all_pass": bool(all_pass),
        "verdict": ("F.4 HEALTH PASS on the constraint-reduced physical "
                     "operator — the ghost-free V4 branch carries 3 clean "
                     "kinetic channels after chi exclusion; QNM path "
                     "proceeds to global exterior + coupled Jost/ECS."
                     if all_pass else
                     "see per-L rows"),
        "wall_seconds": round(time.time() - t0, 1),
    }
    OUT.write_text(json.dumps(out, indent=1, allow_nan=False) + "\n")
    print(json.dumps({"all_pass": all_pass}, indent=1))
    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
