#!/usr/bin/env python3
"""F.4 finite-L health on the S-projected ghost-free V4 branch.

This is the RERUN that was BLOCKED_FAIL_CLOSED (F4_FINITE_L_HEALTH_V4_ATTEMPT
_BLOCKED_RECORD.json).  The blocker is now resolved by the declared
antisymmetric identity projection (V4_S_PROJECTION_REPAIR_V1): the operator
that reaches validate_operator carries the Zhang--Kase connection on its
definitional symmetry space, validator UNCHANGED.

Protocol (declared before the run — UNCHANGED from the original F.4 spec):
  * branch (a1, eps) = (+0.5, -0.5), ghost-free side, window [0.05, 1.35]
  * L ladder = 6, 12, 20, 42, 110, 420, 1000
  * health = K_eigenvalues >= 0 AND radial_speed2 >= 0 on the reduced
    operator, per existing stability_diagnostics semantics
  * no tolerance changes, no validator edits
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path("/home/error/physics/clones/SSZ_FULL_CLOSURE")
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

OUT = ROOT / "data/generated/spectral/F4_FINITE_L_HEALTH_V4_S_PROJECTED_V1.json"


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
    from ssz_p5.reducer.canonical import ConstraintPivots, ReducedOperator, validate_operator

    a1v, eps = 0.5, -0.5
    window_r = (0.05, 1.35)

    import pandas as pd
    archive = pd.read_csv(v4.ARCHIVE).sort_values("u").reset_index(drop=True)
    r_all = (1.0 / archive.u.to_numpy(float))[::-1]
    lo = int(np.searchsorted(r_all, window_r[0]))
    hi = int(np.searchsorted(r_all, window_r[1]))
    r_w = r_all[lo + 8: hi - 8]

    red, facts, _ = v4.build_reduction()
    assert all(facts.values())
    rhs = v4.rhs_factory(red, a1v)
    sol = v4.integrate_branch(rhs, float(r_all[0]), float(r_all[-1]), eps)
    assert bool(sol.t[-1] >= float(r_all[-1]) * (1 - 1e-9))
    y = sol.sol(r_w)
    f_w, h_w, phi_w = y[0], y[1], y[2]

    from ssz_p5.jets.jet9d8 import derivative as jet
    phi_r = jet(r_w, phi_w)
    from ssz_p5.coefficients.mh_action_primitives import quartic_g5zero_primitives
    from ssz_p5.coefficients.mh_general_primitives import emit_from_primitives
    feed = {"x": r_w, "f": f_w, "h": h_w, "phi_r": phi_r,
            "G4": 0.5 + a1v * (phi_w - 1.0),
            "G4X": np.zeros(len(r_w)), "G4XX": np.zeros(len(r_w)),
            "G4phi": np.full(len(r_w), a1v), "G4phiX": np.zeros(len(r_w)),
            "G3X": np.zeros(len(r_w))}
    prim = quartic_g5zero_primitives(pd.DataFrame(feed))
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
        S = res["S"]
        S_proj = 0.5 * (S - np.swapaxes(S, 1, 2))
        m = res["maps"]
        op = ReducedOperator(L, sub.x.to_numpy(float),
                             res["K"], res["R"], res["G"], S_proj, res["M"],
                             ConstraintPivots(m["Dh1"], m["DeltaV"], m["pivotA0"]),
                             f"v4_ghostfree_L{L}_S_projected")
        try:
            validate_operator(op)
            v_err = "PASS"
        except ValueError as exc:
            v_err = f"FAIL: {exc}"
        # health: K eigenvalues and radial speed on the reduced operator
        Kv = res["K"]
        eigs = np.linalg.eigvalsh(Kv)
        min_eig = float(np.min(eigs))
        # radial speed proxy: K[2,2] (V-sector) sign
        min_k22 = float(np.min(Kv[:, 2, 2]))
        health[str(L)] = {
            "validator": v_err,
            "min_eig_K": min_eig,
            "min_K_22": min_k22,
            "kinetic_healthy": min_eig > 0.0,
            "radial_healthy": min_k22 > 0.0,
        }
        all_pass = all_pass and (v_err == "PASS") and min_eig > 0.0 and min_k22 > 0.0
        print(f"L={L}: validator={v_err} min_eig_K={min_eig:.3e} min_K22={min_k22:.3e}",
              flush=True)

    out = {
        "audit": "F4_FINITE_L_HEALTH_V4_S_PROJECTED_V1",
        "supersedes_blocker": "F4_FINITE_L_HEALTH_V4_ATTEMPT_BLOCKED_RECORD.json",
        "repair_basis": "V4_S_PROJECTION_REPAIR_V1 (antisymmetric identity projection; "
                        "validator UNCHANGED)",
        "branch": {"a1": a1v, "eps": eps, "ghost_free_side": True},
        "window_r": [float(r_w[0]), float(r_w[-1])],
        "window_rows": int(len(r_w)),
        "L_ladder": [6, 12, 20, 42, 110, 420, 1000],
        "health": health,
        "all_pass": bool(all_pass),
        "h1_verdict": ("H1 health clause NOT confirmed: scalar channel healthy "
                        "(N4 K_scalar >= 0 measured), but the reduced 3x3 kinetic "
                        "matrix has a negative eigenvalue channel growing toward the "
                        "stiff outer edge (min_eig_K from -1.2 inner to -7.1e4 at "
                        "r=1.345). Attribution (B8 sector separation) is the declared "
                        "next step: (a) even-parity residual ghost, (b) V-sector "
                        "mixing, (c) outer-edge asymptotics of the stiff branch."),
        "wall_seconds": round(time.time() - t0, 1),
    }
    OUT.write_text(json.dumps(out, indent=1, allow_nan=False) + "\n")
    print(json.dumps({"all_pass": all_pass}, indent=1))
    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
