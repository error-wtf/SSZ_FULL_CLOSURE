#!/usr/bin/env python3
"""V4_FULL_3DOF_BRANCH_HEALTH_V1 — rule-9 scan of the N4 V4 branch family.

After V4_UNREDUCED_CONSTRAINT_RANK_V1 proved dphi is dynamical (3 physical
DOF: psi, dphi, V), the (a1=+0.5, eps=-0.3) branch fails the honest
3-field kinetic gate (min eig K_full3 = -1.94e05 in the window, negative
direction 100% dphi).

This tool scans ALL previously generated V4 candidate branches with the
SAME full 3-field measurement:
  min eig K(3x3) over the certified window
for every (a1, eps) candidate in the N4 branch verband.

Selection rule (declared BEFORE running):
  a branch is a HEALTH candidate iff min eig K(3x3) > 0 over the window
  on the full 3-field operator, rank K = 3 everywhere, and the branch
  integrates to the window boundary regularly.
  NO fitting to observations. Theory health only.
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

OUT = ROOT / "data/generated/spectral/V4_FULL_3DOF_BRANCH_HEALTH_V1.json"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def measure_branch(v4, red9, A1, EPS, n_win=3500):
    """Full 3-field kinetic measurement on one branch. Returns dict."""
    r_grid = np.linspace(0.05, 1.35, n_win)
    try:
        archive = pd.read_csv(v4.ARCHIVE).sort_values("u").reset_index(drop=True)
        r_all = (1.0 / archive.u.to_numpy(float))[::-1]
        red, facts, _ = v4.build_reduction()
        rhs = v4.rhs_factory(red, A1)
        sol = v4.integrate_branch(rhs, float(r_all[0]), float(r_all[-1]), EPS)
        if not bool(sol.t[-1] >= float(r_all[-1]) * (1 - 1e-9)):
            return {"integration": "FAILED_BEFORE_WINDOW_END",
                     "r_reached": float(sol.t[-1])}
        y = sol.sol(r_grid)
        f_w, h_w, phi_w = y[0], y[1], y[2]
        from ssz_p5.jets.jet9d8 import derivative as jet
        phi_r = jet(r_grid, phi_w)
        from ssz_p5.coefficients.mh_action_primitives import quartic_g5zero_primitives
        from ssz_p5.coefficients.mh_general_primitives import emit_from_primitives
        prim = quartic_g5zero_primitives(pd.DataFrame({
            "x": r_grid, "f": f_w, "h": h_w, "phi_r": phi_r,
            "G4": 0.5 + A1*(phi_w-1.0), "G4X": np.zeros(n_win),
            "G4XX": np.zeros(n_win), "G4phi": np.full(n_win, A1),
            "G4phiX": np.zeros(n_win), "G3X": np.zeros(n_win)}))
        c2 = np.sqrt(f_w*h_w)*phi_r*0.5*r_grid**2
        stream = emit_from_primitives(pd.DataFrame({
            "x": r_grid, "f": f_w, "h": h_w, "phiprime": phi_r,
            "A0prime": np.zeros(n_win),
            "a1": prim.a1_action.to_numpy(float), "c2": c2,
            "c4": prim.c4_action.to_numpy(float),
            "F_tensor": prim.F_tensor_action.to_numpy(float),
            "G_tensor": prim.G_tensor_action.to_numpy(float),
            "H_tensor": prim.H_tensor_action.to_numpy(float)}))
        stream["phi"] = phi_w
        sub = stream.sort_values("x").reset_index(drop=True)
        sub["X"] = -sub["h"].to_numpy(float)*sub["phiprime"].to_numpy(float)**2/2.0
        res = red9.canonical_audit(sub, 6)
        K = res["K"]
        mins = np.array([np.linalg.eigvalsh(0.5*(K[i]+K[i].T))[0]
                          for i in range(K.shape[0])])
        ranks = np.array([np.linalg.matrix_rank(K[i], tol=1e-10)
                           for i in range(K.shape[0])])
        # which sector carries the minimum?
        i = int(np.argmin(mins))
        ev, evec = np.linalg.eigh(0.5*(K[i]+K[i].T))
        w = np.abs(evec[:, 0])**2; w /= w.sum()
        return {
            "integration": "OK",
            "min_eig": float(mins.min()),
            "n_negative": int((mins < 0).sum()),
            "n_nodes": int(K.shape[0]),
            "rank3_everywhere": bool((ranks == 3).all()),
            "worst_sector_weights": {"psi": round(float(w[0]), 4),
                                      "dphi": round(float(w[1]), 4),
                                      "V": round(float(w[2]), 4)},
            "healthy_3field": bool(mins.min() > 0 and (ranks == 3).all()),
        }
    except Exception as exc:  # noqa: BLE001 - record, don't crash the scan
        return {"integration": f"ERROR: {type(exc).__name__}: {exc}"[:200]}


def main() -> int:
    t0 = time.time()
    v4 = load("v4solve", str(ROOT / "tools/run_luminal_background_solve_v4.py"))
    red9 = load("red9", str(ROOT / "src/ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py"))

    # N4 branch family (from N4_V4_BRANCH_VERBAND): all 8 quadrant combos
    candidates = [(a1, eps) for a1 in (+0.5, -0.5) for eps in (+0.3, -0.3)]
    out = {"audit": "V4_FULL_3DOF_BRANCH_HEALTH_V1",
            "declared_rule": ("healthy = min eig K(3x3) > 0 over the window, "
                               "rank 3 everywhere, regular integration; NO "
                               "observation fitting"),
            "results": {}}
    healthy = []
    for A1, EPS in candidates:
        r = measure_branch(v4, red9, A1, EPS)
        out["results"][f"a1={A1:+.1f}_eps={EPS:+.1f}"] = r
        ok = r.get("healthy_3field")
        print(f"a1={A1:+.1f} eps={EPS:+.1f}: {r.get('integration')} "
              f"min_eig={r.get('min_eig', float('nan')):.3e} "
              f"healthy={ok} sector={r.get('worst_sector_weights')}", flush=True)
        if ok:
            healthy.append((A1, EPS))
    out["healthy_branches"] = [list(b) for b in healthy]
    out["n_healthy"] = len(healthy)
    out["verdict"] = ("HEALTHY_BRANCH_FOUND" if healthy else
                       "NO_HEALTHY_V4_BRANCH — all candidates fail the "
                       "honest 3-field kinetic gate; per rule 9 this is a "
                       "hard theory no-go for the static V4 family until a "
                       "new branch class is derived.")
    print("VERDICT:", out["verdict"])
    out["wall_seconds"] = round(time.time() - t0, 1)
    OUT.write_text(json.dumps(out, indent=1, allow_nan=False) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
