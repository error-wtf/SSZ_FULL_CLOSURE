#!/usr/bin/env python3
"""V4_FULL_3DOF_OPERATOR_V2 — rule-5 export: run the EXISTING action-level
reducer directly (NO post-hoc Schur) on the healthy production branch
(a1=-0.5, eps=-0.3) and export the FULL 3-field K/R/G/S/M matrices in the
(psi, dphi, V) basis over the entire certified domain, for the whole
L ladder 6,12,20,42,110,420,1000.

Also measures (rule 6, no interpretation by label):
  - rank(K) per node
  - eigenvalues(K) per node
  - rank(G) per node
  - constraint pivots
and reports min lambda1/2/3 with locations per L.

Gate (rule 7): rank K = 3 AND min eig K > 0 over the domain for the
L ladder -> branch is three-channel healthy -> solver proceeds with
Psi = (psi, dphi, V)^T.
"""
from __future__ import annotations

import hashlib
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

OUT_NPZ = ROOT / "data/generated/spectral/V4_FULL_3DOF_OPERATOR_V2.npz"
OUT_JSON = ROOT / "data/generated/spectral/V4_FULL_3DOF_OPERATOR_V2.json"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def main() -> int:
    t0 = time.time()
    v4 = load("v4solve", str(ROOT / "tools/run_luminal_background_solve_v4.py"))
    red9 = load("red9", str(ROOT / "src/ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py"))
    from ssz_p5.jets.jet9d8 import derivative as jet

    A1, EPS = -0.5, -0.3          # healthy production branch (rule 9 scan)
    LADDER = [6, 12, 20, 42, 110, 420, 1000]
    r_grid = np.linspace(0.05, 1.35, 3500)

    # ---------- regenerate the branch from action inputs (not stale NPZ)
    archive = pd.read_csv(v4.ARCHIVE).sort_values("u").reset_index(drop=True)
    r_all = (1.0 / archive.u.to_numpy(float))[::-1]
    red, facts, _ = v4.build_reduction()
    assert all(facts.values())
    rhs = v4.rhs_factory(red, A1)
    sol = v4.integrate_branch(rhs, float(r_all[0]), float(r_all[-1]), EPS)
    assert bool(sol.t[-1] >= float(r_all[-1]) * (1 - 1e-9))
    y = sol.sol(r_grid)
    f_w, h_w, phi_w = y[0], y[1], y[2]
    phi_r = jet(r_grid, phi_w)
    from ssz_p5.coefficients.mh_action_primitives import quartic_g5zero_primitives
    from ssz_p5.coefficients.mh_general_primitives import emit_from_primitives
    prim = quartic_g5zero_primitives(pd.DataFrame({
        "x": r_grid, "f": f_w, "h": h_w, "phi_r": phi_r,
        "G4": 0.5 + A1*(phi_w-1.0), "G4X": np.zeros(len(r_grid)),
        "G4XX": np.zeros(len(r_grid)), "G4phi": np.full(len(r_grid), A1),
        "G4phiX": np.zeros(len(r_grid)), "G3X": np.zeros(len(r_grid))}))
    c2 = np.sqrt(f_w*h_w)*phi_r*0.5*r_grid**2
    stream = emit_from_primitives(pd.DataFrame({
        "x": r_grid, "f": f_w, "h": h_w, "phiprime": phi_r,
        "A0prime": np.zeros(len(r_grid)),
        "a1": prim.a1_action.to_numpy(float), "c2": c2,
        "c4": prim.c4_action.to_numpy(float),
        "F_tensor": prim.F_tensor_action.to_numpy(float),
        "G_tensor": prim.G_tensor_action.to_numpy(float),
        "H_tensor": prim.H_tensor_action.to_numpy(float)}))
    stream["phi"] = phi_w
    sub = stream.sort_values("x").reset_index(drop=True)
    sub["X"] = -sub["h"].to_numpy(float)*sub["phiprime"].to_numpy(float)**2/2.0

    # ---------- run the EXISTING action-level reducer per L (no Schur)
    store = {}
    summary = {}
    all_healthy = True
    for L in LADDER:
        res = red9.canonical_audit(sub, L)
        K, G = res["K"], res["G"]
        S = 0.5*(res["S"] - np.swapaxes(res["S"], 1, 2))
        M, Rm = res["M"], res["R"]
        store[f"L{L}"] = {"K": K, "R": Rm, "G": G, "S": S, "M": M,
                           "maps": {k: res["maps"][k]
                                     for k in ("Dh1", "DeltaV", "pivotA0")}}
        # rule-6 measurements
        eigs = np.array([np.linalg.eigvalsh(0.5*(K[i]+K[i].T))
                          for i in range(K.shape[0])])
        ranks = np.array([np.linalg.matrix_rank(K[i], tol=1e-10)
                           for i in range(K.shape[0])])
        rankG = np.array([np.linalg.matrix_rank(G[i], tol=1e-10)
                           for i in range(G.shape[0])])
        i1, i2, i3 = (int(np.argmin(eigs[:, j])) for j in range(3))
        piv = {"Dh1_min_abs": float(np.min(np.abs(res["maps"]["Dh1"]))),
                "DeltaV_min_abs": float(np.min(np.abs(res["maps"]["DeltaV"]))),
                "pivotA0_min_abs": float(np.min(np.abs(res["maps"]["pivotA0"])))}
        healthy = bool(eigs.min() > 0 and (ranks == 3).all())
        all_healthy = all_healthy and healthy
        summary[f"L{L}"] = {
            "rankK3_everywhere": bool((ranks == 3).all()),
            "rankG3_everywhere": bool((rankG == 3).all()),
            "min_lambda1": float(eigs[:, 0].min()),
            "min_lambda1_at_r": float(r_grid[i1]),
            "min_lambda2": float(eigs[:, 1].min()),
            "min_lambda2_at_r": float(r_grid[i2]),
            "min_lambda3": float(eigs[:, 2].min()),
            "min_lambda3_at_r": float(r_grid[i3]),
            "pivots_min_abs": piv,
            "healthy": healthy,
        }
        print(f"L={L}: min eig = {eigs.min():.3e} @ r={r_grid[np.unravel_index(np.argmin(eigs), eigs.shape)[0]]:.3f} "
              f"| rank3={summary[f'L{L}']['rankK3_everywhere']} | healthy={healthy}", flush=True)

    out = {
        "audit": "V4_FULL_3DOF_OPERATOR_V2",
        "branch": {"a1": A1, "eps": EPS},
        "basis": ["psi", "dphi", "V"],
        "no_posthoc_schur": True,
        "reducer": "ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.canonical_audit",
        "domain": [0.05, 1.35], "n_nodes": int(len(r_grid)),
        "L_ladder": LADDER,
        "per_L_summary": summary,
        "all_healthy": bool(all_healthy),
        "gate_rule7": ("rank K = 3 AND min eig K > 0 over the domain/L ladder"
                        if all_healthy else "RULE 7 FAIL"),
        "verdict": ("THREE_CHANNEL_HEALTHY — solver proceeds with "
                     "Psi = (psi, dphi, V)^T, no channel reduction."
                     if all_healthy else
                     "NOT THREE-CHANNEL HEALTHY — see per_L_summary"),
        "wall_seconds": round(time.time() - t0, 1),
    }

    # ---------- export
    arrays = {"r": r_grid, "L_ladder": np.array(LADDER)}
    for L in LADDER:
        for blk in ("K", "R", "G", "S", "M"):
            arrays[f"L{L}_{blk}"] = store[f"L{L}"][blk]
        for piv in ("Dh1", "DeltaV", "pivotA0"):
            arrays[f"L{L}_{piv}"] = store[f"L{L}"]["maps"][piv]
    np.savez_compressed(OUT_NPZ, **arrays)
    h = hashlib.sha256(OUT_NPZ.read_bytes()).hexdigest()
    out["npz_sha256"] = h
    OUT_JSON.write_text(json.dumps(out, indent=1, allow_nan=False) + "\n")
    (OUT_NPZ.parent / (OUT_NPZ.name + ".sha256")).write_text(h + "\n")
    print("VERDICT:", out["verdict"][:120])
    return 0 if all_healthy else 1


if __name__ == "__main__":
    raise SystemExit(main())
