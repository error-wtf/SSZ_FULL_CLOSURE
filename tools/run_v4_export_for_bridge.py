#!/usr/bin/env python3
"""Produce a V4-branch spectral export in the frozen-export schema.

Emits the ghost-free V4 branch (a1=+0.5, eps=-0.3) through the certified
chain: solve -> 41-slot emission -> canonical_audit (L ladder) -> antisym-
metric identity projection on S (declared in V4_S_PROJECTION_REPAIR_V1,
validator UNCHANGED) -> npz in the exact FROZEN_SPECTRAL_OPERATOR_EXPORT_V1
array schema ({L}_r/K/R/G/S/M/Dh1/DeltaV/pivotA0) + provenance JSON with
SHA-256 sidecar.

The export is produced by the CLOSURE repo (producer side).  The bridge
(validator side) consumes it without importing any closure code.
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

OUT_NPZ = ROOT / "data/generated/spectral/V4_GHOSTFREE_BRANCH_EXPORT_V1.npz"
OUT_JSON = ROOT / "data/generated/spectral/V4_GHOSTFREE_BRANCH_EXPORT_V1.json"
LADDER = (6, 12, 20)


def load(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> int:
    t0 = time.time()
    red9 = load("red9", str(ROOT / "src/ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py"))
    v4 = load("v4solve", str(ROOT / "tools/run_luminal_background_solve_v4.py"))
    from ssz_p5.reducer.canonical import validate_operator, ReducedOperator, ConstraintPivots
    from ssz_p5.jets.jet9d8 import derivative as jet

    A1, EPS = 0.5, -0.3
    WINDOW_R = (0.05, 1.35)

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

    arrays = {}
    proj_max = 0.0
    for L in LADDER:
        res = red9.canonical_audit(sub, L)
        S = 0.5 * (res["S"] - np.swapaxes(res["S"], 1, 2))
        proj_max = max(proj_max, float(np.max(np.abs(res["S"] - S))))
        m = res["maps"]
        op = ReducedOperator(L, sub.x.to_numpy(float),
                             res["K"], res["R"], res["G"], S, res["M"],
                             ConstraintPivots(m["Dh1"], m["DeltaV"], m["pivotA0"]),
                             f"v4_ghostfree_L{L}")
        validate_operator(op)
        print(f"L={L}: validator PASS (projection removed {proj_max:.3e} max)", flush=True)
        arrays[f"{L}_r"] = sub.x.to_numpy(float)
        arrays[f"{L}_K"] = res["K"]
        arrays[f"{L}_R"] = res["R"]
        arrays[f"{L}_G"] = res["G"]
        arrays[f"{L}_S"] = S
        arrays[f"{L}_M"] = res["M"]
        arrays[f"{L}_Dh1"] = m["Dh1"]
        arrays[f"{L}_DeltaV"] = m["DeltaV"]
        arrays[f"{L}_pivotA0"] = m["pivotA0"]

    np.savez(OUT_NPZ, **arrays)
    digest = sha256_file(OUT_NPZ)
    prov = {
        "export": "V4_GHOSTFREE_BRANCH_EXPORT_V1",
        "schema": "FROZEN_SPECTRAL_OPERATOR_EXPORT_V1 array schema",
        "branch": {"a1": A1, "eps": EPS, "ghost_free_side": A1 * EPS < 0},
        "window_r": [float(r_w[0]), float(r_w[-1])],
        "window_rows": int(len(r_w)),
        "L_ladder": list(LADDER),
        "S_projection": "antisymmetric identity projection S := (S - S^T)/2, "
                         "declared in V4_S_PROJECTION_REPAIR_V1; validator UNCHANGED",
        "max_S_projection_abs": proj_max,
        "background_residuals": "E00/E11/E22 certified in N2_LUMINAL_SOLVE_RESULT_V4",
        "npz_sha256": digest,
        "producer": "SSZ_FULL_CLOSURE (closure branch spectroscopy-real-data-20261004)",
        "validator_note": "consumer (SSZ-Transport-Bridge) imports NO closure code; "
                          "hash-bound npz only",
        "wall_seconds": round(time.time() - t0, 1),
    }
    OUT_JSON.write_text(json.dumps(prov, indent=1) + "\n")
    print("npz:", OUT_NPZ.name, digest[:16], "| provenance written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
