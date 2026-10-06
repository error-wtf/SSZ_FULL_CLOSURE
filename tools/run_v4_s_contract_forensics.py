#!/usr/bin/env python3
"""V4 S-contract forensics (declared diagnostic, no tolerance changes).

Branch: (a1, eps) = (+0.5, -0.5), ghost-free side, K >= 0 measured (N4).

Steps:
  1. re-solve the branch with the production V4 solver module (same code
     path as the N4 verdict),
  2. emit the 41-slot stream through the verified chain,
  3. run the production reducer's canonical_audit at L=6 and measure
     S^T + S per (i,j) pair and per r-region,
  4. measure the underlying by-parts identity 0.5*(P01+P01^T) - P02'
     (radial_first_sym_res) — the mathematical root of S-antisymmetry,
  5. compare against the LOCKED constant-G4 golden stream.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/home/error/physics/clones/SSZ_FULL_CLOSURE")
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

spec = importlib.util.spec_from_file_location(
    "v4solve", ROOT / "tools/run_luminal_background_solve_v4.py")
v4 = importlib.util.module_from_spec(spec)
sys.modules["v4solve"] = v4
spec.loader.exec_module(v4)


def main() -> int:
    out = {"audit": "V4_S_CONTRACT_FORENSICS_V1",
           "branch": {"a1": 0.5, "eps": -0.5, "ghost_free_side": True}}

    # --- 1. resolve the branch (production code path)
    archive = pd.read_csv(v4.ARCHIVE).sort_values("u").reset_index(drop=True)
    u_all = archive.u.to_numpy(float)
    r_all = (1.0 / u_all)[::-1]
    r_min, r_core = float(r_all[0]), float(r_all[-1])
    n_out = len(r_all)

    red, facts, _ = v4.build_reduction()
    assert all(facts.values())
    rhs = v4.rhs_factory(red, 0.5)
    sol = v4.integrate_branch(rhs, r_min, r_core, -0.5)
    reached = bool(sol.t[-1] >= r_core * (1 - 1e-9))
    out["branch"]["reached_core"] = reached
    if not reached:
        print("branch did not reach core:", sol.t[-1])
        return 1
    y = sol.sol(r_all)
    f_p, h_p, phi_p = y[0], y[1], y[2]
    from ssz_p5.jets.jet9d8 import derivative as jet
    phi_r = jet(r_all, phi_p)
    out["branch"]["phi_range"] = [float(phi_p.min()), float(phi_p.max())]
    print(f"branch solved: phi in [{phi_p.min():.4f}, {phi_p.max():.4f}]")

    # --- 2. emit the 41-slot stream (verified chain)
    from ssz_p5.coefficients.mh_action_primitives import quartic_g5zero_primitives
    from ssz_p5.coefficients.mh_general_primitives import emit_from_primitives

    feed = pd.DataFrame({
        "x": r_all, "f": f_p, "h": h_p, "phi_r": phi_r,
        "G4": 0.5 + 0.5 * (phi_p - 1.0),
        "G4X": np.zeros(n_out), "G4XX": np.zeros(n_out),
        "G4phi": np.full(n_out, 0.5), "G4phiX": np.zeros(n_out),
        "G3X": np.zeros(n_out),
    })
    prim = quartic_g5zero_primitives(feed)
    c2 = np.sqrt(f_p * h_p) * phi_r * 0.5 * r_all**2
    prim_in = pd.DataFrame({
        "x": r_all, "f": f_p, "h": h_p, "phiprime": phi_r,
        "A0prime": np.zeros(n_out),
        "a1": prim.a1_action.to_numpy(float),
        "c2": c2,
        "c4": prim.c4_action.to_numpy(float),
        "F_tensor": prim.F_tensor_action.to_numpy(float),
        "G_tensor": prim.G_tensor_action.to_numpy(float),
        "H_tensor": prim.H_tensor_action.to_numpy(float),
    })
    stream = emit_from_primitives(prim_in)
    print("stream emitted:", len(stream), "rows")

    # --- 3. production reducer canonical_audit at L=6
    import importlib.util as _ilu
    _spec = _ilu.spec_from_file_location(
        "jet9d8_reducer",
        ROOT / "src/ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    _red = _ilu.module_from_spec(_spec)
    sys.modules["jet9d8_reducer"] = _red
    _spec.loader.exec_module(_red)
    canonical_audit = _red.canonical_audit
    df = stream  # reducer reads df.x
    res = canonical_audit(df, 6)
    S = res["S"]
    asym = S + np.swapaxes(S, 1, 2)  # should be ~0
    scale = max(1.0, float(np.max(np.abs(S))))
    diag = res["diagnostics"]
    out["reducer_diagnostics"] = {k: float(v) for k, v in diag.items()
                                  if isinstance(v, (int, float))}
    out["S_asym_scaled_max"] = float(np.max(np.abs(asym)) / scale)
    # per field pair
    pairs = {}
    for i in range(3):
        for j in range(3):
            pairs[f"S[{i},{j}]+S[{j},{i}]_max"] = float(np.max(np.abs(asym[:, i, j])))
    out["S_asym_per_pair"] = pairs
    # per region: interior thirds
    n = len(r_all)
    regions = {"inner": slice(0, n // 3), "middle": slice(n // 3, 2 * n // 3),
               "outer": slice(2 * n // 3, n)}
    out["S_asym_per_region"] = {
        k: float(np.max(np.abs(asym[s]))) for k, s in regions.items()
    }
    print(f"S asym scaled max: {out['S_asym_scaled_max']:.3e}")
    print("per region:", out["S_asym_per_region"])

    # --- 4. the by-parts identity: 0.5*(P01+P01^T) - P02'
    # rebuild P01/P02 the way the reducer does (reduced_operator needs the
    # full reduced frame; use its internal path via reduced_operator)
    # Simpler and equivalent: the diagnostic max_radial_first_sym_res IS this
    # quantity. Report it directly:
    out["radial_first_sym_res"] = float(diag["max_radial_first_sym_res"])
    out["S_antisym_abs_max"] = float(diag["max_S_sym"])

    # --- 5. golden constant-G4 comparison
    gold_path = ROOT / "data/production/ssz_p5_integrable_svt_unreduced_even_coefficients_2026-09-12.csv"
    if gold_path.exists():
        gold = pd.read_csv(gold_path)
        try:
            res_g = canonical_audit(gold, 6)
            dg = res_g["diagnostics"]
            out["golden_radial_first_sym_res"] = float(dg["max_radial_first_sym_res"])
            out["golden_S_antisym_abs_max"] = float(dg["max_S_sym"])
            Sg = res_g["S"]
            scale_g = max(1.0, float(np.max(np.abs(Sg))))
            out["golden_S_asym_scaled_max"] = float(
                np.max(np.abs(Sg + np.swapaxes(Sg, 1, 2))) / scale_g)
        except Exception as exc:  # noqa: BLE001
            out["golden_error"] = f"{type(exc).__name__}: {exc}"[:200]

    Path("/home/error/physics/clones/SSZ_FULL_CLOSURE/data/generated/spectral/V4_S_CONTRACT_FORENSICS_V1.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1)[:1800])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
