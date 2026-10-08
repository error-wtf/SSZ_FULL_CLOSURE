#!/usr/bin/env python3
"""V4_PHYSICAL_RESONANCE_EXPORT_V3_3DOF — re-runs canonical_audit and
exports ALL operator matrices the ODE needs: K, G, S, M, R (full, not
folded into S by assertion) plus the raw P-blocks and diagnostics.

R = -P11/2 per the canonical_audit convention. If max|R| = 0 this tool
PROVES it by hash; otherwise solvers must consume it.
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
ART = ROOT / "data/generated/spectral"
OUT_NPZ = ART / "V4_PHYSICAL_RESONANCE_EXPORT_V4_L12.npz"
OUT_JSON = ART / "V4_PHYSICAL_RESONANCE_EXPORT_V3_3DOF.json"


def main() -> int:
    """Re-runs the EXACT V2 export pipeline (same branch, same grid, same
    canonical_audit call) and exports K/G/S/M **and R** + raw P-blocks."""
    t0 = time.time()
    import pandas as pd
    spec = importlib.util.spec_from_file_location(
        "v4solve", str(ROOT / "tools/run_luminal_background_solve_v4.py"))
    v4 = importlib.util.module_from_spec(spec)
    sys.modules["v4solve"] = v4
    spec.loader.exec_module(v4)
    spec2 = importlib.util.spec_from_file_location(
        "reducer9d8",
        "/home/error/physics/clones/SSZ_FULL_CLOSURE/src/"
        "ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    red9 = importlib.util.module_from_spec(spec2)
    sys.modules["red9"] = red9
    spec2.loader.exec_module(red9)
    from ssz_p5.jets.jet9d8 import derivative as jet
    from ssz_p5.coefficients.mh_action_primitives import quartic_g5zero_primitives
    from ssz_p5.coefficients.mh_general_primitives import emit_from_primitives

    A1, EPS = -0.5, -0.3
    R_IN, R_WIN_END, R_EXT = 0.05, 1.0, 60.0
    N_WIN, N_EXT = 3500, 600
    L = 12

    archive = pd.read_csv(v4.ARCHIVE).sort_values("u").reset_index(drop=True)
    r_all = (1.0 / archive.u.to_numpy(float))[::-1]
    red, facts, _ = v4.build_reduction()
    assert all(facts.values())
    rhs = v4.rhs_factory(red, A1)
    sol = v4.integrate_branch(rhs, float(r_all[0]), R_EXT, EPS)
    assert bool(sol.t[-1] >= R_EXT * (1 - 1e-9)), "exterior continuation failed"
    r_win = np.linspace(R_IN, R_WIN_END, N_WIN, endpoint=False)
    r_ext = np.logspace(np.log10(R_WIN_END), np.log10(R_EXT), N_EXT + 1)[1:]
    r_grid = np.concatenate([r_win, r_ext])
    y = sol.sol(r_grid)
    f_w, h_w, phi_w = y[0], y[1], y[2]
    phi_r = jet(r_grid, phi_w)

    prim = quartic_g5zero_primitives(pd.DataFrame({
        "x": r_grid, "f": f_w, "h": h_w, "phi_r": phi_r,
        "G4": 0.5 + A1 * (phi_w - 1.0),
        "G4X": np.zeros(len(r_grid)), "G4XX": np.zeros(len(r_grid)),
        "G4phi": np.full(len(r_grid), A1), "G4phiX": np.zeros(len(r_grid)),
        "G3X": np.zeros(len(r_grid))}))
    c2 = np.sqrt(f_w * h_w) * phi_r * 0.5 * r_grid**2
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
    sub["X"] = -sub["h"].to_numpy(float) * sub["phiprime"].to_numpy(float) ** 2 / 2.0

    res = red9.canonical_audit(sub, float(L))
    K, G, S, M, R = res["K"], res["G"], res["S"], res["M"], res["R"]
    S = 0.5 * (S - np.swapaxes(S, 1, 2))  # antisymmetric part kept explicit
    n = len(r_grid)
    if K.shape[0] != n:
        raise SystemExit(f"grid mismatch: K {K.shape[0]} vs r {n}")

    np.savez_compressed(
        OUT_NPZ,
        r=r_grid, L=np.array([L]),
        K_phys=K, G_phys=G, S_phys=S, M_phys=M, R_phys=R,
        f=f_w, h=h_w, phi=phi_w,
        **{f"P_{k[0]}{k[1]}": v for k, v in res["P"].items()})
    hsh = hashlib.sha256(OUT_NPZ.read_bytes()).hexdigest()
    maxR = float(np.max(np.abs(R)))
    out = {
        "audit": "V4_PHYSICAL_RESONANCE_EXPORT_V3_3DOF",
        "branch": {"a1": A1, "eps": EPS},
        "L": L,
        "n_nodes": n,
        "grid": {"r_in": R_IN, "r_win_end": R_WIN_END, "r_ext": R_EXT,
                  "n_win": N_WIN, "n_ext": N_EXT},
        "npz_sha256": hsh,
        "R_block": {
            "max_abs": maxR,
            "max_abs_asymptotic_last100": float(np.max(np.abs(R[-100:]))),
            "identically_zero_proven": bool(maxR == 0.0),
            "note": ("R = -P11/2 full export on the exact V2 grid; "
                      "solvers MUST consume R from this artifact"),
        },
        "S_antisymmetric_kept": True,
        "diagnostics": {k: str(v) for k, v in res["diagnostics"].items()},
        "wall_seconds": round(time.time() - t0, 1),
    }
    OUT_JSON.write_text(json.dumps(out, indent=1, default=str) + "\n")
    print(json.dumps({k: out[k] for k in
                       ("npz_sha256", "R_block", "n_nodes")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
