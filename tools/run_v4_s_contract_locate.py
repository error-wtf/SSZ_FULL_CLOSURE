#!/usr/bin/env python3
"""Locate the EXACT grid points where the S-antisymmetry residual exceeds the
policy threshold on the ghost-free V4 branch (a1=+0.5, eps=-0.5, window
[0.05, 1.35]) — the F4 blocker, dissected.

Declared diagnostic. No tolerance changes, no validator edits.
"""
from __future__ import annotations

import importlib.util
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
    window_r = (0.05, 1.35)
    a1v, eps = 0.5, -0.5

    archive = pd.read_csv(v4.ARCHIVE).sort_values("u").reset_index(drop=True)
    r_all = (1.0 / archive.u.to_numpy(float))[::-1]
    lo = int(np.searchsorted(r_all, window_r[0]))
    hi = int(np.searchsorted(r_all, window_r[1]))
    win = slice(lo + 8, hi - 8)  # declared interior margin as in F4 tool
    r_w = r_all[win]

    red, facts, _ = v4.build_reduction()
    rhs = v4.rhs_factory(red, a1v)
    sol = v4.integrate_branch(rhs, float(r_all[0]), float(r_all[-1]), eps)
    y = sol.sol(r_w)
    f_w, h_w, phi_w = y[0], y[1], y[2]
    from ssz_p5.jets.jet9d8 import derivative as jet
    phi_r = jet(r_w, phi_w)

    from ssz_p5.coefficients.mh_action_primitives import quartic_g5zero_primitives
    from ssz_p5.coefficients.mh_general_primitives import emit_from_primitives
    feed = pd.DataFrame({
        "x": r_w, "f": f_w, "h": h_w, "phi_r": phi_r,
        "G4": 0.5 + a1v * (phi_w - 1.0),
        "G4X": np.zeros(len(r_w)), "G4XX": np.zeros(len(r_w)),
        "G4phi": np.full(len(r_w), a1v), "G4phiX": np.zeros(len(r_w)),
        "G3X": np.zeros(len(r_w)),
    })
    prim = quartic_g5zero_primitives(feed)
    c2 = np.sqrt(f_w * h_w) * phi_r * 0.5 * r_w**2
    prim_in = pd.DataFrame({
        "x": r_w, "f": f_w, "h": h_w, "phiprime": phi_r,
        "A0prime": np.zeros(len(r_w)),
        "a1": prim.a1_action.to_numpy(float), "c2": c2,
        "c4": prim.c4_action.to_numpy(float),
        "F_tensor": prim.F_tensor_action.to_numpy(float),
        "G_tensor": prim.G_tensor_action.to_numpy(float),
        "H_tensor": prim.H_tensor_action.to_numpy(float),
    })
    stream = emit_from_primitives(prim_in)
    stream["phi"] = phi_w
    sub = stream.sort_values("x").reset_index(drop=True)
    if "X" not in sub.columns:
        sub["X"] = -sub["h"].to_numpy(float) * sub["phiprime"].to_numpy(float) ** 2 / 2.0

    import importlib.util as ilu
    sp2 = ilu.spec_from_file_location(
        "red9", ROOT / "src/ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    red9 = ilu.module_from_spec(sp2)
    sys.modules["red9"] = red9
    sp2.loader.exec_module(red9)
    res = red9.canonical_audit(sub, 6)
    S = res["S"]
    asym = S + np.swapaxes(S, 1, 2)
    P01 = res["P"].get((0, 1))
    P02 = res["P"].get((0, 2))
    P02p = np.empty_like(P02)
    from ssz_p5.jets.jet9d8 import derivative
    for i in range(3):
        for j in range(3):
            P02p[:, i, j] = derivative(sub.x.to_numpy(float), P02[:, i, j], 1)

    scale = max(1.0, float(np.max(np.abs(S))))
    err = np.max(np.abs(asym), axis=(1, 2)) / scale
    bad = np.where(err > 1e-7)[0]
    print(f"rows={len(sub)}  scaled err max={err.max():.3e}  violations(>1e-7)={len(bad)}")
    if len(bad):
        print("first violation at r =", sub.x.to_numpy(float)[bad[0]])
        print("last  violation at r =", sub.x.to_numpy(float)[bad[-1]])
        print("worst  at r =", sub.x.to_numpy(float)[bad[np.argmax(err[bad])]],
              "err =", err[bad].max())
        # distribution
        for lo_b, hi_b in [(0.0, 0.1), (0.1, 0.3), (0.3, 0.6), (0.6, 1.0), (1.0, 1.4)]:
            m = (sub.x.to_numpy(float) >= lo_b) & (sub.x.to_numpy(float) < hi_b)
            if m.any():
                print(f"  r in [{lo_b},{hi_b}): max err = {err[m].max():.3e}  ({m.sum()} rows)")
        # the by-parts identity residual at the worst point
        wi = bad[np.argmax(err[bad])]
        ident = 0.5 * (P01[wi] + np.swapaxes(P01[wi], 0, 1)) - P02p[wi]
        print("by-parts identity residual at worst point (matrix):")
        print(np.array2string(ident, precision=3, suppress_small=False))
        # S[1,1] diagonal check — S antisymmetry REQUIRES zero diagonal:
        print("S diagonal max |S[i,i]| =", float(np.max(np.abs(np.diagonal(S, axis1=1, axis2=2)))))
        print("asym diag[1,1] at worst:", float(asym[wi, 1, 1]),
              "| P01[1,1]:", float(P01[wi, 1, 1]),
              "| P02'[1,1]:", float(P02p[wi, 1, 1]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
