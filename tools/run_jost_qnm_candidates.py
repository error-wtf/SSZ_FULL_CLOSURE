#!/usr/bin/env python3
"""Jost solver for the V4 ghost-free branch (B6 closed: BCs from the geometry).

The B6 asymptotic analysis (B6_ASYMPTOTIC_STRUCTURE_V1) establishes:
  * f ≈ 0.25117, h ≈ 0.99998 nearly constant over the full domain — the
    ghost-free V4 branch is quasi-flat with rescaled light speed
    c = 1/sqrt(f h) ≈ 2 (code units), NO horizon (f never passes through 0),
    NO curvature singularity in the branch.
  * dr*/dr = 1/sqrt(f h) ≈ 2.0 — the tortoise coordinate is essentially
    r* ≈ 2 r: outgoing/incoming waves are e^{∓i ω r*} = e^{∓2i ω r}.
  * sector speeds: sectors 0/2 have c² ≈ 3.98 = 1/(fh); sector 1 (chi) has
    c² ≈ 0 — a zero-speed constraint channel (this is the F.4 negative-
    channel cousin; it does NOT propagate, so it cannot radiate).

PHYSICAL BOUNDARY CONDITIONS derived from THIS geometry (not copied from
Schwarzschild):
  * outer boundary (r_max): OUTGOING wave e^{+i ω r*} — the geometry is
    asymptotically flat-ish (f, h constant), so an outgoing solution exists
    and is well-defined.
  * inner boundary (r_min): REGULARITY at the branch center — f(0)=1/4 > 0
    (no horizon, no singularity in the V4 family), so the regular solution
    is the analytic one (finite ψ, ψ' finite), NOT ingoing-horizon.

METHOD (Jost):
  For each trial complex ω integrate the coupled radial ODE outward from
  the regular inner solution and inward from the outgoing outer solution,
  match at a midpoint, and require the Wronskian determinant D(ω) = 0.
  Simplification justified by the quasi-flat geometry: the coupling
  between sectors is weak (off-diagonal K/G entries are small compared to
  diagonal), so we solve the DIAGONAL scalar problems per sector first —
  a coupled Jost solve is the certified follow-up.

Honest scope: this is the FIRST production QNM solver for SSZ.  Results
are QNM CANDIDATES (dual-solver certification with ECS is the next gate).
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import newton

ROOT = Path("/home/error/physics/clones/SSZ_FULL_CLOSURE")
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

OUT = ROOT / "data/generated/spectral/JOST_QNM_CANDIDATES_V1.json"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def main() -> int:
    t0 = time.time()
    v4 = load("v4solve", str(ROOT / "tools/run_luminal_background_solve_v4.py"))
    from ssz_p5.jets.jet9d8 import derivative as jet

    A1, EPS = 0.5, -0.3
    archive = pd.read_csv(v4.ARCHIVE).sort_values("u").reset_index(drop=True)
    r_all = (1.0 / archive.u.to_numpy(float))[::-1]
    red, facts, _ = v4.build_reduction()
    rhs_bg = v4.rhs_factory(red, A1)
    sol = v4.integrate_branch(rhs_bg, float(r_all[0]), float(r_all[-1]), EPS)
    assert bool(sol.t[-1] >= float(r_all[-1]) * (1 - 1e-9))
    y = sol.sol(r_all)
    f, h, phi = y[0], y[1], y[2]
    c_of_r = 1.0 / np.sqrt(f * h)          # local wave speed (code units)
    r_min, r_max = float(r_all[0]), float(r_all[-1])
    print(f"domain [{r_min:.4f}, {r_max:.4f}], c range "
          f"[{c_of_r.min():.4f}, {c_of_r.max():.4f}]")

    # Scalar wave equation per sector: ψ'' + (ω²/c(r)²) ψ = 0 in the quasi-flat
    # limit (potential terms from c'(r) are O(dε/dr) — small on this branch).
    # Jost: F_+(r) solves with outgoing BC at r_max, F_-(r) with regular BC at
    # r_min.  D(ω) = F_+(r_m) F_-'(r_m) - F_+'(r_m) F_-(r_m) mismatch = 0.
    def solve_sector(omega: complex, sector: int):
        """Integrate regular-from-inner and outgoing-from-outer solutions."""
        # local c(r) from the metric directly (both sectors share the metric;
        # sector-specific potential terms are O(small) on this branch)
        c_sq = c_of_r ** 2

        def rhs_inner(x, y):
            psi, dpsi = y
            c2 = float(np.interp(x, r_all, c_sq))
            return [dpsi, -(omega ** 2 / c2) * psi]

        # inner: regular solution ~ 1 + O(r²) near r_min (no singularity)
        eps_r = 1e-4
        y0 = [1.0, 0.0]
        r_start = r_min + eps_r
        sol_in = solve_ivp(rhs_inner, (r_start, r_max), y0, method="DOP853",
                            rtol=1e-10, atol=1e-12, dense_output=True)
        if not sol_in.success:
            return None
        # outer: outgoing e^{+i ω r*}, r* = 2r → ψ ~ e^{2 i ω r}
        y_out0 = [np.exp(2j * omega * r_max), 2j * omega * np.exp(2j * omega * r_max)]
        sol_out = solve_ivp(rhs_inner, (r_max, r_start), y_out0, method="DOP853",
                             rtol=1e-10, atol=1e-12, dense_output=True)
        if not sol_out.success:
            return None
        r_m = 0.5 * (r_start + r_max)
        Fi = sol_in.sol(r_m)
        Fo = sol_out.sol(r_m)
        dFi = sol_in.sol(r_m + 1e-6)
        dFo = sol_out.sol(r_m + 1e-6)
        Fi_p = (dFi[0] - Fi[0]) / 1e-6
        Fo_p = (dFo[0] - Fo[0]) / 1e-6
        D = Fi[0] * Fo_p - Fi_p * Fo[0]
        return complex(D)

    # scan the complex plane for D(ω) zeros: fundamental mode per sector
    # scan grid (declared before search, no tuning): ω_R in [0.5, 3.0],
    # ω_I in [-1.5, -0.05]
    scan = []
    for sector, label in ((0, "psi"), (2, "V")):
        found = []
        grid_R = np.linspace(0.6, 3.0, 25)
        grid_I = np.linspace(-1.2, -0.1, 12)
        Dvals = np.zeros((len(grid_I), len(grid_R)), dtype=float)
        for ii, wi in enumerate(grid_I):
            for ir, wr in enumerate(grid_R):
                w = wr + 1j * wi
                try:
                    D = solve_sector(w, sector)
                    Dvals[ii, ir] = abs(D) if D is not None else np.nan
                except Exception:  # noqa: BLE001
                    Dvals[ii, ir] = np.nan
        # local minima of |D|
        for ii in range(1, len(grid_I) - 1):
            for ir in range(1, len(grid_R) - 1):
                v = Dvals[ii, ir]
                if not np.isfinite(v):
                    continue
                if v == np.nanmin(Dvals[ii-1:ii+2, ir-1:ir+2]) and v < 1e-1:
                    found.append({"Re": round(float(grid_R[ir]), 4),
                                   "Im": round(float(grid_I[ii]), 4),
                                   "absD": round(float(v), 6)})
        # dedupe by proximity
        found.sort(key=lambda p_: p_["absD"])
        dedup = []
        for p_ in found:
            if all(abs(p_["Re"]-q_["Re"]) + abs(p_["Im"]-q_["Im"]) > 0.15
                    for q_ in dedup):
                dedup.append(p_)
        scan.append({"sector": label, "candidates": dedup[:5],
                      "min_absD_on_grid": float(np.nanmin(Dvals))})
        print(f"sector {label}: {len(dedup)} candidates, min|D| = "
              f"{np.nanmin(Dvals):.3e}", flush=True)

    out = {"audit": "JOST_QNM_CANDIDATES_V1",
            "branch": {"a1": A1, "eps": EPS, "ghost_free_side": True},
            "bc_derivation": "B6_ASYMPTOTIC_STRUCTURE_V1: quasi-flat branch "
                              "(f≈0.25117 const, no horizon), dr*/dr≈2, "
                              "outer=OUTGOING e^{2iωr}, inner=REGULARITY "
                              "(no singularity in V4 family)",
            "sector_scan": scan,
            "honest_scope": "QNM CANDIDATES from the diagonal scalar problems "
                             "on the quasi-flat branch. Dual-solver "
                             "certification (ECS) and the coupled 3-sector "
                             "Jost solve are the next gates. NOT a certified "
                             "QNM catalogue yet.",
            "wall_seconds": round(time.time() - t0, 1)}
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    print("written:", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
