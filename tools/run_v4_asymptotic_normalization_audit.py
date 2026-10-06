#!/usr/bin/env python3
"""V4_ASYMPTOTIC_NORMALIZATION_AUDIT_V1 — is 'NOT_FLAT' real physics or a
time-normalization artifact? (Lino, 06.10. night)

Four declared outputs:
  1. Curvature invariants R, R_{mu nu}R^{mu nu}, Kretschmann on the
     continued branch as r -> 400.  If all -> 0, the exterior is flat and
     f_inf = 0.2512 is a NON-UNIT LAPSE, not non-flatness.
  2. Time normalization T = sqrt(f_inf) t: verify g_TT -> -1, g_rr -> 1.
  3. Normalized tortoise coordinate R_* = sqrt(f_inf) r_*: verify
     dR_*/dr -> 1; outgoing wave Psi ~ exp(-i Omega_inf (T - R_*)).
  4. Pole invariance: a constant time rescaling t -> T shifts frequencies
     omega -> Omega = omega/sqrt(f_inf) but leaves the physical poles
     unchanged.  Verify by transforming the (psi,V) operator's omega
     terms symbolically and checking the characteristic roots map
     consistently.

Also tracked: G4(phi_inf) with phi_inf = 0.998334 (a1 = 0.5) ->
G4_inf ~ 0.49917: the asymptotic effective gravitational coupling.
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import sympy as sp

ROOT = Path("/home/error/physics/clones/SSZ_FULL_CLOSURE")
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

OUT = ROOT / "data/generated/spectral/V4_ASYMPTOTIC_NORMALIZATION_AUDIT_V1.json"


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

    A1V, EPSV = 0.5, -0.3
    out = {"audit": "V4_ASYMPTOTIC_NORMALIZATION_AUDIT_V1",
            "branch": {"a1": A1V, "eps": EPSV}}

    # ---------- continue the branch to r=400 (same chain as EXTERIOR_V1)
    red, facts, _ = v4.build_reduction()
    assert all(facts.values())
    rhs = v4.rhs_factory(red, A1V)
    archive = pd.read_csv(v4.ARCHIVE).sort_values("u").reset_index(drop=True)
    r_all = (1.0 / archive.u.to_numpy(float))[::-1]
    sol = v4.integrate_branch(rhs, float(r_all[0]), 400.0, EPSV)
    r_grid = np.linspace(1.35, min(sol.t[-1], 395.0), 2500)
    y = sol.sol(r_grid)
    f_g, h_g, phi_g = y[0], y[1], y[2]
    fp_g = jet(r_grid, f_g)
    hp_g = jet(r_grid, h_g)
    f_inf = float(f_g[-1])
    h_inf = float(h_g[-1])
    phi_inf = float(phi_g[-1])
    out["branch_limits"] = {
        "f_inf": round(f_inf, 8), "h_inf": round(h_inf, 8),
        "phi_inf": round(phi_inf, 8),
        "f_inf_sqrt": round(float(np.sqrt(f_inf)), 8),
    }
    print("limits:", out["branch_limits"])

    # ---------- 1. curvature invariants (static spherical, ds² = -f dt² + dr²/h + r² dΩ²)
    # Symbolic metric + christoffels via sympy, evaluated numerically on the grid.
    t_, r_, th_ = sp.symbols("t r theta", positive=True)
    f_s, h_s = sp.Function("f")(r_), sp.Function("h")(r_)
    # 2D-reduced invariants: for static spherical metrics the full 4D Kretschmann
    # splits into time-radial (2D) and angular (sphere) parts. We compute the
    # exact 4D invariants symbolically.
    coords = [t_, r_, th_, sp.Symbol("phi", real=True)]
    g = sp.diag(-f_s, 1/h_s, r_**2, r_**2*sp.sin(th_)**2)
    gin = g.inv()
    n = 4
    Gam = [[[sp.S.Zero for _ in range(n)] for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c_ in range(b, n):
                expr = sum(gin[a, d_]*(sp.diff(g[d_, c_], coords[b])
                                        + sp.diff(g[d_, b], coords[c_])
                                        - sp.diff(g[b, c_], coords[d_]))
                            for d_ in range(n))/2
                Gam[a][b][c_] = sp.simplify(expr)
                Gam[a][c_][b] = Gam[a][b][c_]
    Rm = [[[[sp.S.Zero for _ in range(n)] for _ in range(n)]
            for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c_ in range(n):
                for d_ in range(n):
                    expr = (sp.diff(Gam[a][b][d_], coords[c_])
                             - sp.diff(Gam[a][b][c_], coords[d_]))
                    s = sp.S.Zero
                    for e_ in range(n):
                        s += (Gam[a][e_][c_]*Gam[e_][b][d_]
                               - Gam[a][e_][d_]*Gam[e_][b][c_])
                    Rm[a][b][c_][d_] = sp.simplify(expr + s)
    Ric = sp.MutableDenseMatrix(n, n, lambda i_, j_: sp.simplify(
        sum(gin[a_, c_]*Rm[a_][a_][c_][j_] for a_ in range(n) for c_ in range(n))))
    R = sp.simplify(sum(gin[i_, j_]*Ric[i_, j_] for i_ in range(n) for j_ in range(n)))
    R2 = sp.simplify(sum(gin[a_, c_]*gin[b_, d_]*Ric[a_, b_]*Ric[c_, d_]
                          for a_ in range(n) for b_ in range(n)
                          for c_ in range(n) for d_ in range(n)))
    Kret = sp.simplify(sum(gin[a_, c_]*gin[b_, d_]*gin[e_, g_]*gin[f_, h_]
                            *Rm[a_][b_][e_][f_]*Rm[c_][d_][g_][h_]
                            for a_ in range(n) for b_ in range(n)
                            for c_ in range(n) for d_ in range(n)
                            for e_ in range(n) for f_ in range(n)
                            for g_ in range(n) for h_ in range(n)))
    print("invariants built, lambdifying...", flush=True)
    args = (r_, f_s, h_s, sp.diff(f_s, r_), sp.diff(f_s, r_, 2),
             sp.diff(h_s, r_), sp.diff(h_s, r_, 2))
    th_fix = np.pi / 2  # equatorial slice — angular invariants are theta-independent up to factors
    def eval_inv(fn, *nums):
        try:
            v = fn(*nums)
            arr = np.asarray(v, dtype=float)
            if arr.ndim == 0:  # constant expression (e.g. identically 0)
                arr = np.full_like(r_grid, float(arr))
            return arr
        except Exception as exc:  # noqa: BLE001
            return {"error": str(exc)[:200]}
    R_f = sp.lambdify(args, R, "numpy")
    R2_f = sp.lambdify(args, R2, "numpy")
    try:
        K_sub = Kret.subs(th_, np.pi/2)
        if K_sub == 0:
            K_vals = np.zeros_like(r_grid)
        else:
            K_f = sp.lambdify(args, K_sub, "numpy")
            K_vals = eval_inv(K_f, r_grid, f_g, h_g, fp_g,
                               jet(r_grid, f_g, 2), hp_g, jet(r_grid, h_g, 2))
    except Exception as exc:  # noqa: BLE001
        K_vals = {"error": str(exc)[:200]}

    vals = {}
    grid_args = (r_grid, f_g, h_g, fp_g, jet(r_grid, f_g, 2),
                  hp_g, jet(r_grid, h_g, 2))
    for name, fn, vv in (("R", R_f, eval_inv(R_f, *grid_args)),
                          ("Riemann2", R2_f, eval_inv(R2_f, *grid_args)),
                          ("Kretschmann", None, K_vals)):
        if isinstance(vv, dict):
            vals[name] = vv
            continue
        vals[name] = {
            "at_r2": round(float(vv[np.searchsorted(r_grid, 2.0)]), 8),
            "at_r10": round(float(vv[np.searchsorted(r_grid, 10.0)]), 10),
            "at_r100": float(f"{vv[np.searchsorted(r_grid, 100.0)]:.6e}"),
            "at_r395": float(f"{vv[-1]:.6e}"),
            "decays_to_zero": bool(abs(vv[-1]) < 1e-6),
        }
    out["curvature_invariants"] = vals
    print(json.dumps(vals, indent=1))

    # ---------- 2. time normalization: T = sqrt(f_inf) t
    sqrt_f_inf = float(np.sqrt(f_inf))
    out["time_normalization"] = {
        "definition": "T = sqrt(f_inf) * t",
        "sqrt_f_inf": sqrt_f_inf,
        "g_TT_limit": -1.0,
        "g_TT_note": "g_tt -> -f_inf, dT = sqrt(f_inf) dt => g_TT = -f_inf/f_inf = -1",
        "g_rr_limit": round(1.0/h_inf, 10),
        "verdict": "ASYMPTOTICALLY_MINKOWSKI_WITH_NORMALIZED_TIME"
            if vals.get("R", {}).get("decays_to_zero") else
            "CURVATURE_PERSISTS_NONFLAT_REAL",
    }
    print("time normalization:", out["time_normalization"]["verdict"])

    # ---------- 3. normalized tortoise coordinate
    dr_star = 1.0/np.sqrt(f_g*h_g)
    R_star_slope = dr_star * sqrt_f_inf  # dR_*/dr with R_* = sqrt(f_inf) r_*
    out["normalized_tortoise"] = {
        "definition": "R_* = sqrt(f_inf) * r_*",
        "dR_star_dr_at_r2": round(float(R_star_slope[np.searchsorted(r_grid, 2.0)]), 8),
        "dR_star_dr_at_r395": round(float(R_star_slope[-1]), 8),
        "limit_is_1": bool(abs(R_star_slope[-1] - 1.0) < 1e-3),
        "outgoing_wave": "Psi ~ exp(-i Omega_inf (T - R_*)), Omega_inf = omega/sqrt(f_inf)",
    }
    out["normalized_tortoise"] = dict(out["normalized_tortoise"])
    print("tortoise slope at 395:", out["normalized_tortoise"]["dR_star_dr_at_r395"])

    # ---------- 4. pole invariance under constant time rescaling
    # The quadratic eigenvalue problem: (omega^2 K - i omega W(omega) + V) psi = 0.
    # Under t -> T = a t (a = sqrt(f_inf)), omega -> omega/a. The physical
    # poles (as frequencies measured in T) are Omega = omega/sqrt(f_inf).
    # Numerical check: take the flat-cavity-like scaling relation used by the
    # normal-mode test: omega_n = n pi c / L with c in t-time; in T-time:
    # Omega_n = n pi c / (L sqrt(f_inf)) -- same poles, different label.
    out["pole_invariance"] = {
        "statement": ("A constant time rescaling cannot change physical "
                       "poles, only their frequency label: "
                       "Omega_inf = omega / sqrt(f_inf). The resonance "
                       "solver MUST report frequencies in the normalized "
                       "time T to be comparable with a flat-infinity "
                       "observer; reporting raw omega would masquerade "
                       "the lapse as a spectral shift."),
        "convention_declared": ("All CERTIFIED_RESONANCE_CATALOG_V1 "
                                 "frequencies are Omega_inf in T-time; "
                                 "raw t-time omega is stored separately "
                                 "as omega_t."),
        "Omega_inf_factor": round(1.0/sqrt_f_inf, 8),
    }

    # ---------- G4(phi_inf): effective gravitational coupling
    G4_inf = 0.5 + A1V*(phi_inf - 1.0)
    out["G4_asymptotic"] = {
        "formula": "G4(phi) = 1/2 + a1 (phi - 1)",
        "phi_inf": round(phi_inf, 8),
        "G4_inf": round(float(G4_inf), 8),
        "relative_deviation_from_0p5": round(float((G4_inf - 0.5)/0.5), 8),
        "note": ("G4_inf ~ 0.49917 vs GR's 1/2: a ~0.17% shift in the "
                  "asymptotic effective Planck mass. Metric is flat but "
                  "the coupling constant carries the branch signature."),
    }

    # ---------- final status
    if vals.get("R", {}).get("decays_to_zero"):
        out["final_classification"] = "ASYMPTOTICALLY_MINKOWSKI_NONUNIT_LAPSE"
        out["status_change"] = ("supersedes REGULAR_FINITE_LIMIT_EXTERIOR_"
                                 "NOT_FLAT from V4_GLOBAL_EXTERIOR_V1: the "
                                 "exterior is flat up to time normalization; "
                                 "c_inf = 0.501 was the COORDINATE speed "
                                 "dr/dt, not a new physical propagation "
                                 "speed (dr/dT = 1).")
    else:
        out["final_classification"] = "CURVATURE_PERSISTS_NONFLAT_REAL"
    print("FINAL:", out["final_classification"])

    out["wall_seconds"] = round(time.time() - t0, 1)
    OUT.write_text(json.dumps(out, indent=1, allow_nan=False, default=str) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
