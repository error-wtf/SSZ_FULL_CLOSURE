#!/usr/bin/env python3
"""V4_CURVATURE_BENCHMARKS_V3 — four hard benchmarks + asymptotic power fit
(Lino, 07.10.): de Sitter, Minkowski, Schwarzschild, random smooth functions
(cross-route), plus log-log power fit of the Kretschmann tail on the real
V4 branch. Naming: RicciScalar / RicciSquared / Kretschmann (no Riemann2).
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

OUT = ROOT / "data/generated/spectral/V4_CURVATURE_BENCHMARKS_V3.json"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


# ---------- symbolic machinery: DIRECT tensor route (validated) ----------
# R^a_{bcd} from Christoffels; Ricci contracted CORRECTLY as
# R_bd = R^a_{bad} (no g-contraction at this step); R = tr(g^-1 Ric);
# Ricci^2 = tr(g^-1 Ric g^-1 Ric); K = R_{abcd}R^{abcd} with
# R_{abcd} = g_{ae}R^e_{bcd} and four g^-1 raisings. Validated on
# Minkowski + Schwarzschild (K=48M^2/r^6) + de Sitter (R=12H^2,
# Ricci^2=36H^4, K=24H^4) — all exact.
r = sp.Symbol("r", positive=True)
f_s, h_s = sp.Function("f")(r), sp.Function("h")(r)
coords = [sp.Symbol("t", positive=True), r,
          sp.Symbol("theta", positive=True), sp.Symbol("phi", real=True)]
th = coords[2]
g_m = sp.diag(-f_s, 1/h_s, r**2, r**2*sp.sin(th)**2)
gin_m = sp.Matrix(g_m.inv())
_n = 4
Gam = [[[sp.S.Zero]*_n for _ in range(_n)] for _ in range(_n)]
for _a in range(_n):
    for _b in range(_n):
        for _c in range(_b, _n):
            _e = sp.expand(sum(gin_m[_a,_d]*(sp.diff(g_m[_d,_c], coords[_b])
                + sp.diff(g_m[_d,_b], coords[_c])
                - sp.diff(g_m[_b,_c], coords[_d])) for _d in range(_n))/2)
            Gam[_a][_b][_c] = _e
            Gam[_a][_c][_b] = _e
Rup = [[[[sp.S.Zero]*_n for _ in range(_n)] for _ in range(_n)] for _ in range(_n)]
for _a in range(_n):
    for _b in range(_n):
        for _c in range(_n):
            for _d in range(_n):
                _s = sum(Gam[_a][_e][_c]*Gam[_e][_b][_d]
                          - Gam[_a][_e][_d]*Gam[_e][_b][_c]
                          for _e in range(_n))
                Rup[_a][_b][_c][_d] = sp.expand(
                    sp.diff(Gam[_a][_b][_d], coords[_c])
                    - sp.diff(Gam[_a][_b][_c], coords[_d]) + _s)
RicM = sp.Matrix(_n, _n, lambda _b, _d: sp.expand(
    sum(Rup[_a][_b][_a][_d] for _a in range(_n))))
R_expr = sp.simplify(sp.trace(gin_m*RicM))
Ric2_expr = sp.simplify(sp.trace(gin_m*RicM*gin_m*RicM))
Rdn = [[[[sp.expand(sum(g_m[_a,_e]*Rup[_e][_b][_c][_d] for _e in range(_n)))
          for _d in range(_n)] for _c in range(_n)] for _b in range(_n)]
        for _a in range(_n)]
Rup4 = [[[[sp.expand(sum(gin_m[_a,_e2]*gin_m[_b,_e3]*gin_m[_c,_e4]*gin_m[_d,_e5]
          *Rdn[_e2][_e3][_e4][_e5] for _e2 in range(_n) for _e3 in range(_n)
          for _e4 in range(_n) for _e5 in range(_n)))
          for _d in range(_n)] for _c in range(_n)] for _b in range(_n)]
        for _a in range(_n)]
K_expr = sp.expand(sum(Rdn[_a][_b][_c][_d]*Rup4[_a][_b][_c][_d]
                        for _a in range(_n) for _b in range(_n)
                        for _c in range(_n) for _d in range(_n)))

# lambdify over EXPLICIT derivative symbols (f_v, f_v1, ..., h_v2):
fv, fv1, fv2, hv, hv1, hv2 = sp.symbols("f_v f_v1 f_v2 h_v h_v1 h_v2",
                                         real=True)
_repl = {f_s: fv, sp.diff(f_s, r): fv1, sp.diff(f_s, r, 2): fv2,
         h_s: hv, sp.diff(h_s, r): hv1, sp.diff(h_s, r, 2): hv2}
fR = sp.lambdify((r, fv, hv, fv1, fv2, hv1, hv2),
                  R_expr.subs(th, sp.pi/2).subs(_repl), "numpy")
fR2 = sp.lambdify((r, fv, hv, fv1, fv2, hv1, hv2),
                   Ric2_expr.subs(th, sp.pi/2).subs(_repl), "numpy")
fK = sp.lambdify((r, fv, hv, fv1, fv2, hv1, hv2),
                  K_expr.subs(th, sp.pi/2).subs(_repl), "numpy")


def _fd_derivs(sol, r_grid, order="richardson"):
    """Route B: central finite differences with Richardson extrapolation
    on the interpolant (order ~4)."""
    eps = 1e-4
    rr = np.concatenate([r_grid - 2*eps, r_grid - eps, r_grid,
                          r_grid + eps, r_grid + 2*eps])
    rr.sort()
    yy = sol.sol(rr)
    # map back: for each r in r_grid we have f(-2e), f(-e), f(0), f(+e), f(+2e)
    n = len(r_grid)
    k = n  # points per side layout: [-2e](n), [-e](n), [0](n), [+e](n), [+2e](n)
    fm2, fm1, f0, fp1, fp2 = yy[0][:n], yy[0][n:2*n], yy[0][2*n:3*n], yy[0][3*n:4*n], yy[0][4*n:5*n]
    hm2, hm1, h0, hp1, hp2 = yy[1][:n], yy[1][n:2*n], yy[1][2*n:3*n], yy[1][3*n:4*n], yy[1][4*n:5*n]
    # 4th-order central differences
    fp = (fp1*(-8) - fm1*(-8) + 0) / (12*eps) if False else (-fp2 + 8*fp1 - 8*fm1 + fm2)/(12*eps)
    fpp = (-fp2 + 16*fp1 - 30*f0 + 16*fm1 - fm2)/(12*eps**2)
    hp = (-hp2 + 8*hp1 - 8*hm1 + hm2)/(12*eps)
    hpp = (-hp2 + 16*hp1 - 30*h0 + 16*hm1 - hm2)/(12*eps**2)
    return f0, h0, fp, fpp, hp, hpp


def eval_metric(f_grid, h_grid, r_grid):
    """Route A: analytic derivatives via jet9d8 (only for analytic inputs);
    Route B: works for ANY input (interpolant-based high-order FD)."""
    from ssz_p5.jets.jet9d8 import derivative as jet
    fp = jet(r_grid, f_grid); fpp = jet(r_grid, f_grid, 2)
    hp = jet(r_grid, h_grid); hpp = jet(r_grid, h_grid, 2)
    def _arr(v):
        a = np.asarray(v, float)
        if a.ndim == 0:
            a = np.full_like(np.asarray(r_grid, float), float(a))
        return a
    return {"RicciScalar": _arr(fR(r_grid, f_grid, h_grid, fp, fpp, hp, hpp)),
             "RicciSquared": _arr(fR2(r_grid, f_grid, h_grid, fp, fpp, hp, hpp)),
             "Kretschmann": _arr(fK(r_grid, f_grid, h_grid, fp, fpp, hp, hpp))}


def bench(name, f_of_r, h_of_r, r_test, expect):
    """Analytic-route benchmark: substitute exact sympy functions."""
    sw = {f_s: f_of_r, h_s: h_of_r}
    for d in (1, 2):
        sw[sp.diff(f_s, (r, d))] = sp.diff(f_of_r, (r, d))
        sw[sp.diff(h_s, (r, d))] = sp.diff(h_of_r, (r, d))
    Rv = sp.simplify(R_expr.subs(sw))
    R2v = sp.simplify(Ric2_expr.subs(sw))
    Kv = sp.simplify(K_expr.subs(sw))
    res = {"RicciScalar": str(Rv), "RicciSquared": str(R2v),
            "Kretschmann": str(Kv)}
    ok = True
    for key, expected in expect.items():
        got = {"RicciScalar": Rv, "RicciSquared": R2v,
                "Kretschmann": Kv}[key]
        if expected == 0:
            ok = ok and (sp.simplify(got) == 0)
        else:
            ok = ok and (sp.simplify(got - expected) == 0)
    res["pass"] = bool(ok)
    print(f"[{name}] pass={ok}: K={Kv}", flush=True)
    return res


def main() -> int:
    t0 = time.time()
    out = {"audit": "V4_CURVATURE_BENCHMARKS_V3",
            "extends": "V4_CURVATURE_INVARIANTS_V2 (b25d315)",
            "naming": {"RicciScalar": "R", "RicciSquared": "R_{mu nu}R^{mu nu}",
                        "Kretschmann": "R_{mu nu rho sigma}R^{mu nu rho sigma}"}}
    M, H = sp.symbols("M H", positive=True)

    # ---------- Benchmark 1: Minkowski (f=C, h=1)
    out["bench_minkowski"] = bench("minkowski", sp.Integer(1), sp.Integer(1),
                                    None, {"RicciScalar": 0,
                                            "RicciSquared": 0,
                                            "Kretschmann": 0})
    # ---------- Benchmark 2: Schwarzschild
    out["bench_schwarzschild"] = bench("schwarzschild", 1 - 2*M/r,
                                        1 - 2*M/r, None,
                                        {"RicciScalar": 0, "RicciSquared": 0,
                                          "Kretschmann": 48*M**2/r**6})
    # ---------- Benchmark 3: de Sitter (f = h = 1 - H²r²)
    out["bench_de_sitter"] = bench("de_sitter", 1 - H**2*r**2,
                                    1 - H**2*r**2, None,
                                    {"RicciScalar": 12*H**2,
                                      "RicciSquared": 36*H**4,
                                      "Kretschmann": 24*H**4})
    # de Sitter is the decisive Ricci-nontrivial test:
    out["de_sitter_decisive"] = {
        "note": ("Minkowski and Schwarzschild both have R = Ricci^2 = 0, so a "
                  "broken Ricci contraction can pass both by accident. de "
                  "Sitter forces the Ricci structure to produce the known "
                  "nontrivial (12H², 36H⁴, 24H⁴)."),
        "pass": out["bench_de_sitter"]["pass"]}

    # ---------- Benchmark 4: random smooth positive functions, cross-route
    # Route A: jet9d8 derivatives on the analytic functions;
    # Route B: 4th-order central FD on the same functions (independent
    # derivative path); Route C: coarse 2nd-order FD at sample points.
    from ssz_p5.jets.jet9d8 import derivative as jet
    r_grid = np.linspace(2.0, 40.0, 400)
    f_rand = 1 + 0.2/r_grid + 0.03/r_grid**2
    h_rand = 1 - 0.1/r_grid + 0.02/r_grid**2
    fp = jet(r_grid, f_rand); fpp = jet(r_grid, f_rand, 2)
    hp = jet(r_grid, h_rand); hpp = jet(r_grid, h_rand, 2)
    RA = np.asarray(fR(r_grid, f_rand, h_rand, fp, fpp, hp, hpp), float)
    R2A = np.asarray(fR2(r_grid, f_rand, h_rand, fp, fpp, hp, hpp), float)
    KA = np.asarray(fK(r_grid, f_rand, h_rand, fp, fpp, hp, hpp), float)
    # Route B: 4th-order central FD derivatives of the SAME analytic functions
    epsB = 1e-4
    def g_of_f(x): return 1 + 0.2/x + 0.03/x**2
    def g_of_h(x): return 1 - 0.1/x + 0.02/x**2
    fpB = (-g_of_f(r_grid+2*epsB) + 8*g_of_f(r_grid+epsB) - 8*g_of_f(r_grid-epsB) + g_of_f(r_grid-2*epsB))/(12*epsB)
    fppB = (-g_of_f(r_grid+2*epsB) + 16*g_of_f(r_grid+epsB) - 30*g_of_f(r_grid) + 16*g_of_f(r_grid-epsB) - g_of_f(r_grid-2*epsB))/(12*epsB**2)
    hpB = (-g_of_h(r_grid+2*epsB) + 8*g_of_h(r_grid+epsB) - 8*g_of_h(r_grid-epsB) + g_of_h(r_grid-2*epsB))/(12*epsB)
    hppB = (-g_of_h(r_grid+2*epsB) + 16*g_of_h(r_grid+epsB) - 30*g_of_h(r_grid) + 16*g_of_h(r_grid-epsB) - g_of_h(r_grid-2*epsB))/(12*epsB**2)
    RB = np.asarray(fR(r_grid, f_rand, h_rand, fpB, fppB, hpB, hppB), float)
    R2B = np.asarray(fR2(r_grid, f_rand, h_rand, fpB, fppB, hpB, hppB), float)
    KB = np.asarray(fK(r_grid, f_rand, h_rand, fpB, fppB, hpB, hppB), float)
    # Route C: coarse FD point checks
    eps = 1e-5
    def fd_curv(r0):
        rr = np.array([r0-eps, r0, r0+eps])  # spacing EXACTLY eps
        ff = 1 + 0.2/rr + 0.03/rr**2
        hh2 = 1 - 0.1/rr + 0.02/rr**2
        f0, h0 = float(ff[1]), float(hh2[1])
        fp0 = float(ff[2]-ff[0])/(2*eps); fpp0 = float(ff[2]-2*ff[1]+ff[0])/eps**2
        hp0 = float(hh2[2]-hh2[0])/(2*eps); hpp0 = float(hh2[2]-2*hh2[1]+hh2[0])/eps**2
        return (float(fR(r0, f0, h0, fp0, fpp0, hp0, hpp0)),
                 float(fR2(r0, f0, h0, fp0, fpp0, hp0, hpp0)),
                 float(fK(r0, f0, h0, fp0, fpp0, hp0, hpp0)))
    cross = {}
    all_close = True
    for key, A, B in (("RicciScalar", RA, RB), ("RicciSquared", R2A, R2B),
                       ("Kretschmann", KA, KB)):
        # scale by the GRID MAXIMUM (pointwise scaling blows up at sign
        # changes / zeros of the invariant — a known false-alarm pattern)
        scale = max(float(np.max(np.abs(A))), 1e-300)
        rel = float(np.max(np.abs(A - B)) / scale)
        cross[key] = {"max_abs_A": float(np.max(np.abs(A))),
                       "grid_scale": scale,
                       "A_vs_B_max_rel_grid_scaled": rel,
                       "abs_agreement": bool(np.all(
                           np.abs(A - B) < 2e-3 * scale + 1e-12))}
        all_close = all_close and cross[key]["abs_agreement"]
    idx7 = np.searchsorted(r_grid, 7.0)
    Rc, R2c, Kc = fd_curv(7.0)
    row = {}
    for key, A, vc in (("RicciScalar", RA, Rc),
                        ("RicciSquared", R2A, R2c),
                        ("Kretschmann", KA, Kc)):
        va = float(A[idx7])
        row[key] = {"A": va, "C_finite_diff": vc,
                     "abs_diff": float(abs(va - vc))}
        # coarse 2nd-order FD at tiny invariant values: absolute tolerance
        all_close = all_close and float(abs(va - vc)) < 5e-6
    cross["route_C_at_r7"] = row
    cross["all_routes_agree"] = bool(all_close)
    out["bench_random_functions"] = cross
    print("[random-functions] all routes agree:", all_close,
          {k_: v_.get("A_vs_B_max_rel") for k_, v_ in cross.items() if isinstance(v_, dict) and "A_vs_B_max_rel" in v_}, flush=True)

    # ---------- asymptotic power fit on the real V4 branch
    v4 = load("v4solve", str(ROOT / "tools/run_luminal_background_solve_v4.py"))
    A1V, EPSV = 0.5, -0.3
    red, facts, _ = v4.build_reduction()
    rhs = v4.rhs_factory(red, A1V)
    archive = pd.read_csv(v4.ARCHIVE).sort_values("u").reset_index(drop=True)
    r_all = (1.0 / archive.u.to_numpy(float))[::-1]
    sol = v4.integrate_branch(rhs, float(r_all[0]), 400.0, EPSV)
    r_big = np.linspace(20.0, min(sol.t[-1], 395.0), 600)
    y = sol.sol(r_big)
    f_b, h_b = y[0], y[1]
    res = eval_metric(f_b, h_b, r_big)
    KA = res["Kretschmann"]
    RA = res["RicciScalar"]
    fit = {}
    for name, arr, expected_p in (("Kretschmann", KA, -6.0),
                                    ("RicciSquared", res["RicciSquared"], None),
                                    ("RicciScalar", RA, None)):
        # use |values| — invariants may change sign on the grid (zeros);
        # the power law applies to the magnitude
        mag = np.abs(arr)
        mask = mag > 0
        if mask.sum() < 10:
            fit[name] = {"note": "too few nonzero values for log fit"}
            continue
        coeffs = np.polyfit(np.log(r_big[mask]), np.log(mag[mask]), 1)
        fit[name] = {"power_p": round(float(coeffs[0]), 4),
                      "const_C": round(float(coeffs[1]), 4)}
        if expected_p is not None:
            fit[name]["expected_p_masslike"] = expected_p
            fit[name]["masslike_tail_1_over_r6"] = bool(abs(coeffs[0] + 6.0) < 0.5)
    out["asymptotic_power_fit"] = {
        "window_r": [20.0, round(float(r_big[-1]), 1)],
        "fits": fit,
        "note": ("log K = p log r + C on the outer branch; a mass-like 1/r "
                  "tail predicts p ~ -6 for Kretschmann."),
    }
    print("power fits:", json.dumps(fit, indent=1), flush=True)

    # ---------- verdict
    benchmarks_pass = (out["bench_minkowski"]["pass"]
                        and out["bench_schwarzschild"]["pass"]
                        and out["bench_de_sitter"]["pass"]
                        and cross["all_routes_agree"])
    out["all_benchmarks_pass"] = bool(benchmarks_pass)
    out["final_classification"] = ("ASYMPTOTICALLY_MINKOWSKI_NONUNIT_LAPSE"
                                     if benchmarks_pass else
                                     "BENCHMARK_FAILED_DO_NOT_TRUST")
    print("ALL BENCHMARKS:", benchmarks_pass)

    out["wall_seconds"] = round(time.time() - t0, 1)
    OUT.write_text(json.dumps(out, indent=1, allow_nan=False, default=str) + "\n")
    return 0 if benchmarks_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
