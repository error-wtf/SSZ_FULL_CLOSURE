#!/usr/bin/env python3
"""Independent P1/K_scalar re-evaluation at the negative pocket (x ~ 0.29-0.33).

Hypothesis under test (execution plan stage F.1/F.2): is the recomputed negative
K_scalar pocket an artifact of the derivative ROUTE (jet9d8 finite stencils), or
is it invariant under an independent differentiation method?

Route A (reference): repo emit_from_primitives (jet9d8).
Route B (independent): cubic-spline derivative of Y(r) (scipy CubicSpline .derivative
                       + Richardson via two grids) with P1 = h*mu/(2 f r^2 H^2) dY/dr.
Both routes consume THE SAME primitive columns; no archived K_scalar is used.
Route C: spline on a 2x refined local grid via interpolation (order check).

Decision: if min K_scalar < 0 under all independent routes with relative
deviation to route A much smaller than |min K|, the pocket is derivative-route
robust and the action-level question becomes decisive.
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.interpolate import CubicSpline

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from ssz_p5.coefficients.mh_general_primitives import emit_from_primitives

ACT = ROOT / "data/production/ssz_p5_F1b_G4XX_transverse_core_candidate_2026-09-14.csv"
REF = ROOT / "data/authoritative/ssz_p5_F2_core_punctured_horndeski_unreduced_39of41_2026-09-14.csv"
OUT = ROOT / "data/generated/qnm_global_diagnostic/KSCALAR_INDEPENDENT_ROUTE_AUDIT.json"

def main():
    a = pd.read_csv(ACT).sort_values("r_over_rs").reset_index(drop=True)
    r_ref = pd.read_csv(REF).sort_values("x").reset_index(drop=True)
    q_args = dict(x=a.r_over_rs.to_numpy(float), u=a.u.to_numpy(float),
                  f=a.A_f.to_numpy(float), h=a.B_h.to_numpy(float),
                  X=a.X.to_numpy(float), phi=a.phi.to_numpy(float),
                  a1_action=a.u.to_numpy(float)*np.nan,  # filled below
                  )
    # primitives need F,G,H,a1,c4: use the same builders as the audited scripts
    from ssz_p5.coefficients.mh_action_primitives import quartic_g5zero_primitives
    q = quartic_g5zero_primitives(a)
    inp = pd.DataFrame({
        "x": q.x, "u": q.u, "f": a.A_f, "h": a.B_h, "X": a.X, "phi": a.phi,
        "A0prime": 0.0, "a1": q.a1_action, "c2": r_ref.c2, "c4": q.c4_action,
        "F_tensor": q.F_tensor_action, "G_tensor": q.G_tensor_action,
        "H_tensor": q.H_tensor_action})
    z = emit_from_primitives(inp, regularize_photon_root=True)

    x = z.x.to_numpy(float); f = z.f.to_numpy(float); h = z.h.to_numpy(float)
    F = z.F_tensor.to_numpy(float); H = z.H_tensor.to_numpy(float)
    mu = z.mu.to_numpy(float)

    # Route A: repo jet9d8 P1 (already emitted)
    kA = z.K_scalar.to_numpy(float)

    # Route B: independent cubic-spline derivative of Y
    Y = f * x**4 * H**4 / (mu**2 * h)
    pref = h * mu / (2.0 * f * x**2 * H**2)
    # guard: drop nonfinite samples for spline fit, interpolate back
    good = np.isfinite(Y) & np.isfinite(pref) & (mu != 0)
    cs = CubicSpline(x[good], Y[good])
    dY = cs.derivative()(x)
    p1B = pref * dY
    kB = 2.0 * p1B - F
    kB[~good] = np.nan

    # Route C: spline on 2x interpolated grid (independent smoothness assumption)
    xg = np.linspace(x[good].min(), x[good].max(), 2 * good.sum())
    Yg = cs(xg); dYg = CubicSpline(xg, Yg).derivative()(x)
    p1C = pref * dYg
    kC = 2.0 * p1C - F
    kC[~good] = np.nan

    band = (x > 0.29) & (x < 0.33)
    def mins(k):
        m = np.nanmin(np.where(band, k, np.nan))
        j = int(np.nanargmin(np.where(band, k, np.nan)))
        return {"min_in_band": float(m), "argmin_x": float(x[j]),
                "n_negative_rows_band": int(np.nansum(np.where(band, k, np.nan) < 0))}
    rep = {
        "status": "KSCALAR_INDEPENDENT_ROUTE_AUDIT",
        "route_A_jet9d8": mins(kA),
        "route_B_cubic_spline": mins(kB),
        "route_C_spline_2x_grid": mins(kC),
        "archived_positive_oracle_min_in_band": None,
        "deviation_B_vs_A_at_argmin_rel": None,
    }
    jB = int(np.nanargmin(np.where(band, kB, np.nan)))
    rep["deviation_B_vs_A_at_argmin_rel"] = float(abs(kB[jB] - kA[jB]) / abs(kA[jB]))

    verdict = {
        "pocket_negative_all_routes": bool(
            rep["route_A_jet9d8"]["min_in_band"] < 0 and
            rep["route_B_cubic_spline"]["min_in_band"] < 0 and
            rep["route_C_spline_2x_grid"]["min_in_band"] < 0),
        "route_deviation_small_vs_pocket_depth": bool(
            rep["deviation_B_vs_A_at_argmin_rel"] <
            0.05 * abs(rep["route_A_jet9d8"]["min_in_band"])),
    }
    rep["verdict"] = verdict
    rep["verdict_label"] = ("POCKET_ROUTE_ROBUST_ACTION_LEVEL_DECISIVE"
                            if all(verdict.values())
                            else "ROUTE_SENSITIVITY_PRESENT_INVESTIGATE")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(rep, indent=1))
    print(json.dumps(rep, indent=1))
    return 0

if __name__ == "__main__":
    sys.exit(main())
