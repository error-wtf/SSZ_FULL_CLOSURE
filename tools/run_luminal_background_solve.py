#!/usr/bin/env python3
"""N2: luminal background solve via least-squares relaxation.

G4 = 1/2 constant (declared minimal-consistency law; the archived
profile is numerically ~0.5).  Solves for (f, h, phi)(u) such that the
luminal constraints hold: G2 is eliminated analytically from E00 (it
is linear), and E00-G2 == E11-G2 (identical G2 in both equations)
becomes the imposed residual.  Outer asymptotics anchored
(f -> 1/4, h -> 1, phi -> 1; large-u end in the repo convention).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / ("data/production/"
                  "ssz_p5_F1b_G4XX_transverse_core_candidate_2026-09-14.csv")

G4 = 0.5
G4PHI = 0.0
G2X = 1.0


def main() -> int:
    d = pd.read_csv(ARCHIVE).sort_values("u")
    u = d["u"].to_numpy(float)
    n = len(u)
    coarse = np.unique(np.round(np.linspace(0, n - 1, 200)).astype(int))
    uc = u[coarse]
    y0 = np.empty(3 * len(uc))
    y0[0::3] = d["A_f"].to_numpy(float)[coarse]
    y0[1::3] = d["B_h"].to_numpy(float)[coarse]
    y0[2::3] = d["phi"].to_numpy(float)[coarse]

    def derivs(y):
        f = y[0::3]
        h = y[1::3]
        phi = y[2::3]
        return (f, h, phi, np.gradient(f, uc), np.gradient(h, uc),
                np.gradient(phi, uc), np.gradient(np.gradient(phi, uc), uc))

    def residuals(y):
        f, h, phi, fp, hp, phip, phipp = derivs(y)
        r = 1.0 / uc
        e11 = (-G2X * h * phi**2 + 2 * G4 * h / r**2 - 2 * G4 / r**2
               + 2 * G4 * fp * h / (f * r) + 4 * G4PHI * h * phi / r
               + G4PHI * fp * h * phi / f)
        g2_e00 = (2 * G4 * h / r**2 + 2 * G4 * hp / r - 2 * G4 / r**2
                  - 2 * G4PHI * h * phi**2 - 4 * G4PHI * h * phi / r
                  - 2 * G4PHI * h * phipp - G4PHI * hp * phi)
        cons = g2_e00 - e11
        w = 10.0
        return np.concatenate([
            cons,
            w * (f[-3:] - 0.25),
            w * (h[-3:] - 1.0),
            w * (phi[-3:] - 1.0),
        ])

    sol = least_squares(residuals, y0, method="trf", max_nfev=20000)
    _, _, _, _, _, _, _ = derivs(sol.x)
    res = residuals(sol.x)

    out = {
        "audit": "N2_LUMINAL_SOLVE_V1",
        "g4_law": "constant 1/2 (declared)",
        "grid_points": int(len(uc)),
        "u_range": [float(uc[0]), float(uc[-1])],
        "cost": float(sol.cost),
        "optimality": float(sol.optimality),
        "status": int(sol.status),
        "max_abs_residual": float(np.max(np.abs(res))),
        "median_abs_residual": float(np.median(np.abs(res))),
        "note": ("G2 eliminated analytically from E00 (linear); the "
                 "E00-G2 == E11-G2 consistency is the imposed constraint. "
                 "Boundary anchors weighted 10x.  G4 = 1/2 constant."),
    }
    print(json.dumps({k: v for k, v in out.items() if k != "note"}, indent=1))
    out_path = ROOT / "data/generated/spectral/N2_LUMINAL_SOLVE_RESULT.json"
    out_path.write_text(json.dumps(out, indent=1) + "\n")
    print(f"written: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
