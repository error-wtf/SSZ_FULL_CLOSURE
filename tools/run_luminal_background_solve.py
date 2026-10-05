#!/usr/bin/env python3
"""N2 V2: luminal background solve — full G2-law version.

G4 = 1/2 constant, G4phi = 0 (declared minimal-consistency law).
G2(X) = G2X * X (declared; G2(0)=0 is the consistency choice — verified:
with G2(0)=-1 the flat-end residuals stay O(1e-3), with G2(0)=0 they
reach O(1e-8)).

Imposed residuals: E00 = 0 AND E11 = 0 with the complete luminal
expressions (G2 through the law, not eliminated).  Outer asymptotics
anchored (f -> 1/4, h -> 1, phi -> 1; large-u end).
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
G2C = 0.0


def main() -> int:
    d = pd.read_csv(ARCHIVE).sort_values("u")
    u = d["u"].to_numpy(float)
    coarse = np.unique(np.round(np.linspace(0, len(u) - 1, 200)).astype(int))
    uc = u[coarse]
    y0 = np.empty(3 * len(uc))
    y0[0::3] = d["A_f"].to_numpy(float)[coarse]
    y0[1::3] = d["B_h"].to_numpy(float)[coarse]
    y0[2::3] = d["phi"].to_numpy(float)[coarse]

    def derivs(y):
        f, h, phi = y[0::3], y[1::3], y[2::3]
        return (f, h, phi, np.gradient(f, uc), np.gradient(h, uc),
                np.gradient(phi, uc),
                np.gradient(np.gradient(phi, uc), uc))

    def residuals(y):
        f, h, phi, fp, hp, phip, phipp = derivs(y)
        r = 1.0 / uc
        X = phip**2 / (2.0 * np.maximum(f, 1e-12))
        G2 = G2X * X + G2C
        e00 = (G2 - 2 * G4 * h / r**2 - 2 * G4 * hp / r + 2 * G4 / r**2
               - 2 * G4PHI * h * phi**2 - 4 * G4PHI * h * phi / r
               - 2 * G4PHI * h * phipp - G4PHI * hp * phi)
        e11 = (-G2 + 2 * G4 * h / r**2 - 2 * G4 / r**2 + 2 * G4 * fp * h / (f * r)
               + 4 * G4PHI * h * phi / r + G4PHI * fp * h * phi / f)
        w = 10.0
        return np.concatenate([
            e00, e11,
            w * (f[-3:] - 0.25), w * (h[-3:] - 1.0), w * (phi[-3:] - 1.0),
        ])

    sol = least_squares(residuals, y0, method="trf", max_nfev=20000)
    f, h, phi, fp, hp, phip, phipp = derivs(sol.x)
    res = residuals(sol.x)

    out = {
        "audit": "N2_LUMINAL_SOLVE_V2",
        "g4_law": "constant 1/2 (declared)",
        "g2_law": "G2 = G2X*X + G2C with G2X=1, G2C=0 (declared; verified "
                  "consistency choice)",
        "grid_points": int(len(uc)),
        "u_range": [float(uc[0]), float(uc[-1])],
        "cost": float(sol.cost),
        "optimality": float(sol.optimality),
        "status": int(sol.status),
        "max_abs_residual": float(np.max(np.abs(res))),
        "median_abs_residual": float(np.median(np.abs(res))),
        "profile_check": {
            "f_range": [float(f.min()), float(f.max())],
            "h_range": [float(h.min()), float(h.max())],
            "phi_range": [float(phi.min()), float(phi.max())],
        },
        "next": "N3: quartic_g5zero_primitives on THIS solution, then N4 "
                "decision against the declared rules.",
    }
    print(json.dumps({k: v for k, v in out.items()
                      if k not in ("note", "next")}, indent=1))
    out_path = ROOT / "data/generated/spectral/N2_LUMINAL_SOLVE_RESULT_V2.json"
    out_path.write_text(json.dumps(out, indent=1) + "\n")
    print(f"written: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
