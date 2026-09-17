"""Authoritative regional cover for the full-SVT production member.

The strong-H carrier is a principal/control witness only.  It is never used to
classify the central genuine-SVT lobe.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

CENTRAL = "data/production/ssz_p5_integrable_svt_unreduced_even_coefficients_2026-09-12.csv"


def classify_u(u: float) -> str:
    if u < 0.5515230871346237:
        return "weak_exterior_H"
    if u < 0.61:
        return "outer_same_action_H_SVT"
    if u < 0.71:
        return "central_exact_SVT"
    if u < 0.715:
        return "inner_same_action_SVT_H"
    return "punctured_H_core"


def audit_light_ring_regions(root: Path) -> dict:
    d = pd.read_csv(Path(root) / CENTRAL)
    out = []
    for x in (1.5999977950251312, 1.41616064):
        row = d.iloc[int(np.argmin(np.abs(d.x.to_numpy(float) - x)))]
        out.append({
            "x_requested": x,
            "x": float(row.x),
            "u": float(row.u),
            "region": classify_u(float(row.u)),
            "A0prime": float(row.A0prime),
            "pure_H_null_vector_gate": "NOT_APPLICABLE",
        })
    return {
        "cover": {
            "weak_exterior_H": "u < 0.5515230871346237",
            "outer_same_action_H_SVT": "0.5515230871346237 <= u < 0.61",
            "central_exact_SVT": "0.61 <= u < 0.71",
            "inner_same_action_SVT_H": "0.71 <= u < 0.715",
            "punctured_H_core": "u >= 0.715",
        },
        "witnesses": out,
        "strong_H_carrier_role": "principal/control witness; never global production background",
    }
