#!/usr/bin/env python3
"""Audit the zero-vector branch at the P5 light rings and emit an electric seed."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from ssz_p5.action.light_ring_electric import light_ring_witnesses, required_electric_q_profile
from ssz_p5.jets.jet9d8 import derivative
from ssz_p5.provenance.manifest import sha256

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/production/ssz_p5_horndeski_carrier_through_light_rings_to_core_2026-09-12.csv"
OUTDIR = ROOT / "data/generated/absolute_attempt_2026-09-19"
REPORT = OUTDIR / "ZERO_VECTOR_LIGHT_RING_AUDIT.json"
SEED = OUTDIR / "ELECTRIC_HYBRID_EQ85_SEED.csv"


def main() -> int:
    d = pd.read_csv(SOURCE)
    r = d["r_over_rs"].to_numpy(float)
    f = d["f"].to_numpy(float)
    h = d["h"].to_numpy(float)
    H = d["H_equals_F_equals_G"].to_numpy(float)
    a4 = d["a4_Horndeski"].to_numpy(float)
    fp = derivative(r, f)
    fpp = derivative(r, f, 2)
    a4p = derivative(r, a4)

    witnesses = light_ring_witnesses(r, f, h, H, a4, fp, fpp)
    if len(witnesses) != 2:
        raise RuntimeError(f"expected two resolved P5 light rings, got {len(witnesses)}")
    witnesses = sorted(witnesses, key=lambda w: w.radius)
    inner, outer = witnesses

    q = required_electric_q_profile(r, f, h, H, a4, a4p, fp, fpp)
    seed = pd.DataFrame({
        "r_over_rs": r,
        "u": d["u"].to_numpy(float),
        "f": f,
        "h": h,
        "tensor_F_target": H,
        "tensor_H": H,
        "a4": a4,
        "eq85_required_A0prime2_v8": q,
    })
    OUTDIR.mkdir(parents=True, exist_ok=True)
    seed.to_csv(SEED, index=False)

    inner_no_go = bool(inner.a4 > 0.0 and inner.tensor_h > 0.0 and inner.zero_vector_required_tensor_f < 0.0)
    outer_consistent = bool(abs(outer.zero_vector_required_f_over_h - 1.0) < 1e-4)
    report = {
        "status": "PASS_DIAGNOSTIC_ZERO_VECTOR_EXCLUDED_AT_INNER_LIGHT_RING" if inner_no_go else "NO_EXCLUSION",
        "scope": "Local necessary on-shell Eq.85 diagnostic; not a no-go theorem for electric SVT members or P5 geometry",
        "source_csv": str(SOURCE.relative_to(ROOT)),
        "source_sha256": sha256(SOURCE),
        "derivatives": "JET9D8 window=9 degree=8",
        "light_rings": [w.__dict__ for w in witnesses],
        "inner_zero_vector_no_go": inner_no_go,
        "outer_zero_vector_consistency_check": outer_consistent,
        "logic": (
            "At a light ring 2f-rf'=0, Eq.85 loses the a4' control. "
            "For A0prime=0 the inner P5 light ring requires tensor_F<0 while "
            "tensor_H>0 implies a4>0. The electric term -2rfh A0prime^2 v8 "
            "is therefore a necessary available control for a stable hybrid branch."
        ),
        "electric_seed": {
            "path": str(SEED.relative_to(ROOT)),
            "definition": "pointwise Eq.85 solution for q=A0prime^2*v8 with tensor_F target fixed to tensor_H of the carrier witness",
            "min": float(np.min(q)),
            "max": float(np.max(q)),
            "median": float(np.median(q)),
            "positive_fraction": float(np.mean(q > 0.0)),
            "inner_light_ring_q_for_F_equals_H": inner.electric_q_for_f_equals_h,
            "outer_light_ring_q_for_F_equals_H": outer.electric_q_for_f_equals_h,
            "warning": "This seed solves only the necessary background identity for a chosen tensor target; it is not a holonomic action reconstruction.",
        },
    }
    REPORT.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0 if inner_no_go and outer_consistent else 2


if __name__ == "__main__":
    raise SystemExit(main())
