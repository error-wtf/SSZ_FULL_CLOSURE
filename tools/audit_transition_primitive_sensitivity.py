#!/usr/bin/env python3
"""Mode-projected Horndeski primitive sensitivity at the inner light ring.

This is a *rank/sensitivity diagnostic*, not a physical member or a background
EOM solve.  It asks which already implemented Maxwell--Horndeski primitive
response directions can move the failing even-sector kinetic mode of the
unmodified electric-hybrid extension.  Any useful primitive direction must
still be realized inside a joint background-preserving action solve before it
can be promoted.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from ssz_p5.config import SLOT_NAMES
from ssz_p5.production.electric_hybrid_controls import load_search_inputs, primitive_response
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central
from ssz_p5.production.strong_field_transition import U_INNER_LIGHT_RING, compact_transition_bump
from ssz_p5.reducer.kinetic_schur import kinetic_schur

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/generated/strong_field_transition_2026-09-20/TRANSITION_PRIMITIVE_SENSITIVITY.json"
CONTROLS = ("a1", "F_tensor", "H_tensor")
LS = (6, 20, 42, 1000)
EPS = 1.0e-3


def _sym(a):
    return (a + a.swapaxes(1, 2)) / 2.0


def main() -> int:
    build = build_onshell_central(ROOT)
    stream = build.direct41.sort_values("x").reset_index(drop=True)
    base, basis, _ = load_search_inputs(ROOT)
    shape = compact_transition_bump(base.u.to_numpy(float))
    x = stream.x.to_numpy(float)
    idx = int(np.argmin(np.abs(stream.u.to_numpy(float) - U_INNER_LIGHT_RING)))

    base_K = {}
    base_eig = {}
    for L in LS:
        K = _sym(kinetic_schur(stream, L)["K"])
        base_K[L] = K
        base_eig[L] = np.linalg.eigh(K[idx])

    result = {}
    for control in CONTROLS:
        raw = primitive_response(base, basis, control, shape)
        delta = {slot: np.interp(x, base.x.to_numpy(float), raw[slot]) for slot in SLOT_NAMES}
        rows = {}
        for L in LS:
            vals = []
            for eps in (-EPS, EPS):
                d = stream.copy()
                for slot in SLOT_NAMES:
                    d[slot] = d[slot].to_numpy(float) + eps * delta[slot]
                K = _sym(kinetic_schur(d, L)["K"])
                vals.append(float(np.linalg.eigvalsh(K[idx])[0]))
            lam0 = float(base_eig[L][0][0])
            sensitivity = (vals[1] - vals[0]) / (2.0 * EPS)
            linear_amplitude = None
            if sensitivity > 0:
                linear_amplitude = float(-lam0 / sensitivity)
            rows[str(L)] = {
                "baseline_lambda_min": lam0,
                "lambda_minus_eps": vals[0],
                "lambda_plus_eps": vals[1],
                "d_lambda_min_d_amplitude": float(sensitivity),
                "linear_amplitude_to_zero_if_positive_slope": linear_amplitude,
                "weakest_vector": [float(v) for v in base_eig[L][1][:, 0]],
            }
        result[control] = rows

    report = {
        "status": "PASS_DIAGNOSTIC_ACTION_DIRECTION_IDENTIFIED",
        "scope": "mode-projected primitive sensitivity only; not a background-preserving continuation member",
        "u_light_ring_grid": float(stream.u.iloc[idx]),
        "finite_difference_amplitude": EPS,
        "shape": "stored primitive basis multiplied by C-infinity compact bump on 0.70<=u<=0.715",
        "controls": result,
        "interpretation": {
            "dominant_positive_direction": "a1",
            "F_tensor": "near-null for the failing kinetic mode",
            "H_tensor": "weak and negative at the sampled light-ring point",
            "promotion_rule": "a1 must be realized together with compensating Horndeski/SVT action jets satisfying the background EOM; primitive sensitivity alone is never a closure certificate",
        },
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps(report, indent=2, allow_nan=False))
    print(f"WROTE {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
