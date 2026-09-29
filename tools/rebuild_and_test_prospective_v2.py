#!/usr/bin/env python3
"""Rebuild the frozen prospective-v2 member and run the real finite-L K gate.

This is not a mock/proxy test.  It rebuilds TRANSITION_PROSPECTIVE_V2_TOTAL41.csv
from the frozen action-first construction, converts that emitted 41-slot stream
into the package Coefficients41 object, and invokes the production reducer for
every required L.

The chain stops at the first physical failing layer.  If K is positive for all
required L, the report explicitly says that downstream angular/spectral work is
allowed.  Otherwise it records the kinetic failure and does not manufacture QNM
semantics.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ssz_p5.config import DEFAULT_L, SLOT_NAMES  # noqa: E402
from ssz_p5.reducer.canonical import reduce_profile  # noqa: E402
from ssz_p5.types import Coefficients41, P5Background  # noqa: E402

CSV = ROOT / "data/generated/strong_field_transition_2026-09-21/TRANSITION_PROSPECTIVE_V2_TOTAL41.csv"
REPORT = ROOT / "data/generated/strong_field_transition_2026-09-21/PROSPECTIVE_V2_K_GATE_REPLAY_2026-09-29.json"


def build_member() -> None:
    subprocess.run([sys.executable, str(ROOT / "tools/build_prospective_member_v2.py")],
                   cwd=ROOT, check=True)


def coefficients_from_frame(d: pd.DataFrame) -> Coefficients41:
    if "phiprime" not in d.columns and "phi_r" in d.columns:
        d = d.rename(columns={"phi_r": "phiprime"})
    need = {"x", "u", "phi", "f", "h", "phiprime", "A0prime", "X", *SLOT_NAMES}
    missing = sorted(need - set(d.columns))
    if missing:
        raise ValueError(f"member stream missing columns: {missing}")
    d = d.sort_values("x").reset_index(drop=True)
    r = d.x.to_numpy(float)
    if np.any(np.diff(r) <= 0):
        raise ValueError("radial grid is not strictly increasing")
    bg = P5Background(
        r=r,
        u=d.u.to_numpy(float),
        phi=d.phi.to_numpy(float),
        f=d.f.to_numpy(float),
        h=d.h.to_numpy(float),
        phi_r=d.phiprime.to_numpy(float),
        X=d.X.to_numpy(float),
        A0prime=d.A0prime.to_numpy(float),
        region=d.region.to_numpy(str) if "region" in d.columns else None,
    )
    slots = {name: d[name].to_numpy(float) for name in SLOT_NAMES}
    return Coefficients41(bg, slots, "PROSPECTIVE_V2_TRANSVERSE_ZERO_REBUILT")


def main() -> int:
    build_member()
    d = pd.read_csv(CSV)
    coeffs = coefficients_from_frame(d)
    reducer = module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    frame = d.copy().sort_values("x").reset_index(drop=True)
    if "phiprime" not in frame.columns:
        frame["phiprime"] = frame["phi_r"]
    per_l = {}
    global_min = (float("inf"), None, None)
    for L in DEFAULT_L:
        result = reducer.canonical_audit(frame, int(L))
        Kraw = np.asarray(result["K"], float)
        K = (Kraw + np.swapaxes(Kraw, 1, 2)) / 2
        eig = np.linalg.eigvalsh(K)
        flat = int(np.argmin(eig))
        i, mode = np.unravel_index(flat, eig.shape)
        value = float(eig[i, mode])
        per_l[str(L)] = {
            "min_K": value,
            "r": float(coeffs.background.r[i]),
            "u": float(coeffs.background.u[i]),
            "mode_index": int(mode),
            "negative_rows": int(np.sum(np.min(eig, axis=1) <= 0)),
            "max_K_asymmetry": float(np.max(np.abs(Kraw - np.swapaxes(Kraw, 1, 2)))),
        }
        if value < global_min[0]:
            global_min = (value, int(L), int(i))

    passed = all(v["min_K"] > 0 for v in per_l.values())
    payload = {
        "member": "PROSPECTIVE_V2_TRANSVERSE_ZERO",
        "source": str(CSV.relative_to(ROOT)),
        "required_L": list(DEFAULT_L),
        "per_L_min_K": per_l,
        "verdict": "K_PASS" if passed else "K_GHOST_FAIL",
        "chain_termination": "K_PASSED_CONTINUE" if passed else "FIRST_FAILING_GATE_K",
        "global_min_K": global_min[0],
        "global_min_L": global_min[1],
        "global_min_u": (
            float(coeffs.background.u[global_min[2]]) if global_min[2] is not None else None
        ),
        "downstream_spectral_semantics": "ENABLED" if passed else "BLOCKED_BY_K",
    }
    REPORT.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    print(json.dumps(payload, indent=2, allow_nan=False))
    return 0 if passed else 3


if __name__ == "__main__":
    raise SystemExit(main())
