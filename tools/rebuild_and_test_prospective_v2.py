#!/usr/bin/env python3
"""Rebuild the frozen prospective-v2 member and run the real finite-L K gate.

The script rebuilds TRANSITION_PROSPECTIVE_V2_TOTAL41.csv from the frozen
action-first construction, runs the canonical production reducer for every
required L, and cross-checks the trusted dense-domain minima against the frozen
PROSPECTIVE_V2_K_GATE.json witness.

No QNM or spectral semantics are manufactured if the kinetic gate fails.
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
from ssz_p5.numerics import module  # noqa: E402
from ssz_p5.types import Coefficients41, P5Background  # noqa: E402

OUT = ROOT / "data/generated/strong_field_transition_2026-09-21"
CSV = OUT / "TRANSITION_PROSPECTIVE_V2_TOTAL41.csv"
REPORT = OUT / "PROSPECTIVE_V2_K_GATE_REPLAY_2026-09-29.json"
FROZEN_GATE = OUT / "PROSPECTIVE_V2_K_GATE.json"

TRUST_U_MIN = 0.7002175544
TRUST_U_MAX = 0.7079894973743436


def build_member() -> None:
    subprocess.run(
        [sys.executable, str(ROOT / "tools/build_prospective_member_v2.py")],
        cwd=ROOT,
        check=True,
    )


def coefficients_from_frame(d: pd.DataFrame) -> Coefficients41:
    if "phiprime" not in d.columns and "phi_r" in d.columns:
        d = d.rename(columns={"phi_r": "phiprime"})

    required = {
        "x",
        "u",
        "phi",
        "f",
        "h",
        "phiprime",
        "A0prime",
        "X",
        *SLOT_NAMES,
    }
    missing = sorted(required - set(d.columns))
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

    frame = pd.read_csv(CSV).sort_values("x").reset_index(drop=True)
    coeffs = coefficients_from_frame(frame)
    if "phiprime" not in frame.columns:
        frame["phiprime"] = frame["phi_r"]

    frozen = json.loads(FROZEN_GATE.read_text())
    trusted = (
        (coeffs.background.u >= TRUST_U_MIN)
        & (coeffs.background.u <= TRUST_U_MAX + 1e-15)
    )
    if not np.any(trusted):
        raise ValueError("trusted dense domain is empty")

    reducer = module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    per_l: dict[str, dict] = {}
    frozen_replay: dict[str, dict] = {}
    global_min = (float("inf"), None, None)

    for L in DEFAULT_L:
        result = reducer.canonical_audit(frame, int(L))
        k_raw = np.asarray(result["K"], float)
        k_sym = (k_raw + np.swapaxes(k_raw, 1, 2)) / 2
        eig = np.linalg.eigvalsh(k_sym)

        flat = int(np.argmin(eig))
        i, mode = np.unravel_index(flat, eig.shape)
        value = float(eig[i, mode])

        trusted_eig = eig[trusted]
        trusted_u = coeffs.background.u[trusted]
        trusted_flat = int(np.argmin(trusted_eig))
        ti, tmode = np.unravel_index(trusted_flat, trusted_eig.shape)
        trusted_value = float(trusted_eig[ti, tmode])

        expected = float(frozen["per_L_min_K"][str(L)]["min_K"])
        expected_u = float(frozen["per_L_min_K"][str(L)]["u_at_min"])

        frozen_replay[str(L)] = {
            "min_K": trusted_value,
            "u_at_min": float(trusted_u[ti]),
            "mode_index": int(tmode),
            "expected_min_K": expected,
            "expected_u_at_min": expected_u,
            "abs_delta_min_K": abs(trusted_value - expected),
            "abs_delta_u": abs(float(trusted_u[ti]) - expected_u),
        }

        per_l[str(L)] = {
            "min_K": value,
            "r": float(coeffs.background.r[i]),
            "u": float(coeffs.background.u[i]),
            "mode_index": int(mode),
            "negative_rows": int(np.sum(np.min(eig, axis=1) <= 0)),
            "max_K_asymmetry": float(
                np.max(np.abs(k_raw - np.swapaxes(k_raw, 1, 2)))
            ),
        }

        if value < global_min[0]:
            global_min = (value, int(L), int(i))

    passed = all(v["min_K"] > 0 for v in per_l.values())
    replay_pass = all(
        v["abs_delta_min_K"] <= 1e-9 and v["abs_delta_u"] <= 1e-12
        for v in frozen_replay.values()
    )

    payload = {
        "member": "PROSPECTIVE_V2_TRANSVERSE_ZERO",
        "source": str(CSV.relative_to(ROOT)),
        "required_L": list(DEFAULT_L),
        "per_L_min_K": per_l,
        "trusted_dense_domain": [TRUST_U_MIN, TRUST_U_MAX],
        "frozen_gate_replay": frozen_replay,
        "frozen_gate_replay_pass": replay_pass,
        "verdict": "K_PASS" if passed else "K_GHOST_FAIL",
        "chain_termination": (
            "K_PASSED_CONTINUE" if passed else "FIRST_FAILING_GATE_K"
        ),
        "global_min_K": global_min[0],
        "global_min_L": global_min[1],
        "global_min_u": (
            float(coeffs.background.u[global_min[2]])
            if global_min[2] is not None
            else None
        ),
        "downstream_spectral_semantics": "ENABLED" if passed else "BLOCKED_BY_K",
    }

    REPORT.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    print(json.dumps(payload, indent=2, allow_nan=False))

    if not replay_pass:
        return 4
    return 0 if passed else 3


if __name__ == "__main__":
    raise SystemExit(main())
