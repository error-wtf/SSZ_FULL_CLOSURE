#!/usr/bin/env python3
"""Deep source-origin and numerical-robustness audit of the u~0.7013 K crossing.

Goals
-----
1. Verify whether the tail ghost is already present in the archived selected
   central 41-stream, independently of the action replay.
2. Cross-check the canonical reducer against the independent kinetic Schur
   implementation.
3. Test derivative-stencil robustness of the zero crossing.
4. Record constraint-pivot behavior and weakest eigenvectors around the crossing.

No repair, fitting, optimization, re-freeze or member selection is performed.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ssz_p5.config import DEFAULT_L  # noqa: E402
from ssz_p5.numerics import module  # noqa: E402
from ssz_p5.production.central_action import central_action_inputs, replay_central_action  # noqa: E402
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central  # noqa: E402
from ssz_p5.production.regional_coefficients import central_selected  # noqa: E402
from ssz_p5.reducer.kinetic_schur import kinetic_schur  # noqa: E402

OUT = ROOT / "data/generated/spectral/CURRENT_MEMBER_SOURCE_ORIGIN_ROBUSTNESS_AUDIT.json"
STENCILS = ((7, 6), (9, 8), (11, 8))


def first_zero(u: np.ndarray, y: np.ndarray) -> float | None:
    order = np.argsort(u)
    uu = np.asarray(u, float)[order]
    yy = np.asarray(y, float)[order]
    for i in range(len(uu) - 1):
        if yy[i] == 0:
            return float(uu[i])
        if np.isfinite(yy[i]) and np.isfinite(yy[i + 1]) and yy[i] * yy[i + 1] < 0:
            t = -yy[i] / (yy[i + 1] - yy[i])
            return float(uu[i] + t * (uu[i + 1] - uu[i]))
    return None


def canonical_K(d, reducer, L: int, window: int, degree: int):
    a = reducer.canonical_audit(d, int(L), window, degree)
    K = np.asarray(a["K"], float)
    Ks = (K + K.swapaxes(1, 2)) / 2.0
    return a, Ks, np.linalg.eigvalsh(Ks)


def stream_audit(d, reducer, L: int, window: int, degree: int) -> dict:
    d = d.sort_values("x").reset_index(drop=True)
    u = d.u.to_numpy(float)
    a, Ks, ev = canonical_K(d, reducer, L, window, degree)
    kmin = ev[:, 0]
    avail = (u >= 0.61) & (u < 0.71)
    prod = (u > 0.62) & (u < 0.70)
    tail = (u >= 0.70) & (u < 0.71)

    # independent Schur replay
    schur = kinetic_schur(d, int(L), window=window, degree=degree)
    Ki = np.asarray(schur["K"], float)
    Kis = (Ki + Ki.swapaxes(1, 2)) / 2.0
    scaled_agreement = float(
        np.max(
            np.abs(Ks[avail] - Kis[avail])
            / np.maximum(1.0, np.abs(Kis[avail]))
        )
    )

    z = first_zero(u[avail], kmin[avail])
    if z is None:
        iz = int(np.argmin(kmin[avail]))
        idx = int(np.flatnonzero(avail)[iz])
    else:
        idx = int(np.argmin(np.abs(u - z)))
    w, V = np.linalg.eigh(Ks[idx])

    Dh1 = np.asarray(schur["Dh1"], float)
    auxdet = np.asarray(schur["auxiliary_determinant"], float)
    return {
        "production_min_K": float(np.min(kmin[prod])),
        "tail_min_K": float(np.min(kmin[tail])),
        "tail_negative_rows": int(np.sum(kmin[tail] <= 0)),
        "first_K_zero_u": z,
        "canonical_vs_independent_schur_max_scaled": scaled_agreement,
        "crossing_grid_u": float(u[idx]),
        "crossing_eigenvalues_K": [float(x) for x in w],
        "weakest_eigenvector": [float(x) for x in V[:, 0]],
        "Dh1_at_crossing": float(Dh1[idx]),
        "auxiliary_determinant_at_crossing": float(auxdet[idx]),
        "min_abs_Dh1_available": float(np.min(np.abs(Dh1[avail]))),
        "min_abs_auxiliary_determinant_available": float(np.min(np.abs(auxdet[avail]))),
    }


def main() -> int:
    reducer = module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")

    raw_selected = central_selected(ROOT).sort_values("x").reset_index(drop=True)
    action_inputs = central_action_inputs(ROOT).sort_values("x").reset_index(drop=True)
    action_replay = replay_central_action(action_inputs).sort_values("x").reset_index(drop=True)
    current = build_onshell_central(ROOT).direct41.sort_values("x").reset_index(drop=True)

    streams = {
        "R0_archived_selected_central_41": raw_selected,
        "R1_same_source_action_replay": action_replay,
        "R2_current_locked_electric": current,
    }

    results = {}
    for sname, d in streams.items():
        per_stencil = {}
        for window, degree in STENCILS:
            key = f"w{window}d{degree}"
            per_stencil[key] = {
                str(L): stream_audit(d, reducer, int(L), window, degree)
                for L in DEFAULT_L
            }
        results[sname] = per_stencil

    # Zero-crossing spread across stencils, per stream/L.
    robustness = {}
    for sname, bystencil in results.items():
        robustness[sname] = {}
        for L in DEFAULT_L:
            vals = [
                bystencil[k][str(L)]["first_K_zero_u"]
                for k in bystencil
                if bystencil[k][str(L)]["first_K_zero_u"] is not None
            ]
            robustness[sname][str(L)] = {
                "n_crossings": len(vals),
                "min_zero_u": float(min(vals)) if vals else None,
                "max_zero_u": float(max(vals)) if vals else None,
                "spread": float(max(vals) - min(vals)) if vals else None,
            }

    # Direct raw-vs-action-replay 41-slot comparison on the common grid.
    common = min(len(raw_selected), len(action_replay))
    slot_diff = {}
    for col in raw_selected.columns:
        if col not in action_replay.columns:
            continue
        try:
            a = raw_selected[col].to_numpy(float)[:common]
            b = action_replay[col].to_numpy(float)[:common]
        except Exception:
            continue
        denom = np.maximum(1.0, np.maximum(np.abs(a), np.abs(b)))
        slot_diff[col] = float(np.max(np.abs(a - b) / denom))
    largest_slot_diffs = sorted(slot_diff.items(), key=lambda kv: kv[1], reverse=True)[:20]

    r0_fail = any(
        results["R0_archived_selected_central_41"]["w9d8"][str(L)]["tail_negative_rows"] > 0
        for L in DEFAULT_L
    )
    r1_fail = any(
        results["R1_same_source_action_replay"]["w9d8"][str(L)]["tail_negative_rows"] > 0
        for L in DEFAULT_L
    )
    r2_fail = any(
        results["R2_current_locked_electric"]["w9d8"][str(L)]["tail_negative_rows"] > 0
        for L in DEFAULT_L
    )

    if r0_fail:
        diagnosis = "TAIL_GHOST_EXISTS_IN_ARCHIVED_SELECTED_CENTRAL_41_STREAM"
        earliest_tracked_stage = "ARCHIVED_SELECTED_CENTRAL_MEMBER"
    elif r1_fail:
        diagnosis = "TAIL_GHOST_INTRODUCED_BY_SOURCE_ACTION_REPLAY"
        earliest_tracked_stage = "CENTRAL_ACTION_REPLAY"
    elif r2_fail:
        diagnosis = "TAIL_GHOST_INTRODUCED_AFTER_SOURCE_REPLAY"
        earliest_tracked_stage = "CURRENT_ELECTRIC_CHAIN"
    else:
        diagnosis = "NO_TAIL_GHOST_REPRODUCED"
        earliest_tracked_stage = None

    payload = {
        "scope": "deep source-origin + numerical robustness; diagnostic only",
        "stencils": [{"window": w, "degree": d} for w, d in STENCILS],
        "streams": results,
        "zero_crossing_robustness": robustness,
        "raw_vs_action_replay_largest_scaled_slot_diffs": largest_slot_diffs,
        "tail_fail_flags_w9d8": {"R0": r0_fail, "R1": r1_fail, "R2": r2_fail},
        "diagnosis": diagnosis,
        "earliest_tracked_stage": earliest_tracked_stage,
        "interpretation_guard": (
            "An archived selected 41-stream failure localizes the problem to or before "
            "that historical member; it does not identify a unique fundamental-action cause."
        ),
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    print(json.dumps(payload, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
