#!/usr/bin/env python3
"""Localize the current-member ghost relative to the Electric A2/background solve.

This audit compares:
P0  the registered central exact-SVT action inputs replayed directly from their
    source action (pre current Electric A2 re-solve);
P1  the current Electric action after the A2/background solve with the later
    principal Hessian restoration exactly undone and with no G2XX lift.

If P0 is healthy in the tail and P1 is not, the first reproducible appearance
of the ghost is the Electric A2/background re-solve block.  If P0 already
fails, the origin lies in the earlier central source action.

No parameter search, no fitting and no member mutation are performed.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ssz_p5.config import DEFAULT_L  # noqa: E402
from ssz_p5.numerics import module  # noqa: E402
from ssz_p5.production.central_action import (  # noqa: E402
    background_residuals,
    central_action_inputs,
    replay_central_action,
)
from ssz_p5.production.electric_hybrid_onshell_central import (  # noqa: E402
    build_onshell_central,
)
from ssz_p5.production.full_action_lower import (  # noqa: E402
    complete_total_action_jets,
    emit_lower_slots,
)

OUT = ROOT / "data/generated/spectral/CURRENT_MEMBER_A2_ORIGIN_AUDIT.json"


def first_zero(u: np.ndarray, y: np.ndarray) -> float | None:
    order = np.argsort(u)
    uu, yy = np.asarray(u)[order], np.asarray(y)[order]
    for i in range(len(uu) - 1):
        if yy[i] == 0:
            return float(uu[i])
        if np.isfinite(yy[i]) and np.isfinite(yy[i + 1]) and yy[i] * yy[i + 1] < 0:
            t = -yy[i] / (yy[i + 1] - yy[i])
            return float(uu[i] + t * (uu[i + 1] - uu[i]))
    return None


def audit_stream(d: pd.DataFrame, reducer, L: int) -> dict:
    d = d.sort_values("x").reset_index(drop=True)
    u = d.u.to_numpy(float)
    a = reducer.canonical_audit(d, int(L))
    K = np.asarray(a["K"], float)
    Ks = (K + K.swapaxes(1, 2)) / 2.0
    k = np.linalg.eigvalsh(Ks)[:, 0]
    prod = (u > 0.62) & (u < 0.70)
    tail = (u >= 0.70) & (u < 0.71)
    avail = (u >= 0.61) & (u < 0.71)
    return {
        "production_min_K": float(np.min(k[prod])),
        "production_negative_K_rows": int(np.sum(k[prod] <= 0)),
        "tail_min_K": float(np.min(k[tail])),
        "tail_negative_K_rows": int(np.sum(k[tail] <= 0)),
        "first_K_zero_u": first_zero(u[avail], k[avail]),
    }


def emit_from_action(action: pd.DataFrame) -> pd.DataFrame:
    zk = module("ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py")
    d = complete_total_action_jets(action.copy().sort_values("x").reset_index(drop=True))
    lower, d = emit_lower_slots(d)
    out = zk.emit(
        d,
        selected_v5=lower.v5.to_numpy(float),
        selected_c3=lower.c3.to_numpy(float),
        selected_e3=lower.e3.to_numpy(float),
        v6_phi_selector="action",
    )
    r = out.x.to_numpy(float)
    out["a5"] = (
        zk.dr(r, out.a2.to_numpy(float), 1, 9, 8)
        - zk.dr(r, out.a1.to_numpy(float), 2, 9, 8)
        - zk.dr(r, out.A0prime.to_numpy(float) * out.v4.to_numpy(float) / 2.0, 1, 9, 8)
        + out.A0prime.to_numpy(float) * out.v5.to_numpy(float) / 2.0
    )
    return out


def nearest_row(d: pd.DataFrame, u0: float) -> int:
    return int(np.argmin(np.abs(d.u.to_numpy(float) - u0)))


def scaled_change(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.max(np.abs(a - b) / np.maximum(1.0, np.maximum(np.abs(a), np.abs(b)))))


def main() -> int:
    reducer = module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")

    # P0: source central exact-SVT action, no current Electric A2 re-solve.
    p0_action = central_action_inputs(ROOT).sort_values("x").reset_index(drop=True)
    p0_stream = replay_central_action(p0_action)

    # P1: current Electric action after A2/background solve, but undo the later
    # principal Hessian restore so this comparison is upstream-clean.
    build = build_onshell_central(ROOT)
    p1_action = build.action.copy().sort_values("x").reset_index(drop=True)
    for target, delta in (
        ("f2XX", "principal_delta_f2XX"),
        ("f2XF", "principal_delta_f2XF"),
        ("f2FF", "principal_delta_f2FF"),
    ):
        p1_action[target] = p1_action[target].to_numpy(float) - p1_action[delta].to_numpy(float)
    p1_stream = emit_from_action(p1_action)

    per_l = {
        "P0_source_central_action": {
            str(L): audit_stream(p0_stream, reducer, int(L)) for L in DEFAULT_L
        },
        "P1_post_A2_background_solve_pre_hessian": {
            str(L): audit_stream(p1_stream, reducer, int(L)) for L in DEFAULT_L
        },
    }

    p0_fail = any(v["tail_negative_K_rows"] > 0 for v in per_l["P0_source_central_action"].values())
    p1_fail = any(v["tail_negative_K_rows"] > 0 for v in per_l["P1_post_A2_background_solve_pre_hessian"].values())

    if (not p0_fail) and p1_fail:
        diagnosis = "GHOST_FIRST_APPEARS_IN_ELECTRIC_A2_BACKGROUND_RESOLVE_BLOCK"
        next_target = "DECOMPOSE_A2_V6_F3_F2F_VS_METRIC_F2_F2X_RESOLVE"
    elif p0_fail:
        diagnosis = "GHOST_ALREADY_PRESENT_IN_SOURCE_CENTRAL_ACTION"
        next_target = "EARLIER_CENTRAL_ACTION_INPUTS"
    elif not p1_fail:
        diagnosis = "A2_BLOCK_NOT_REPRODUCING_TAIL_GHOST"
        next_target = "RECHECK_STAGE_PROVENANCE"
    else:
        diagnosis = "UNCLASSIFIED"
        next_target = "MANUAL_REVIEW"

    # Background equations on the two action-level stages.
    bg = {}
    for name, action in {
        "P0_source_central_action": p0_action,
        "P1_post_A2_background_solve_pre_hessian": p1_action,
    }.items():
        vals = background_residuals(action)
        u = action.u.to_numpy(float)
        region = (u >= 0.61) & (u < 0.71)
        bg[name] = {
            key: float(np.max(np.abs(np.asarray(val)[region])))
            for key, val in vals.items()
        }

    # Action-field differences on the common overlap, including a witness around
    # the first current ghost crossing.
    p0 = p0_action.copy()
    p1 = p1_action.copy()
    common_u = p1.u.to_numpy(float)
    action_changes = {}
    for col in ("f2", "f2X", "f2F", "f3", "f3X", "f4", "f4X", "f4XX", "A0prime"):
        if col not in p0.columns or col not in p1.columns:
            continue
        a = np.interp(common_u, p0.u.to_numpy(float), p0[col].to_numpy(float))
        b = p1[col].to_numpy(float)
        action_changes[col] = {
            "max_scaled_change": scaled_change(a, b),
            "max_abs_change": float(np.max(np.abs(a - b))),
        }

    ghost_us = [
        v["first_K_zero_u"]
        for v in per_l["P1_post_A2_background_solve_pre_hessian"].values()
        if v["first_K_zero_u"] is not None
    ]
    witness_u = float(np.median(ghost_us)) if ghost_us else 0.7013
    i0 = nearest_row(p0, witness_u)
    i1 = nearest_row(p1, witness_u)
    witness = {
        "u_target": witness_u,
        "P0_u": float(p0.u.iloc[i0]),
        "P1_u": float(p1.u.iloc[i1]),
        "fields": {},
    }
    for col in action_changes:
        witness["fields"][col] = {
            "P0": float(p0[col].iloc[i0]),
            "P1": float(p1[col].iloc[i1]),
            "delta": float(p1[col].iloc[i1] - p0[col].iloc[i0]),
        }

    payload = {
        "scope": "A2/background-solve origin audit; diagnostic only; no fitting/repair",
        "per_L": per_l,
        "tail_fail_flags": {"P0": p0_fail, "P1": p1_fail},
        "background_residual_max_abs": bg,
        "action_field_changes": action_changes,
        "ghost_crossing_witness": witness,
        "diagnosis": diagnosis,
        "next_target": next_target,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    print(json.dumps(payload, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
