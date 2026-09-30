#!/usr/bin/env python3
"""Diagnose the u>0.70 instability of the locked current electric member.

No repair or fitting is performed.  The script decomposes the already frozen
member into its pre-G2XX stream and the registered background-null G2XX lift,
then measures the finite-L kinetic/radial response and localizes the first
zero crossing of lambda_min(K).  A scale scan of the *existing fixed lift
direction* is diagnostic only and must not be used to select a new member.
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
from ssz_p5.production.electric_hybrid_onshell_central import (  # noqa: E402
    build_onshell_central,
)

OUT = ROOT / "data/generated/spectral/CURRENT_MEMBER_TAIL_FAILURE_DIAGNOSTIC.json"


def first_zero(u: np.ndarray, y: np.ndarray) -> float | None:
    """Linear interpolation of the first sign change in ascending u."""
    order = np.argsort(u)
    uu = u[order]
    yy = y[order]
    for i in range(len(uu) - 1):
        if yy[i] == 0:
            return float(uu[i])
        if yy[i] * yy[i + 1] < 0:
            t = -yy[i] / (yy[i + 1] - yy[i])
            return float(uu[i] + t * (uu[i + 1] - uu[i]))
    return None


def stream_at_scale(pre: pd.DataFrame, full: pd.DataFrame, scale: float) -> pd.DataFrame:
    """Interpolate only the registered G2XX primitive response in 41-space."""
    out = pre.copy()
    for name in ("c2", "c6", "e2"):
        out[name] = pre[name].to_numpy(float) + scale * (
            full[name].to_numpy(float) - pre[name].to_numpy(float)
        )
    return out


def audit_stream(d: pd.DataFrame, reducer, L: int) -> dict:
    u = d.u.to_numpy(float)
    a = reducer.canonical_audit(d, int(L))
    K = np.asarray(a["K"], float)
    G = np.asarray(a["G"], float)
    Ks = (K + K.swapaxes(1, 2)) / 2.0
    Gs = (G + G.swapaxes(1, 2)) / 2.0
    evals = np.linalg.eigvalsh(Ks)
    kmin = evals[:, 0]

    cr2 = np.full(len(d), np.nan)
    for i in range(len(d)):
        w, U = np.linalg.eigh(Ks[i])
        if np.min(w) <= 0:
            continue
        invsqrt = U @ np.diag(1.0 / np.sqrt(w)) @ U.T
        C = invsqrt @ Gs[i] @ invsqrt
        cr2[i] = np.linalg.eigvalsh((C + C.T) / 2.0)[0]

    prod = (u > 0.62) & (u < 0.70)
    tail = u >= 0.70
    return {
        "global_min_K": float(np.min(kmin)),
        "u_at_global_min_K": float(u[int(np.argmin(kmin))]),
        "production_min_K": float(np.min(kmin[prod])),
        "tail_min_K": float(np.min(kmin[tail])),
        "tail_negative_K_rows": int(np.sum(kmin[tail] <= 0)),
        "first_K_zero_u": first_zero(u, kmin),
        "min_cr2_where_defined": (
            float(np.nanmin(cr2)) if np.any(np.isfinite(cr2)) else None
        ),
        "first_cr2_zero_u": first_zero(u[np.isfinite(cr2)], cr2[np.isfinite(cr2)]),
    }


def main() -> int:
    build = build_onshell_central(ROOT)
    pre = build.pre_lift41.sort_values("x").reset_index(drop=True)
    full = build.direct41.sort_values("x").reset_index(drop=True)
    np.testing.assert_allclose(pre.u, full.u, rtol=0, atol=0)

    reducer = module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    scales = (0.0, 0.25, 0.5, 0.75, 1.0)
    per_l = {}
    for L in DEFAULT_L:
        rows = {}
        for scale in scales:
            rows[f"{scale:.2f}"] = audit_stream(
                stream_at_scale(pre, full, scale), reducer, int(L)
            )
        base = rows["0.00"]
        locked = rows["1.00"]
        per_l[str(L)] = {
            "scale_scan": rows,
            "pre_lift_tail_K_pass": base["tail_negative_K_rows"] == 0,
            "locked_tail_K_pass": locked["tail_negative_K_rows"] == 0,
            "delta_tail_min_K_locked_minus_pre": (
                locked["tail_min_K"] - base["tail_min_K"]
            ),
        }

    # exact 41-stream response provenance
    deltas = {}
    for name in full.columns:
        if name in pre.columns:
            try:
                a = pre[name].to_numpy(float)
                b = full[name].to_numpy(float)
            except Exception:
                continue
            diff = np.max(np.abs(b - a))
            if diff > 0:
                deltas[name] = float(diff)

    payload = {
        "scope": (
            "Diagnostic decomposition of the locked current member only; "
            "no repair, no parameter selection, no physical re-freeze."
        ),
        "registered_lift_direction": "background-null G2XX primitive",
        "changed_41_slots_pre_to_locked": deltas,
        "scale_scan_semantics": (
            "counterfactual sensitivity along the already registered lift direction; "
            "NOT a search and NOT admissible as a new member-selection rule"
        ),
        "per_L": per_l,
        "all_pre_lift_tail_K_pass": all(
            x["pre_lift_tail_K_pass"] for x in per_l.values()
        ),
        "all_locked_tail_K_pass": all(
            x["locked_tail_K_pass"] for x in per_l.values()
        ),
    }

    if payload["all_pre_lift_tail_K_pass"] and not payload["all_locked_tail_K_pass"]:
        payload["diagnosis"] = "TAIL_GHOST_INDUCED_BY_REGISTERED_G2XX_LIFT_RESPONSE"
    elif not payload["all_pre_lift_tail_K_pass"]:
        payload["diagnosis"] = "TAIL_GHOST_ALREADY_PRESENT_BEFORE_G2XX_LIFT"
    else:
        payload["diagnosis"] = "NO_TAIL_K_GHOST_REPRODUCED"

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    print(json.dumps(payload, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
