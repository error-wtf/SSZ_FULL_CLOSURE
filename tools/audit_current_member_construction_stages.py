#!/usr/bin/env python3
"""Stage-by-stage origin audit for the current Electric-member tail ghost.

Purpose
-------
Locate the earliest *reproducible construction stage* at which the finite-L
kinetic eigenvalue becomes negative near u~0.7013.

This is a falsification/diagnostic tool only:
* no optimization;
* no fitting;
* no new member selection;
* no re-freeze;
* every counterfactual is explicitly labeled.

Stages tested
-------------
S0  archived/current reference streams, if they expose a reducer-compatible
    41-slot table;
S1  current Electric action after A2/background solve with the stored principal
    Hessian restore exactly undone, then holonomically re-completed and freshly
    re-emitted from the action;
S2  locked current pre-G2XX 41-stream (post Hessian restore + lower slots);
S3  locked current full stream (post registered G2XX lift).

S1->S2 localizes whether the principal Hessian restoration is responsible.
S2->S3 was already tested separately; this script reproduces that conclusion.
S0 comparisons are provenance/context only and are NOT exact one-variable
counterfactuals unless explicitly reported as such.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ssz_p5.config import DEFAULT_L, SLOT_NAMES  # noqa: E402
from ssz_p5.numerics import module  # noqa: E402
from ssz_p5.production.electric_hybrid_onshell_central import (  # noqa: E402
    build_onshell_central,
)
from ssz_p5.production.full_action_lower import (  # noqa: E402
    complete_total_action_jets,
    emit_lower_slots,
)

OUT = ROOT / "data/generated/spectral/CURRENT_MEMBER_CONSTRUCTION_STAGE_AUDIT.json"

REFERENCE_FILES = {
    "historical_central_search_baseline": ROOT
    / "data/generated/absolute_attempt_2026-09-19/ELECTRIC_HYBRID_CENTRAL_SEARCH_BASELINE.csv",
    "historical_principal_c6_corrected": ROOT
    / "data/generated/absolute_attempt_2026-09-20/ELECTRIC_HYBRID_PRINCIPAL_LATE_RAMP_FEASIBILITY_C6_CORRECTED.csv",
    "archival_integrable_unreduced_even": ROOT
    / "archive/full_working_snapshot/ssz_p5_integrable_svt_unreduced_even_coefficients_2026-09-12.csv",
}


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


def audit_stream(d: pd.DataFrame, reducer, L: int) -> dict:
    d = d.sort_values("x").reset_index(drop=True)
    u = d.u.to_numpy(float)
    a = reducer.canonical_audit(d, int(L))
    K = np.asarray(a["K"], float)
    G = np.asarray(a["G"], float)
    Ks = (K + K.swapaxes(1, 2)) / 2.0
    Gs = (G + G.swapaxes(1, 2)) / 2.0
    kvals = np.linalg.eigvalsh(Ks)
    kmin = kvals[:, 0]

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
    avail = (u >= 0.61) & (u < 0.71)
    return {
        "rows": int(len(d)),
        "u_min": float(u.min()),
        "u_max": float(u.max()),
        "production_min_K": float(np.min(kmin[prod])) if np.any(prod) else None,
        "production_negative_K_rows": int(np.sum(kmin[prod] <= 0)) if np.any(prod) else None,
        "tail_min_K": float(np.min(kmin[tail])) if np.any(tail) else None,
        "tail_negative_K_rows": int(np.sum(kmin[tail] <= 0)) if np.any(tail) else None,
        "available_min_K": float(np.min(kmin[avail])) if np.any(avail) else None,
        "first_K_zero_u": first_zero(u[avail], kmin[avail]) if np.any(avail) else None,
        "min_cr2_where_K_positive": (
            float(np.nanmin(cr2[avail])) if np.any(np.isfinite(cr2[avail])) else None
        ),
        "first_cr2_zero_u": (
            first_zero(u[avail][np.isfinite(cr2[avail])], cr2[avail][np.isfinite(cr2[avail])])
            if np.any(np.isfinite(cr2[avail])) else None
        ),
    }


def emit_from_action(action: pd.DataFrame) -> pd.DataFrame:
    """Freshly complete and emit a same-action 41-stream."""
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


def compatible_reference(path: Path) -> tuple[bool, str | None, pd.DataFrame | None]:
    if not path.exists():
        return False, "missing_file", None
    try:
        d = pd.read_csv(path)
    except Exception as exc:
        return False, f"read_error:{type(exc).__name__}", None
    required = {"u", "x"} | set(SLOT_NAMES)
    missing = sorted(required - set(d.columns))
    if missing:
        return False, "missing_columns:" + ",".join(missing[:20]), None
    return True, None, d


def main() -> int:
    build = build_onshell_central(ROOT)
    current_action = build.action.copy().sort_values("x").reset_index(drop=True)

    # S1: exactly undo only the stored principal Hessian restore.
    required_deltas = (
        ("f2XX", "principal_delta_f2XX"),
        ("f2XF", "principal_delta_f2XF"),
        ("f2FF", "principal_delta_f2FF"),
    )
    for target, delta in required_deltas:
        if target not in current_action or delta not in current_action:
            raise RuntimeError(f"missing stored principal restore provenance: {target}/{delta}")

    pre_hessian_action = current_action.copy()
    undo_summary = {}
    for target, delta in required_deltas:
        dv = pre_hessian_action[delta].to_numpy(float)
        pre_hessian_action[target] = pre_hessian_action[target].to_numpy(float) - dv
        undo_summary[target] = {
            "max_abs_undone": float(np.max(np.abs(dv))),
            "rms_undone": float(np.sqrt(np.mean(dv * dv))),
        }

    stages = {
        "S1_current_after_A2_before_principal_hessian_restore": emit_from_action(pre_hessian_action),
        "S2_current_post_hessian_pre_G2XX": build.pre_lift41.copy(),
        "S3_current_locked_post_G2XX": build.direct41.copy(),
    }

    reducer = module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    results = {}
    for name, stream in stages.items():
        results[name] = {
            str(L): audit_stream(stream, reducer, int(L)) for L in DEFAULT_L
        }

    refs = {}
    for name, path in REFERENCE_FILES.items():
        ok, reason, d = compatible_reference(path)
        item = {
            "path": str(path.relative_to(ROOT)),
            "reducer_compatible_41_stream": bool(ok),
            "reason_if_not": reason,
        }
        if ok and d is not None:
            item["per_L"] = {
                str(L): audit_stream(d, reducer, int(L)) for L in DEFAULT_L
            }
        refs[name] = item

    def tail_fail(stage: str) -> bool:
        return any(
            row["tail_negative_K_rows"] is not None and row["tail_negative_K_rows"] > 0
            for row in results[stage].values()
        )

    s1_fail = tail_fail("S1_current_after_A2_before_principal_hessian_restore")
    s2_fail = tail_fail("S2_current_post_hessian_pre_G2XX")
    s3_fail = tail_fail("S3_current_locked_post_G2XX")

    if s1_fail:
        diagnosis = "GHOST_PRESENT_BEFORE_PRINCIPAL_HESSIAN_RESTORE"
        next_target = "A2_BACKGROUND_SOLVE_OR_EARLIER_ACTION_INPUT"
    elif s2_fail:
        diagnosis = "GHOST_FIRST_APPEARS_WITH_PRINCIPAL_HESSIAN_RESTORE"
        next_target = "PRINCIPAL_HESSIAN_RESTORE"
    elif s3_fail:
        diagnosis = "GHOST_FIRST_APPEARS_WITH_G2XX_LIFT"
        next_target = "REGISTERED_G2XX_LIFT"
    else:
        diagnosis = "NO_CURRENT_STAGE_REPRODUCES_TAIL_GHOST"
        next_target = "NUMERICAL_OR_PROVENANCE_MISMATCH"

    payload = {
        "scope": "stage-localization diagnostic; no fitting, no repair, no member selection",
        "undo_principal_hessian_restore": undo_summary,
        "stage_semantics": {
            "S1": "current Electric action after A2/background solve with stored principal Hessian restore undone; holonomically re-completed and freshly emitted",
            "S2": "locked pre-G2XX current stream",
            "S3": "locked full current stream",
            "S0_references": "historical/provenance context only; not exact one-variable counterfactuals",
        },
        "current_stages": results,
        "reference_streams": refs,
        "tail_fail_flags": {"S1": s1_fail, "S2": s2_fail, "S3": s3_fail},
        "diagnosis": diagnosis,
        "next_target": next_target,
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    print(json.dumps(payload, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
