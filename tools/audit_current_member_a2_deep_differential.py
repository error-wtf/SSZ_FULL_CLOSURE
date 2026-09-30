#!/usr/bin/env python3
"""Deep A2/background differential audit for the current Electric member.

Compares the actual pre-A2 archival source stream used by the current builder
with the freshly reconstructed current pre-Hessian stream.  It then performs
41-slot leave-one-back substitutions on the *current grid* to identify which
coefficient channels carry the newly introduced high-L tail ghost.

This is diagnostic only.  No fitted repair, no new member, no promotion.
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
from ssz_p5.production.full_action_lower import complete_total_action_jets, emit_lower_slots  # noqa: E402

OUT = ROOT / "data/generated/spectral/CURRENT_MEMBER_A2_DEEP_DIFFERENTIAL_AUDIT.json"
RAW_PATH = ROOT / "archive/full_working_snapshot/ssz_p5_integrable_svt_unreduced_even_coefficients_2026-09-12.csv"


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


def current_pre_hessian_stream() -> tuple[pd.DataFrame, pd.DataFrame]:
    build = build_onshell_central(ROOT)
    action = build.action.copy().sort_values("x").reset_index(drop=True)
    for target, delta in (
        ("f2XX", "principal_delta_f2XX"),
        ("f2XF", "principal_delta_f2XF"),
        ("f2FF", "principal_delta_f2FF"),
    ):
        action[target] = action[target].to_numpy(float) - action[delta].to_numpy(float)
    return action, emit_from_action(action)


def interp_column(src: pd.DataFrame, name: str, xnew: np.ndarray) -> np.ndarray:
    s = src.sort_values("x")
    return np.interp(xnew, s.x.to_numpy(float), s[name].to_numpy(float))


def on_current_grid(raw: pd.DataFrame, cur: pd.DataFrame) -> pd.DataFrame:
    out = cur.copy()
    x = cur.x.to_numpy(float)
    for col in cur.columns:
        if col in raw.columns:
            try:
                out[col] = interp_column(raw, col, x)
            except Exception:
                pass
    # Preserve current coordinate columns exactly to eliminate grid-coordinate noise.
    for col in ("x", "u", "phi", "f", "h", "phiprime", "A0prime"):
        if col in cur.columns:
            out[col] = cur[col].to_numpy()
    return out


def metrics(d: pd.DataFrame, reducer, L: int) -> dict:
    d = d.sort_values("x").reset_index(drop=True)
    u = d.u.to_numpy(float)
    a = reducer.canonical_audit(d, int(L))
    K = np.asarray(a["K"], float)
    Ks = (K + K.swapaxes(1, 2)) / 2
    kmin = np.linalg.eigvalsh(Ks)[:, 0]
    prod = (u > 0.62) & (u < 0.70)
    tail = (u >= 0.70) & (u < 0.71)
    j = np.flatnonzero(tail)[np.argmin(kmin[tail])]
    return {
        "production_min_K": float(np.min(kmin[prod])),
        "tail_min_K": float(kmin[j]),
        "u_at_tail_min_K": float(u[j]),
        "tail_negative_rows": int(np.sum(kmin[tail] <= 0)),
    }


def main() -> int:
    raw = pd.read_csv(RAW_PATH)
    action, cur = current_pre_hessian_stream()
    cur = cur.sort_values("x").reset_index(drop=True)

    required = {"x", "u"} | set(SLOT_NAMES)
    missing = sorted(required - set(raw.columns))
    if missing:
        raise RuntimeError(f"raw pre-A2 source lacks reducer slots: {missing}")

    rawg = on_current_grid(raw, cur)
    reducer = module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")

    baseline = {}
    for L in DEFAULT_L:
        baseline[str(L)] = {
            "pre_A2_archival_source_on_current_grid": metrics(rawg, reducer, int(L)),
            "post_A2_pre_hessian_current": metrics(cur, reducer, int(L)),
        }

    # Quantify changed slots over the physical comparison domain.
    mask = (cur.u.to_numpy(float) >= 0.61) & (cur.u.to_numpy(float) < 0.71)
    slot_deltas = {}
    changed_slots = []
    for slot in SLOT_NAMES:
        a = rawg[slot].to_numpy(float)
        b = cur[slot].to_numpy(float)
        diff = b - a
        maxabs = float(np.max(np.abs(diff[mask])))
        rms = float(np.sqrt(np.mean(diff[mask] ** 2)))
        scale = float(np.max(np.maximum(1.0, np.abs(a[mask]))))
        slot_deltas[slot] = {
            "max_abs": maxabs,
            "rms": rms,
            "max_abs_scaled": maxabs / scale,
        }
        if maxabs > 1e-12:
            changed_slots.append(slot)

    # Leave-one-back: current stream with exactly one 41-slot restored to pre-A2.
    # This is a response diagnostic, not an action-consistent member.
    one_slot = {}
    for slot in changed_slots:
        d = cur.copy()
        d[slot] = rawg[slot].to_numpy(float)
        per_L = {}
        for L in DEFAULT_L:
            base = baseline[str(L)]["post_A2_pre_hessian_current"]["tail_min_K"]
            m = metrics(d, reducer, int(L))
            per_L[str(L)] = {
                **m,
                "delta_tail_min_K_vs_current": float(m["tail_min_K"] - base),
                "repairs_tail_K": bool(m["tail_min_K"] > 0),
            }
        one_slot[slot] = per_L

    # Rank by improvement at high L and by cross-L mean improvement.
    ranking = []
    for slot, perL in one_slot.items():
        improvements = [perL[str(L)]["delta_tail_min_K_vs_current"] for L in DEFAULT_L]
        ranking.append({
            "slot": slot,
            "mean_improvement": float(np.mean(improvements)),
            "improvement_L1000": float(perL["1000"]["delta_tail_min_K_vs_current"]),
            "repairs_count": int(sum(perL[str(L)]["repairs_tail_K"] for L in DEFAULT_L)),
        })
    ranking.sort(key=lambda r: (r["repairs_count"], r["improvement_L1000"], r["mean_improvement"]), reverse=True)

    # Restore cumulative top-N channels, to detect distributed rather than single-slot causation.
    cumulative = {}
    top = [r["slot"] for r in ranking[:12]]
    d = cur.copy()
    for n, slot in enumerate(top, start=1):
        d[slot] = rawg[slot].to_numpy(float)
        cumulative[str(n)] = {
            "restored_slots": top[:n],
            "per_L": {str(L): metrics(d, reducer, int(L)) for L in DEFAULT_L},
        }

    # Action/source-variable differences when directly comparable.
    source_vars = {}
    for name in ("f2", "f2X", "f2F", "f3", "v6", "v6_A2_resolved", "f2XX", "f2XF", "f2FF"):
        if name not in action.columns:
            continue
        raw_name = name if name in raw.columns else ("v6" if name == "v6_A2_resolved" and "v6" in raw.columns else None)
        if raw_name is None:
            continue
        x = action.x.to_numpy(float)
        old = interp_column(raw, raw_name, x)
        new = action[name].to_numpy(float)
        dm = new - old
        source_vars[name] = {
            "raw_column": raw_name,
            "max_abs_delta": float(np.max(np.abs(dm))),
            "rms_delta": float(np.sqrt(np.mean(dm * dm))),
            "delta_near_u_0p7013": float(dm[int(np.argmin(np.abs(action.u.to_numpy(float) - 0.7013)))]),
        }

    highL_pre = baseline["1000"]["pre_A2_archival_source_on_current_grid"]["tail_min_K"]
    highL_post = baseline["1000"]["post_A2_pre_hessian_current"]["tail_min_K"]
    diagnosis = (
        "A2_RECONSTRUCTION_INTRODUCES_HIGH_L_TAIL_GHOST"
        if highL_pre > 0 and highL_post < 0
        else "HIGH_L_GHOST_NOT_CLEANLY_LOCALIZED_TO_A2_RECONSTRUCTION"
    )

    payload = {
        "scope": "deep differential audit; no fitting, no repair, no promotion",
        "raw_source": str(RAW_PATH.relative_to(ROOT)),
        "comparison_grid": "current pre-Hessian stream; archival slots linearly interpolated to same x grid",
        "baseline": baseline,
        "changed_slots": changed_slots,
        "slot_deltas": slot_deltas,
        "single_slot_leave_back": one_slot,
        "single_slot_ranking": ranking,
        "cumulative_top12_leave_back": cumulative,
        "action_source_variable_deltas": source_vars,
        "diagnosis": diagnosis,
        "caveat": (
            "single-slot and cumulative leave-back streams are not action-consistent members; "
            "they localize reducer sensitivity only. Baseline pre/post comparison is the physically "
            "relevant construction-stage comparison."
        ),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    print(json.dumps({
        "diagnosis": diagnosis,
        "L1000_pre_A2_tail_min_K": highL_pre,
        "L1000_post_A2_tail_min_K": highL_post,
        "top12": ranking[:12],
        "source_var_deltas": source_vars,
    }, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
