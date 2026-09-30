#!/usr/bin/env python3
"""Action-variable causality audit for the new u~0.7013 high-L kinetic ghost.

The 41-slot differential audit identifies e1 as the only single-slot leave-back
that repairs all tested L.  Appendix-A gives

    e1 = (a2 - 2 r h a6 - A0' v4/2)/(f h phi')

with
    a2 = ... + A0' v4/2 + r A0' v6/(2 phi')
    a6 = ... - A0' v6/(4 h phi') + ...
    v6 = 4 alpha2
       = 2 h^(3/2) A0'/(r sqrt(f))
         * (r phi' f3 - 4 f4 + h phi'^2 (f4X + 2 tf4)).

Therefore v4 cancels identically from e1 and, with the background/f4 sector
held fixed, the Electric reconstruction can change e1 through f3.

This audit tests that implication directly at action level.  It does not fit or
promote a repaired member.  Counterfactual substitutions are explicitly
labelled diagnostic-only.
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
from ssz_p5.jets.jet9d8 import profile_derivative  # noqa: E402
from ssz_p5.numerics import module  # noqa: E402
from ssz_p5.production.central_action import central_action_inputs  # noqa: E402
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central  # noqa: E402
from ssz_p5.production.full_action_lower import complete_total_action_jets, emit_lower_slots  # noqa: E402

OUT = ROOT / "data/generated/spectral/CURRENT_MEMBER_ACTION_VARIABLE_CAUSALITY_AUDIT.json"


def interp(src: pd.DataFrame, col: str, xnew: np.ndarray) -> np.ndarray:
    s = src.sort_values("x")
    return np.interp(xnew, s.x.to_numpy(float), s[col].to_numpy(float))


def prepare_current_pre_hessian() -> pd.DataFrame:
    b = build_onshell_central(ROOT)
    d = b.action.copy().sort_values("x").reset_index(drop=True)
    for target, delta in (
        ("f2XX", "principal_delta_f2XX"),
        ("f2XF", "principal_delta_f2XF"),
        ("f2FF", "principal_delta_f2FF"),
    ):
        d[target] = d[target].to_numpy(float) - d[delta].to_numpy(float)
    return d


def refresh_f2phi(d: pd.DataFrame) -> pd.DataFrame:
    out = d.copy().sort_values("x").reset_index(drop=True)
    r = out.x.to_numpy(float)
    f, h, A = (out[c].to_numpy(float) for c in ("f", "h", "A0prime"))
    X = out.X.to_numpy(float)
    ph = out.phiprime.to_numpy(float)
    F = h * A**2 / (2.0 * f)
    Xp = profile_derivative(r, X, 1, 9, 8)
    Fp = profile_derivative(r, F, 1, 9, 8)
    f2p = profile_derivative(r, out.f2.to_numpy(float), 1, 9, 8)
    f2Y = out.f2Y.to_numpy(float) if "f2Y" in out else np.zeros(len(out))
    # Y is present but f2Y=0 on the current production branch.
    Y = 4.0 * X * F
    Yp = profile_derivative(r, Y, 1, 9, 8)
    out["f2phi"] = (
        f2p
        - out.f2X.to_numpy(float) * Xp
        - out.f2F.to_numpy(float) * Fp
        - f2Y * Yp
    ) / ph
    return out


def emit_action(d: pd.DataFrame) -> pd.DataFrame:
    zk = module("ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py")
    d = refresh_f2phi(d)
    d = complete_total_action_jets(d)
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


def metrics(stream: pd.DataFrame, reducer, L: int) -> dict:
    d = stream.sort_values("x").reset_index(drop=True)
    u = d.u.to_numpy(float)
    a = reducer.canonical_audit(d, int(L))
    K = np.asarray(a["K"], float)
    Ks = (K + K.swapaxes(1, 2)) / 2.0
    kmin = np.linalg.eigvalsh(Ks)[:, 0]
    prod = (u > 0.62) & (u < 0.70)
    tail = (u >= 0.70) & (u < 0.71)
    return {
        "production_min_K": float(np.min(kmin[prod])),
        "tail_min_K": float(np.min(kmin[tail])),
        "tail_negative_rows": int(np.sum(kmin[tail] <= 0)),
    }


def e1_identity_residual(stream: pd.DataFrame) -> dict:
    d = stream.sort_values("x").reset_index(drop=True)
    r = d.x.to_numpy(float)
    f, h, ph, A = (d[c].to_numpy(float) for c in ("f", "h", "phiprime", "A0prime"))
    rhs = (d.a2.to_numpy(float) - 2*r*h*d.a6.to_numpy(float) - A*d.v4.to_numpy(float)/2) / (f*h*ph)
    err = d.e1.to_numpy(float) - rhs
    return {
        "max_abs": float(np.max(np.abs(err))),
        "max_scaled": float(np.max(np.abs(err)/np.maximum(1.0, np.abs(rhs)))),
    }


def main() -> int:
    current = prepare_current_pre_hessian()
    source = central_action_inputs(ROOT).sort_values("x").reset_index(drop=True)
    x = current.x.to_numpy(float)

    pre = {}
    for col in ("f3", "f2F", "f2", "f2X"):
        pre[col] = interp(source, col, x)

    variants: dict[str, pd.DataFrame] = {"current_post_A2": current.copy()}

    def make(name: str, cols: tuple[str, ...]):
        d = current.copy()
        for col in cols:
            d[col] = pre[col]
        variants[name] = d

    make("restore_f3_only", ("f3",))
    make("restore_f2F_only", ("f2F",))
    make("restore_f2_only", ("f2",))
    make("restore_f2X_only", ("f2X",))
    make("restore_f3_f2F", ("f3", "f2F"))
    make("restore_f2_f2X", ("f2", "f2X"))
    make("restore_f3_f2F_f2_f2X", ("f3", "f2F", "f2", "f2X"))

    reducer = module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    emitted = {name: emit_action(d) for name, d in variants.items()}
    result = {}
    for name, stream in emitted.items():
        result[name] = {
            "per_L": {str(L): metrics(stream, reducer, int(L)) for L in DEFAULT_L},
            "e1_identity": e1_identity_residual(stream),
        }

    # Directly test the analytic e1 <- v6 <- f3 channel at the common ghost location.
    u = emitted["current_post_A2"].u.to_numpy(float)
    i = int(np.argmin(np.abs(u - 0.7013062423913266)))
    witnesses = {}
    for name, stream in emitted.items():
        row = stream.iloc[i]
        witnesses[name] = {
            "u": float(row.u),
            "f3": float(variants[name].iloc[i].f3),
            "v6": float(row.v6),
            "a2": float(row.a2),
            "a6": float(row.a6),
            "v4": float(row.v4),
            "e1": float(row.e1),
        }

    # Does restoring f3 alone restore the pre-A2 high-L sign for every tested L?
    f3_repairs_all = all(
        result["restore_f3_only"]["per_L"][str(L)]["tail_min_K"] > 0
        for L in DEFAULT_L
    )
    full_restore_repairs_all = all(
        result["restore_f3_f2F_f2_f2X"]["per_L"][str(L)]["tail_min_K"] > 0
        for L in DEFAULT_L
    )

    current1000 = result["current_post_A2"]["per_L"]["1000"]["tail_min_K"]
    f3_1000 = result["restore_f3_only"]["per_L"]["1000"]["tail_min_K"]

    if f3_repairs_all:
        diagnosis = "F3_CHANNEL_IS_SUFFICIENT_TO_REMOVE_CURRENT_ALL_L_TAIL_GHOST_IN_COUNTERFACTUAL_ACTION_REPLAY"
        next_target = "DERIVE_AND_TEST_ACTION_CONSISTENT_F3_BRANCH_FROM_A2_JA_SYSTEM"
    elif full_restore_repairs_all:
        diagnosis = "DISTRIBUTED_A2_ACTION_CHANNELS_REQUIRED_FOR_TAIL_GHOST_REMOVAL"
        next_target = "DECOMPOSE_COUPLED_F3_F2F_F2_F2X_RESPONSE"
    else:
        diagnosis = "A2_ACTION_VARIABLE_CAUSE_NOT_ISOLATED_BY_TESTED_RESTORES"
        next_target = "EXPAND_ACTION_VARIABLE_DECOMPOSITION"

    payload = {
        "scope": "action-variable causality diagnostic; no fitting, no promoted repair",
        "analytic_dependency": {
            "e1_formula": "(a2 - 2*r*h*a6 - A0prime*v4/2)/(f*h*phiprime)",
            "v4_cancels_after_substitution": True,
            "v6_depends_on": "f3,f4,f4X,tf4 with fixed background",
            "current_chain_variable": "f3 is changed by Electric A2 resolve; f4/f4X/tf4 held fixed here",
        },
        "variants": result,
        "witness_near_common_zero": witnesses,
        "f3_restore_repairs_all_tested_L": f3_repairs_all,
        "full_A2_lower_restore_repairs_all_tested_L": full_restore_repairs_all,
        "L1000_current_tail_min_K": current1000,
        "L1000_restore_f3_tail_min_K": f3_1000,
        "diagnosis": diagnosis,
        "next_target": next_target,
        "guard": (
            "Restored variants are diagnostic counterfactuals.  A physical replacement must "
            "re-solve the same action/background equations and all dependent jets together."
        ),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    print(json.dumps({
        "diagnosis": diagnosis,
        "next_target": next_target,
        "f3_restore_repairs_all_tested_L": f3_repairs_all,
        "full_A2_lower_restore_repairs_all_tested_L": full_restore_repairs_all,
        "L1000_current_tail_min_K": current1000,
        "L1000_restore_f3_tail_min_K": f3_1000,
        "witness": witnesses,
    }, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
