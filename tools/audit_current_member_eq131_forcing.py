#!/usr/bin/env python3
"""Test whether the new high-L ghost is specifically forced by the Eq.131/A2 ODE.

Constructs a diagnostic branch using the archived/source f3 profile, then
recomputes f2F from JA=0 and f2,f2X from E00=E11=0 on the same current
background.  This preserves the algebraic background equations but does NOT
impose Eq.131.  We then compare kinetic health and the undivided Eq.131
residual against the current A2-resolved member.

No fit, optimization, re-freeze, or promotion is performed.
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
from ssz_p5.production.central_action import central_action_inputs, background_residuals  # noqa: E402
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central  # noqa: E402
from ssz_p5.production.full_action_lower import complete_total_action_jets, emit_lower_slots  # noqa: E402

OUT = ROOT / "data/generated/spectral/CURRENT_MEMBER_EQ131_FORCING_AUDIT.json"


def emit_action(d: pd.DataFrame) -> pd.DataFrame:
    zk = module("ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py")
    d = complete_total_action_jets(d.copy().sort_values("x").reset_index(drop=True))
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
    G = np.asarray(a["G"], float)
    Ks = (K + K.swapaxes(1, 2))/2
    Gs = (G + G.swapaxes(1, 2))/2
    ke = np.linalg.eigvalsh(Ks)[:, 0]
    tail = (u >= 0.70) & (u < 0.71)
    prod = (u > 0.62) & (u < 0.70)
    cr2 = np.full(len(d), np.nan)
    for i in np.flatnonzero(np.linalg.eigvalsh(Ks)[:,0] > 0):
        w,U = np.linalg.eigh(Ks[i])
        inv = U @ np.diag(1/np.sqrt(w)) @ U.T
        C = inv @ Gs[i] @ inv
        cr2[i] = np.linalg.eigvalsh((C+C.T)/2)[0]
    return {
        "production_min_K": float(np.min(ke[prod])),
        "tail_min_K": float(np.min(ke[tail])),
        "tail_negative_rows": int(np.sum(ke[tail] <= 0)),
        "tail_min_cr2_where_defined": (
            float(np.nanmin(cr2[tail])) if np.any(np.isfinite(cr2[tail])) else None
        ),
    }


def maxabs(res: dict[str, np.ndarray], mask: np.ndarray) -> dict[str, float]:
    return {k: float(np.max(np.abs(np.asarray(v)[mask]))) for k,v in res.items()}


def main() -> int:
    zk = module("ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py")
    reducer = module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")

    build = build_onshell_central(ROOT)
    cur = build.action.copy().sort_values("x").reset_index(drop=True)
    # Undo later Hessian restore so we stay at the A2 stage.
    for target, delta in (
        ("f2XX", "principal_delta_f2XX"),
        ("f2XF", "principal_delta_f2XF"),
        ("f2FF", "principal_delta_f2FF"),
    ):
        cur[target] = cur[target].to_numpy(float) - cur[delta].to_numpy(float)

    src = central_action_inputs(ROOT).sort_values("x").reset_index(drop=True)
    r = cur.x.to_numpy(float)
    u = cur.u.to_numpy(float)
    f,h,ph,ap = (cur[c].to_numpy(float) for c in ("f","h","phiprime","A0prime"))
    fp = zk.dr(r,f,1,9,8)
    hp = zk.dr(r,h,1,9,8)

    def interp(col: str) -> np.ndarray:
        return np.interp(r, src.x.to_numpy(float), src[col].to_numpy(float))

    alt = cur.copy()
    alt["f3"] = interp("f3")

    # Keep the current fixed f4-sector. Recompute JA=0 -> f2F.
    f3 = alt.f3.to_numpy(float)
    f4 = alt.f4.to_numpy(float)
    f4x = alt.f4X.to_numpy(float)
    tf4 = alt.tf4.to_numpy(float)
    f2y = alt.f2Y.to_numpy(float)
    alt["f2F"] = (
        2*h*ph**2*f2y
        - (
            4*r*h*ph*f3
            + 8*(1-h)*f4
            + 2*h**2*ph**2*(f4x + 2*tf4)
        )/r**2
    )

    # Solve E00=0 algebraically for f2.
    f2F = alt.f2F.to_numpy(float)
    f3x = alt.f3X.to_numpy(float)
    f4xx = alt.f4XX.to_numpy(float)
    tf3 = alt.tf3.to_numpy(float)
    tf4 = alt.tf4.to_numpy(float)

    rest00 = (
        f*(1-h)
        - r**2*h*ap**2*(f2F - 2*h*ph**2*f2y)
        - 2*r*h**2*ph*ap**2*f3
        + h*ap**2*(4*(h-1)*f4 - h**2*ph**2*(f4x + 2*tf4))
    )
    alt["f2"] = (r*f*hp - rest00)/(r**2*f)

    # Solve E11=0 algebraically for f2X.
    f2 = alt.f2.to_numpy(float)
    const11 = (
        f*(1-h)
        + r**2*(f*f2 - h*ap**2*(f2F - 4*h*ph**2*f2y))
        - 2*r*h**2*ph*ap**2*(3*f3 - h*ph**2*f3x)
        + h*ap**2*(
            4*(3*h-1)*f4
            - h*(9*h-4)*ph**2*f4x
            + h**3*ph**4*f4xx
            - 10*h**2*ph**2*tf4
        )
    )
    alt["f2X"] = (r*h*fp - const11)/(r**2*f*h*ph**2)

    # Refresh f2phi from total on-curve chain rule before completing mixed jets.
    X = alt.X.to_numpy(float)
    F = h*ap**2/(2*f)
    Xp = zk.dr(r,X,1,9,8)
    Fp = zk.dr(r,F,1,9,8)
    Y = 4*X*F
    Yp = zk.dr(r,Y,1,9,8)
    alt["f2phi"] = (
        zk.dr(r,alt.f2.to_numpy(float),1,9,8)
        - alt.f2X.to_numpy(float)*Xp
        - alt.f2F.to_numpy(float)*Fp
        - alt.f2Y.to_numpy(float)*Yp
    )/ph

    current_stream = emit_action(cur)
    alt_stream = emit_action(alt)

    region = (u >= 0.61) & (u < 0.71)
    prod = (u > 0.62) & (u < 0.70)
    tail = (u >= 0.70) & (u < 0.71)

    bg_cur = background_residuals(complete_total_action_jets(cur.copy()))
    bg_alt = background_residuals(complete_total_action_jets(alt.copy()))

    perL = {
        "current_A2": {str(L): metrics(current_stream,reducer,int(L)) for L in DEFAULT_L},
        "source_f3_plus_JA_E00_E11_resolved": {str(L): metrics(alt_stream,reducer,int(L)) for L in DEFAULT_L},
    }

    # High-L were healthy before A2, so use L>=110 to isolate the new ghost.
    highL = (110,420,1000)
    repairs_new_highL = all(
        perL["source_f3_plus_JA_E00_E11_resolved"][str(L)]["tail_min_K"] > 0
        for L in highL
    )

    alt_bg = maxabs(bg_alt, region)
    cur_bg = maxabs(bg_cur, region)
    alt_eq131_prod = float(np.max(np.abs(np.asarray(bg_alt["eq131_residual"])[prod])))
    alt_eq131_tail = float(np.max(np.abs(np.asarray(bg_alt["eq131_residual"])[tail])))
    cur_eq131_tail = float(np.max(np.abs(np.asarray(bg_cur["eq131_residual"])[tail])))

    if repairs_new_highL and alt_eq131_tail > max(1e-8, 100*cur_eq131_tail):
        diagnosis = "NEW_HIGH_L_GHOST_IS_TIED_TO_EQ131_A2_ENFORCEMENT_NOT_JA_E00_E11"
        next_target = "TEST_EQ131_BOUNDARY_BRANCH_AND_INITIAL_CONDITION_DEPENDENCE"
    elif repairs_new_highL:
        diagnosis = "SOURCE_F3_ALGEBRAIC_BRANCH_REPAIRS_HIGH_L_WITHOUT_LARGE_EQ131_PENALTY"
        next_target = "CHECK_FULL_SAME_ACTION_CONSISTENCY_AND_BOUNDARY_PROVENANCE"
    else:
        diagnosis = "ALGEBRAIC_SOURCE_F3_BRANCH_DOES_NOT_REMOVE_NEW_HIGH_L_GHOST"
        next_target = "EXPAND_COUPLED_ACTION_SECTOR"

    payload = {
        "scope": "Eq131 forcing audit; diagnostic only; no fit/no promotion",
        "construction": "source f3; recompute f2F by JA=0, f2 by E00=0, f2X by E11=0; fixed geometry/f4 sector",
        "per_L": perL,
        "background_residual_max_abs_region": {
            "current_A2": cur_bg,
            "source_f3_algebraic_branch": alt_bg,
        },
        "eq131_residual_max_abs": {
            "current_tail": cur_eq131_tail,
            "source_f3_branch_production": alt_eq131_prod,
            "source_f3_branch_tail": alt_eq131_tail,
        },
        "repairs_new_highL_110_420_1000": repairs_new_highL,
        "diagnosis": diagnosis,
        "next_target": next_target,
        "guard": "Alternative branch is algebraically on-shell for JA/E00/E11 only; Eq131 is intentionally left as the discriminator.",
    }
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(payload,indent=2,allow_nan=False)+"\n")
    print(json.dumps(payload,indent=2,allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
