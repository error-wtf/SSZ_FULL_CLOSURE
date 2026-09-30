#!/usr/bin/env python3
"""Localize the inner-tail radial-characteristic failure of the current member.

The full-domain audit found a second pathology distinct from the outer kinetic
ghost: for u<=0.62 the kinetic matrix stays positive, but the minimum
characteristic eigenvalue of K^{-1/2} G K^{-1/2} becomes negative.

This audit:
* locates the first c_r^2 zero nearest the production boundary,
* tracks the failing eigenvector and its physical-basis composition,
* compares construction stages S1/S2/S3,
* compares the archived pre-A2 source stream where reducer-compatible,
* tests whether the same characteristic branch is shared across L.

No repair or parameter search is performed.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ssz_p5.config import DEFAULT_L, SLOT_NAMES  # noqa:E402
from ssz_p5.numerics import module  # noqa:E402
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central  # noqa:E402
from ssz_p5.production.full_action_lower import complete_total_action_jets, emit_lower_slots  # noqa:E402

OUT = ROOT / "data/generated/spectral/CURRENT_MEMBER_INNER_RADIAL_ORIGIN_AUDIT.json"
RAW = ROOT / "archive/full_working_snapshot/ssz_p5_integrable_svt_unreduced_even_coefficients_2026-09-12.csv"


def first_zero_near_boundary(u, y, boundary=0.62):
    order = np.argsort(u)
    uu = np.asarray(u)[order]
    yy = np.asarray(y)[order]
    hits = []
    for i in range(len(uu)-1):
        if not np.isfinite(yy[i:i+2]).all():
            continue
        if yy[i] == 0:
            hits.append(float(uu[i]))
        elif yy[i]*yy[i+1] < 0:
            t = -yy[i]/(yy[i+1]-yy[i])
            hits.append(float(uu[i] + t*(uu[i+1]-uu[i])))
    return min(hits, key=lambda z: abs(z-boundary)) if hits else None


def emit_from_action(action):
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
        - zk.dr(r, out.A0prime.to_numpy(float)*out.v4.to_numpy(float)/2, 1, 9, 8)
        + out.A0prime.to_numpy(float)*out.v5.to_numpy(float)/2
    )
    return out


def audit_stream(d, red, L):
    d = d.sort_values("x").reset_index(drop=True)
    u = d.u.to_numpy(float)
    a = red.canonical_audit(d, int(L))
    K = np.asarray(a["K"], float)
    G = np.asarray(a["G"], float)
    Ks = (K + K.swapaxes(1,2))/2
    Gs = (G + G.swapaxes(1,2))/2

    crmin = np.full(len(d), np.nan)
    comps = np.full((len(d),3), np.nan)
    kmin = np.linalg.eigvalsh(Ks)[:,0]
    for i in range(len(d)):
        w,U = np.linalg.eigh(Ks[i])
        if np.min(w) <= 0:
            continue
        inv = U @ np.diag(1/np.sqrt(w)) @ U.T
        C = (inv @ Gs[i] @ inv)
        C = (C + C.T)/2
        ev,V = np.linalg.eigh(C)
        crmin[i] = ev[0]
        # map whitened characteristic eigenvector back to physical Y basis
        y = inv @ V[:,0]
        norm = np.linalg.norm(y)
        if norm:
            comps[i] = np.abs(y/norm)**2

    inner = u <= 0.62
    j = np.flatnonzero(inner)[np.nanargmin(crmin[inner])]
    z = first_zero_near_boundary(u[u<=0.625], crmin[u<=0.625])
    return {
        "min_K_inner": float(np.min(kmin[inner])),
        "min_cr2_inner": float(crmin[j]),
        "u_at_min_cr2": float(u[j]),
        "first_cr2_zero_near_0p62": z,
        "negative_cr2_rows_inner": int(np.sum(crmin[inner] < 0)),
        "failing_eigenvector_component_power": {
            "psi": float(comps[j,0]),
            "dphi": float(comps[j,1]),
            "V": float(comps[j,2]),
        },
    }


def raw_on_grid(raw, cur):
    out = cur.copy()
    xr = raw.x.to_numpy(float)
    order = np.argsort(xr)
    xs = xr[order]
    for c in cur.columns:
        if c in raw.columns:
            try:
                arr = raw[c].to_numpy(float)[order]
                if np.isfinite(arr).all():
                    out[c] = np.interp(cur.x.to_numpy(float), xs, arr)
            except Exception:
                pass
    for c in ("x","u","f","h","phi","phiprime","A0prime"):
        if c in cur.columns:
            out[c] = cur[c].to_numpy()
    return out


def main():
    build = build_onshell_central(ROOT)
    action = build.action.copy().sort_values("x").reset_index(drop=True)
    pre_h = action.copy()
    for target, delta in (
        ("f2XX","principal_delta_f2XX"),
        ("f2XF","principal_delta_f2XF"),
        ("f2FF","principal_delta_f2FF"),
    ):
        pre_h[target] = pre_h[target].to_numpy(float) - pre_h[delta].to_numpy(float)

    stages = {
        "S1_pre_hessian": emit_from_action(pre_h),
        "S2_pre_G2XX": build.pre_lift41.copy(),
        "S3_locked_current": build.direct41.copy(),
    }

    raw = pd.read_csv(RAW)
    required = {"x","u"} | set(SLOT_NAMES)
    raw_compatible = required.issubset(raw.columns)
    if raw_compatible:
        stages["S0_archival_source_on_current_grid"] = raw_on_grid(raw, build.direct41.copy())

    red = module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    results = {
        name: {str(L): audit_stream(stream, red, int(L)) for L in DEFAULT_L}
        for name, stream in stages.items()
    }

    def fails(name):
        return any(results[name][str(L)]["min_cr2_inner"] < 0 for L in DEFAULT_L)

    flags = {name: fails(name) for name in stages}

    if raw_compatible and flags.get("S0_archival_source_on_current_grid", False):
        diagnosis = "INNER_RADIAL_FAILURE_ALREADY_PRESENT_IN_ARCHIVAL_PRE_A2_SOURCE"
        next_target = "TRACE_ARCHIVAL_PRINCIPAL_G_CHANNEL"
    elif flags["S1_pre_hessian"]:
        diagnosis = "INNER_RADIAL_FAILURE_PRESENT_BEFORE_HESSIAN_RESTORE_IN_CURRENT_A2_BRANCH"
        next_target = "TRACE_A2_OR_EARLIER_G_CHANNEL"
    elif flags["S2_pre_G2XX"]:
        diagnosis = "INNER_RADIAL_FAILURE_FIRST_APPEARS_WITH_HESSIAN_RESTORE"
        next_target = "TRACE_HESSIAN_TO_G_RESPONSE"
    elif flags["S3_locked_current"]:
        diagnosis = "INNER_RADIAL_FAILURE_FIRST_APPEARS_WITH_G2XX_LIFT"
        next_target = "TRACE_G2XX_TO_G_RESPONSE"
    else:
        diagnosis = "NO_INNER_RADIAL_FAILURE_REPRODUCED"
        next_target = "CHECK_AUDIT_PROVENANCE"

    # quantify cross-L commonality of zero location and eigenvector composition
    locked = results["S3_locked_current"]
    zeros = [locked[str(L)]["first_cr2_zero_near_0p62"] for L in DEFAULT_L]
    zeros_f = [z for z in zeros if z is not None]
    comp = np.array([
        list(locked[str(L)]["failing_eigenvector_component_power"].values())
        for L in DEFAULT_L
    ], float)

    report = {
        "scope": "inner-tail radial characteristic origin audit; no repair/no fitting",
        "raw_source_reducer_compatible": raw_compatible,
        "stage_fail_flags": flags,
        "stages": results,
        "locked_cross_L": {
            "zero_u_mean": float(np.mean(zeros_f)) if zeros_f else None,
            "zero_u_range": float(np.ptp(zeros_f)) if len(zeros_f)>1 else 0.0 if zeros_f else None,
            "mean_failing_component_power": {
                "psi": float(np.mean(comp[:,0])),
                "dphi": float(np.mean(comp[:,1])),
                "V": float(np.mean(comp[:,2])),
            },
            "component_power_std": {
                "psi": float(np.std(comp[:,0])),
                "dphi": float(np.std(comp[:,1])),
                "V": float(np.std(comp[:,2])),
            },
        },
        "diagnosis": diagnosis,
        "next_target": next_target,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    print(json.dumps({
        "diagnosis":diagnosis,
        "next_target":next_target,
        "stage_fail_flags":flags,
        "locked_cross_L":report["locked_cross_L"],
    }, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
