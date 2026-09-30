#!/usr/bin/env python3
"""Test whether the u~0.7013 kinetic ghost coincides with an SSZ sign/branch flip.

This is diagnostic only. It:
1) locates K zero crossings for all required L;
2) searches the current same-action background/action profiles for sign flips,
   zeros, non-finite values, and sharp derivative/extremum changes near the
   common crossing;
3) decomposes the electric A2 ODE numerator term-by-term near the crossing;
4) reports whether any explicit sign flip is colocated with the ghost.

No parameter is changed and no member is re-frozen.
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
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central  # noqa: E402

OUT = ROOT / "data/generated/spectral/CURRENT_MEMBER_SIGN_BRANCH_DIAGNOSTIC.json"
WINDOW = (0.698, 0.704)


def first_zero(u: np.ndarray, y: np.ndarray) -> float | None:
    order = np.argsort(u)
    u = np.asarray(u, float)[order]
    y = np.asarray(y, float)[order]
    for i in range(len(u) - 1):
        if y[i] == 0:
            return float(u[i])
        if np.isfinite(y[i]) and np.isfinite(y[i+1]) and y[i] * y[i+1] < 0:
            t = -y[i] / (y[i+1] - y[i])
            return float(u[i] + t * (u[i+1] - u[i]))
    return None


def sign_flips(u: np.ndarray, y: np.ndarray, lo: float, hi: float) -> list[float]:
    order = np.argsort(u)
    u = np.asarray(u, float)[order]
    y = np.asarray(y, float)[order]
    out = []
    for i in range(len(u)-1):
        if u[i+1] < lo or u[i] > hi:
            continue
        if not (np.isfinite(y[i]) and np.isfinite(y[i+1])):
            continue
        if y[i] == 0:
            out.append(float(u[i]))
        elif y[i] * y[i+1] < 0:
            t = -y[i] / (y[i+1]-y[i])
            out.append(float(u[i] + t*(u[i+1]-u[i])))
    return out


def extrema_near(u: np.ndarray, y: np.ndarray, lo: float, hi: float) -> list[dict]:
    order = np.argsort(u)
    u = np.asarray(u, float)[order]
    y = np.asarray(y, float)[order]
    dy = np.gradient(y, u)
    out = []
    for i in range(len(u)-1):
        if u[i+1] < lo or u[i] > hi:
            continue
        if np.isfinite(dy[i]) and np.isfinite(dy[i+1]) and dy[i]*dy[i+1] < 0:
            out.append({"u": float((u[i]+u[i+1])/2), "kind": "d/du sign change"})
    return out[:20]


def main() -> int:
    build = build_onshell_central(ROOT)
    d = build.action.copy().sort_values("u").reset_index(drop=True)
    op = build.direct41.copy().sort_values("u").reset_index(drop=True)
    reducer = module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    zk = module("ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py")

    u = op.u.to_numpy(float)

    # 1) K zero crossings
    k_cross = {}
    k_at = {}
    for L in DEFAULT_L:
        a = reducer.canonical_audit(op.sort_values("x").reset_index(drop=True), int(L))
        # canonical_audit order follows x-sorted stream; map u accordingly
        ox = op.sort_values("x").reset_index(drop=True)
        uu = ox.u.to_numpy(float)
        K = np.asarray(a["K"], float)
        Ks = (K + K.swapaxes(1,2))/2
        kmin = np.linalg.eigvalsh(Ks)[:,0]
        z = first_zero(uu[(uu>=0.61)&(uu<0.71)], kmin[(uu>=0.61)&(uu<0.71)])
        k_cross[str(L)] = z
        # nearest sampled values around crossing
        if z is not None:
            ii = np.argsort(np.abs(uu-z))[:4]
            k_at[str(L)] = [
                {"u": float(uu[j]), "lambda_min_K": float(kmin[j])}
                for j in sorted(ii, key=lambda j: uu[j])
            ]

    finite_cross = np.array([x for x in k_cross.values() if x is not None], float)
    common = float(np.median(finite_cross))

    # 2) profile/action sign flips and extrema
    candidate_cols = [
        "f","h","phi","phiprime","A0prime","X","Fbg_action","Ybg_action",
        "f2","f2X","f2F","f2Y","f2phi","f2phiphi","f3","f3X",
        "f4","f4X","f4XX","tf4","G2XX_lift","G2Xphi_lift","G2phiphi_lift",
        "principal_delta_f2XX","principal_delta_f2XF","principal_delta_f2FF",
        "v6_A2_resolved",
    ]
    profiles = {}
    for col in candidate_cols:
        if col not in d.columns:
            continue
        y = d[col].to_numpy(float)
        m = (d.u.to_numpy(float)>=WINDOW[0])&(d.u.to_numpy(float)<=WINDOW[1])
        profiles[col] = {
            "min_window": float(np.nanmin(y[m])),
            "max_window": float(np.nanmax(y[m])),
            "sign_flips_window": sign_flips(d.u.to_numpy(float), y, *WINDOW),
            "extrema_window": extrema_near(d.u.to_numpy(float), y, *WINDOW),
            "nonfinite_window": int(np.sum(~np.isfinite(y[m]))),
        }

    # Also check direct 41 coefficients themselves for sign flips near crossing.
    slot_flips = {}
    for col in op.columns:
        try:
            y = op[col].to_numpy(float)
        except Exception:
            continue
        flips = sign_flips(u, y, *WINDOW)
        if flips:
            slot_flips[col] = flips

    # 3) Electric A2 ODE term decomposition on action grid.
    # Use action/direct41 quantities already produced by the exact current build.
    # Derivatives are recomputed with the same JET9D8 operator.
    a = d.sort_values("x").reset_index(drop=True)
    oo = op.sort_values("x").reset_index(drop=True)
    r = a.x.to_numpy(float)
    uu = a.u.to_numpy(float)
    f = a.f.to_numpy(float)
    h = a.h.to_numpy(float)
    ph = a.phiprime.to_numpy(float)
    aa = a.A0prime.to_numpy(float)
    vv = a.v6_A2_resolved.to_numpy(float)
    f3 = a.f3.to_numpy(float)
    f2F = a.f2F.to_numpy(float)
    f4 = a.f4.to_numpy(float)
    f4x = a.f4X.to_numpy(float)
    fp = zk.dr(r, f, 1, 9, 8)
    fpp = zk.dr(r, f, 2, 9, 8)
    hp = zk.dr(r, h, 1, 9, 8)
    app = zk.dr(r, aa, 1, 9, 8)

    # a4 is a direct coefficient from the same emitted stream.
    a4 = oo.a4.to_numpy(float)
    v10 = -np.sqrt(f*h)/(2*r) * (
        r*f2F + 2*h*ph*f3 + (h*fp/f)*(r*ph*f3 - 4*f4 + h*ph**2*f4x)
    )
    a7 = (1 - (4*h*aa**2/f)*f4)/(4*r**2*np.sqrt(f*h))

    terms = {
        "A4_geom_hp_h": a4 * (2*f**2*(r*hp + 2*h)),
        "A4_geom_fp2": a4 * (h*r**2*fp**2),
        "A4_geom_fpp_mix": -a4 * (f*r*(r*fp*hp + 2*h*(r*fpp + fp))),
        "electric_v10": -f*h*r**2 * (aa*(4*v10*aa)),
        "electric_v6_fp": -f*h*r**2 * (aa*(vv*fp)),
        "electric_v6_app": -f*h*r**2 * (f*vv*app),
        "a7_term": -f*h*r**2 * (8*a7*f**2),
    }
    den = f**2*h*r**2*aa
    num = np.sum(np.vstack(list(terms.values())), axis=0)
    rhs = num/den

    a2 = {
        "denominator_sign_flips_window": sign_flips(uu, den, *WINDOW),
        "numerator_sign_flips_window": sign_flips(uu, num, *WINDOW),
        "rhs_sign_flips_window": sign_flips(uu, rhs, *WINDOW),
        "term_sign_flips_window": {
            name: sign_flips(uu, arr, *WINDOW) for name, arr in terms.items()
        },
        "samples_around_common_K_zero": [],
    }
    near = np.argsort(np.abs(uu-common))[:8]
    for j in sorted(near, key=lambda j: uu[j]):
        row = {
            "u": float(uu[j]), "r": float(r[j]), "num": float(num[j]),
            "den": float(den[j]), "dv6_dr_rhs": float(rhs[j]),
            "v6": float(vv[j]), "f": float(f[j]), "h": float(h[j]),
            "phiprime": float(ph[j]), "A0prime": float(aa[j]), "X": float(a.X.iloc[j]),
        }
        row.update({name: float(arr[j]) for name, arr in terms.items()})
        a2["samples_around_common_K_zero"].append(row)

    explicit_profile_flips = {
        k:v["sign_flips_window"] for k,v in profiles.items() if v["sign_flips_window"]
    }
    a2_flips = {
        "denominator": a2["denominator_sign_flips_window"],
        "numerator": a2["numerator_sign_flips_window"],
        "rhs": a2["rhs_sign_flips_window"],
        "terms": {k:v for k,v in a2["term_sign_flips_window"].items() if v},
    }

    tol = 5e-4
    colocated = []
    for name, flips in explicit_profile_flips.items():
        for z in flips:
            if abs(z-common) <= tol:
                colocated.append({"kind":"profile", "name":name, "u":z})
    for name, flips in slot_flips.items():
        for z in flips:
            if abs(z-common) <= tol:
                colocated.append({"kind":"41slot", "name":name, "u":z})
    for name, flips in [("A2_den",a2["denominator_sign_flips_window"]),("A2_num",a2["numerator_sign_flips_window"]),("A2_rhs",a2["rhs_sign_flips_window"])]:
        for z in flips:
            if abs(z-common) <= tol:
                colocated.append({"kind":"A2", "name":name, "u":z})
    for name, flips in a2["term_sign_flips_window"].items():
        for z in flips:
            if abs(z-common) <= tol:
                colocated.append({"kind":"A2_term", "name":name, "u":z})

    payload = {
        "scope": "test explicit/implicit sign or branch change near the common kinetic zero; no modification",
        "window_u": list(WINDOW),
        "K_zero_crossings": k_cross,
        "common_K_zero_median_u": common,
        "K_samples": k_at,
        "profile_action_diagnostics": profiles,
        "direct_41_slot_sign_flips_window": slot_flips,
        "A2_ODE_decomposition": a2,
        "explicit_profile_sign_flips_window": explicit_profile_flips,
        "A2_sign_flip_summary": a2_flips,
        "colocated_with_K_zero_within_5e-4": colocated,
        "verdict": (
            "COLOCATED_SIGN_OR_BRANCH_EVENT_FOUND"
            if colocated else
            "NO_EXPLICIT_SIGN_FLIP_COLOCATED_WITH_K_ZERO"
        ),
        "interpretation_guard": (
            "Absence of an explicit sign flip does not exclude a smooth branch inconsistency, "
            "conditioning failure, or wrong upstream action branch. It only tests direct zero/sign events."
        ),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2, allow_nan=False)+"\n")
    print(json.dumps(payload, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
