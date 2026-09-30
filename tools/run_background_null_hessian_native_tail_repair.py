#!/usr/bin/env python3
"""Action-consistent background-null Hessian repair trial for the native tails.

This uses the exact smooth f2(phi,X,F,Y) background-null Hessian manifold:

    H @ (phi',X',F',Y') = 0

implemented by holonomic_hessian_y. The deformation changes only second
action jets, leaves the background value/first jets unchanged, and therefore is
an admissible same-background action freedom rather than an arbitrary 41-slot
patch.

We use twelve smooth controls: six transverse Hessian directions multiplied by
one inner-tail envelope and one outer-tail envelope. The exact Appendix-A
41-stream response is linear in these controls.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ssz_p5.config import DEFAULT_L  # noqa:E402
from ssz_p5.numerics import module  # noqa:E402
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central  # noqa:E402
from ssz_p5.production.holonomic_hessian_y import (  # noqa:E402
    TRANSVERSE,
    complete_hessian,
    response_from_hessian,
    chain_residual,
)

OUT = ROOT / "data/generated/spectral/BACKGROUND_NULL_HESSIAN_NATIVE_TAIL_REPAIR.json"
STREAM_OUT = ROOT / "data/generated/spectral/BACKGROUND_NULL_HESSIAN_NATIVE_TAIL_REPAIR_DIRECT41.csv"

LS_OPT = (6, 42, 1000)
LS_VERIFY = tuple(int(x) for x in DEFAULT_L)
EPS = 2e-4
MAX_IT = 10
TRUST = 0.4


def cinf_interval(u, lo, hi, edge=7e-4):
    u = np.asarray(u, float)
    out = np.zeros_like(u)
    core = (u >= lo + edge) & (u <= hi - edge)
    out[core] = 1.0
    left = (u > lo) & (u < lo + edge)
    right = (u > hi - edge) & (u < hi)
    if np.any(left):
        t = (u[left] - lo) / edge
        out[left] = t**3 * (10 - 15*t + 6*t*t)
    if np.any(right):
        t = (hi - u[right]) / edge
        out[right] = t**3 * (10 - 15*t + 6*t*t)
    return out


def build_response_columns(action, direct):
    u = action.u.to_numpy(float)
    envs = {
        "inner": cinf_interval(u, 0.6100, 0.6200, edge=9e-4),
        "outer": cinf_interval(u, 0.7000, float(u.max()) + 1e-12, edge=7e-4),
    }
    cols = []
    meta = []
    for region, env in envs.items():
        for j, name in enumerate(TRANSVERSE):
            q = np.zeros((len(action), 6), float)
            q[:, j] = env
            H, _completed = complete_hessian(action, q)
            cres, cnorm = chain_residual(action, H)
            resp = response_from_hessian(action, H)
            delta = {slot: np.zeros(len(action), float) for slot in direct.columns if slot not in ("u", "x")}
            for slot in ("v5", "c3", "e3", "v1", "v4", "c2"):
                delta[slot] = resp[slot].to_numpy(float)
            cols.append(delta)
            meta.append({
                "region": region,
                "transverse": name,
                "max_abs_chain_residual": float(np.max(np.abs(cres))),
                "max_normalized_chain_residual": float(np.max(cnorm)),
                "support_nonzero_rows": int(np.sum(env != 0)),
            })
    return cols, meta


def stream_from(base, cols, p):
    d = base.copy()
    for slot in ("v5", "c3", "e3", "v1", "v4", "c2"):
        arr = d[slot].to_numpy(float).copy()
        for j, a in enumerate(p):
            if a:
                arr += float(a) * cols[j][slot]
        d[slot] = arr
    return d


def point_metrics(stream, reducer, L, idx):
    a = reducer.canonical_audit(stream, int(L))
    K = np.asarray(a["K"], float)
    G = np.asarray(a["G"], float)
    out = []
    for i in idx:
        Ks = (K[i] + K[i].T) / 2
        Gs = (G[i] + G[i].T) / 2
        w, U = np.linalg.eigh(Ks)
        kmin = float(w[0])
        cr2 = np.nan
        if kmin > 0:
            inv = U @ np.diag(1 / np.sqrt(w)) @ U.T
            C = inv @ Gs @ inv
            cr2 = float(np.linalg.eigvalsh((C + C.T) / 2)[0])
        out.append((kmin, cr2))
    return out


def objective_vector(stream, reducer, u):
    inner_all = np.flatnonzero(u <= 0.62)
    outer_all = np.flatnonzero(u >= 0.70)
    inner = inner_all[np.unique(np.linspace(0, len(inner_all)-1, 9).round().astype(int))]
    outer = outer_all[np.unique(np.linspace(0, len(outer_all)-1, 13).round().astype(int))]
    vals = []
    tags = []
    for L in LS_OPT:
        met = point_metrics(stream, reducer, L, outer)
        for i, (k, _c) in zip(outer, met):
            vals.append(k); tags.append(("outer_K", L, int(i)))
    for L in LS_OPT:
        met = point_metrics(stream, reducer, L, inner)
        for i, (k, c) in zip(inner, met):
            vals.append(c); tags.append(("inner_cr2", L, int(i)))
            vals.append(k); tags.append(("inner_K_guard", L, int(i)))
    return np.asarray(vals, float), tags


def worst_by_kind(vals, tags, u):
    out = {}
    for kind in ("outer_K", "inner_cr2", "inner_K_guard"):
        ids = [i for i,t in enumerate(tags) if t[0] == kind]
        q = ids[int(np.argmin(vals[ids]))]
        tag = tags[q]
        out[kind] = {"value": float(vals[q]), "L": int(tag[1]), "u": float(u[tag[2]]), "row": int(tag[2])}
    return out


def verify_full(stream, reducer, u):
    inner = np.flatnonzero(u <= 0.62)
    outer = np.flatnonzero(u >= 0.70)
    prod = np.flatnonzero((u > 0.62) & (u < 0.70))
    report = {}
    for L in LS_VERIFY:
        a = reducer.canonical_audit(stream, int(L))
        K = np.asarray(a["K"], float); G = np.asarray(a["G"], float)
        Ks = (K + K.transpose(0,2,1))/2
        Gs = (G + G.transpose(0,2,1))/2
        ke = np.linalg.eigvalsh(Ks)[:,0]
        cr = np.full(len(u), np.nan)
        for i in np.flatnonzero(ke > 0):
            w,U = np.linalg.eigh(Ks[i])
            inv = U @ np.diag(1/np.sqrt(w)) @ U.T
            C = inv @ Gs[i] @ inv
            cr[i] = np.linalg.eigvalsh((C+C.T)/2)[0]
        report[str(L)] = {
            "production_min_K": float(np.min(ke[prod])),
            "outer_min_K": float(np.min(ke[outer])),
            "outer_negative_K_rows": int(np.sum(ke[outer] <= 0)),
            "inner_min_K": float(np.min(ke[inner])),
            "inner_min_cr2": float(np.nanmin(cr[inner])),
            "inner_negative_cr2_rows": int(np.sum(cr[inner] < 0)),
            "full_min_cr2_where_K_positive": float(np.nanmin(cr)),
            "full_negative_cr2_rows_where_K_positive": int(np.sum(cr < 0)),
        }
    return report


def main():
    build = build_onshell_central(ROOT)
    action = build.action.sort_values("x").reset_index(drop=True)
    base = build.direct41.sort_values("x").reset_index(drop=True)
    if not np.allclose(action.x, base.x, rtol=0, atol=0):
        raise RuntimeError("action/direct grid mismatch")
    u = base.u.to_numpy(float)
    reducer = module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")

    cols, meta = build_response_columns(action, base)
    n = len(cols)
    p = np.zeros(n)
    stream = stream_from(base, cols, p)
    vals, tags = objective_vector(stream, reducer, u)
    initial = worst_by_kind(vals, tags, u)
    hist = []

    for it in range(MAX_IT):
        G = np.zeros((len(vals), n))
        for j in range(n):
            pp = p.copy(); pp[j] += EPS
            vv, _ = objective_vector(stream_from(base, cols, pp), reducer, u)
            vv = np.where(np.isfinite(vv), vv, -1e6)
            G[:,j] = (vv - np.where(np.isfinite(vals), vals, -1e6)) / EPS

        c = np.zeros(n+1); c[-1] = -1
        Aub = np.column_stack([-G, np.ones(len(vals))])
        bub = np.where(np.isfinite(vals), vals, -1e6)
        sol = linprog(c, A_ub=Aub, b_ub=bub,
                      bounds=[(-TRUST, TRUST)]*n + [(None,None)], method="highs")
        if not sol.success:
            hist.append({"iteration":it,"status":"LP_FAIL","message":sol.message})
            break
        dp = sol.x[:-1]
        accepted = False
        best = None
        oldw = worst_by_kind(vals,tags,u)
        oldscore = min(oldw["outer_K"]["value"], oldw["inner_cr2"]["value"])
        for frac in (1.0, .5, .25, .125, .0625):
            pt = p + frac*dp
            st = stream_from(base, cols, pt)
            vt, tt = objective_vector(st, reducer, u)
            if not np.all(np.isfinite(vt)):
                continue
            w = worst_by_kind(vt, tt, u)
            score = min(w["outer_K"]["value"], w["inner_cr2"]["value"])
            if w["inner_K_guard"]["value"] <= 0:
                continue
            rec = {"iteration":it,"fraction":frac,"score":float(score),"worst":w,
                   "step_max_abs":float(np.max(np.abs(frac*dp)))}
            if best is None or score > best[0]:
                best = (score, pt, st, vt, tt, rec)
            if score > oldscore + 1e-5:
                p,stream,vals,tags = pt,st,vt,tt
                rec["status"]="ACCEPT"; hist.append(rec); accepted=True
                break
        if not accepted:
            if best and best[0] > oldscore:
                _,p,stream,vals,tags,rec=best
                rec["status"]="ACCEPT_BEST"; hist.append(rec); accepted=True
            else:
                hist.append({"iteration":it,"status":"STALL"})
                break
        w = worst_by_kind(vals,tags,u)
        if w["outer_K"]["value"] > 0 and w["inner_cr2"]["value"] > 0:
            break

    final_sample = worst_by_kind(vals, tags, u)
    full = verify_full(stream, reducer, u)
    all_outer = all(v["outer_min_K"] > 0 for v in full.values())
    all_inner = all(v["inner_min_cr2"] > 0 for v in full.values())
    all_prod = all(v["production_min_K"] > 0 for v in full.values())

    STREAM_OUT.parent.mkdir(parents=True,exist_ok=True)
    stream.to_csv(STREAM_OUT,index=False)

    if all_outer and all_inner and all_prod:
        status="BACKGROUND_NULL_HESSIAN_NATIVE_TAIL_REPAIR_PASS"
        nxt="RUN_FULL_NATIVE_RADIAL_AND_SPECTROSCOPY_PRECERTIFICATE"
    elif min(final_sample["outer_K"]["value"],final_sample["inner_cr2"]["value"]) > min(initial["outer_K"]["value"],initial["inner_cr2"]["value"]):
        status="BACKGROUND_NULL_HESSIAN_NATIVE_TAIL_REPAIR_IMPROVED_NOT_CLOSED"
        nxt="EXPAND_HESSIAN_BASIS_TO_MULTIPLE_LOCAL_ENVELOPES"
    else:
        status="BACKGROUND_NULL_HESSIAN_NATIVE_TAIL_REPAIR_STALLED"
        nxt="USE_MULTIPLE_LOCAL_HESSIAN_ENVELOPES_OR_ADDITIONAL_BACKGROUND_NULL_SECTORS"

    report={
        "status":status,
        "next_target":nxt,
        "control_semantics":"exact background-null smooth f2 Hessian; action-consistent Appendix-A response",
        "response_columns":meta,
        "initial_sample_worst":initial,
        "iterations":hist,
        "final_sample_worst":final_sample,
        "full_verification":full,
        "all_outer_K_pass":all_outer,
        "all_inner_radial_pass":all_inner,
        "all_production_K_pass":all_prod,
        "parameters":[float(x) for x in p],
        "stream_csv":str(STREAM_OUT.relative_to(ROOT)),
    }
    OUT.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps(report,indent=2,allow_nan=False))
    return 0

if __name__=="__main__":
    raise SystemExit(main())
