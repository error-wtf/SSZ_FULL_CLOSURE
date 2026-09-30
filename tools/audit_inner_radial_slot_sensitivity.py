#!/usr/bin/env python3
"""41-slot sensitivity audit for the common inner radial c_r^2 crossing.

The inner-tail audit localized the current-member radial characteristic zero to
u~=0.61902123 and found an almost L-independent failing eigenvector that is
~99.6% psi, ~0.4% dphi, negligible V.  The failure is already present in the
archival pre-A2 source, so the next diagnostic cut is the reducer-level channel.

This script perturbs each of the 41 emitted coefficient slots with the same
compact smooth bump centered on the crossing and measures the finite-difference
response of
    lambda_min(K)
and
    c_r,min^2 = lambda_min(K^{-1/2} G K^{-1/2})
at the center.

The perturbations are diagnostic and are NOT action-consistent candidate
members.  They only identify the coefficient channels to trace upstream.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))

from ssz_p5.config import SLOT_NAMES  # noqa:E402
from ssz_p5.numerics import module  # noqa:E402
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central  # noqa:E402

OUT=ROOT/"data/generated/spectral/CURRENT_MEMBER_INNER_RADIAL_SLOT_SENSITIVITY.json"
LS=(6,110,1000)
U0=0.6190212327
WIDTH=0.0012
REL_EPS=2e-6


def bump(u):
    z=(np.asarray(u,float)-U0)/WIDTH
    q=np.zeros_like(z)
    m=np.abs(z)<1
    q[m]=np.exp(-1/(1-z[m]**2))/np.exp(-1)
    return q


def at_center(d,red,L):
    u=d.u.to_numpy(float)
    i=int(np.argmin(np.abs(u-U0)))
    a=red.canonical_audit(d,int(L))
    K=np.asarray(a["K"],float);G=np.asarray(a["G"],float)
    Ks=(K[i]+K[i].T)/2;Gs=(G[i]+G[i].T)/2
    kw,U=np.linalg.eigh(Ks)
    kmin=float(kw[0])
    if kmin<=0:
        return {"u":float(u[i]),"kmin":kmin,"cr2":None}
    inv=U@np.diag(1/np.sqrt(kw))@U.T
    C=inv@Gs@inv;C=(C+C.T)/2
    ev,V=np.linalg.eigh(C)
    y=inv@V[:,0];p=np.abs(y/np.linalg.norm(y))**2
    return {
        "u":float(u[i]),"kmin":kmin,"cr2":float(ev[0]),
        "component_power":{"psi":float(p[0]),"dphi":float(p[1]),"V":float(p[2])},
    }


def main():
    full=build_onshell_central(ROOT).direct41.sort_values("x").reset_index(drop=True)
    # generous local window for profile-aware 9-point derivatives
    d=full[(full.u>=0.6135)&(full.u<=0.6245)].sort_values("x").reset_index(drop=True)
    u=d.u.to_numpy(float);sh=bump(u)
    red=module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")

    base={str(L):at_center(d,red,L) for L in LS}
    per_slot={}
    for slot in SLOT_NAMES:
        arr=d[slot].to_numpy(float)
        ic=int(np.argmin(np.abs(u-U0)))
        scale=max(1.0,abs(float(arr[ic])))
        eps=REL_EPS*scale
        p=d.copy();p[slot]=arr+eps*sh
        m=d.copy();m[slot]=arr-eps*sh
        rows={}
        for L in LS:
            bp=at_center(p,red,L);bm=at_center(m,red,L)
            if bp["cr2"] is None or bm["cr2"] is None:
                rows[str(L)]={"available":False}
                continue
            rows[str(L)]={
                "available":True,
                "dcr2_dslot":float((bp["cr2"]-bm["cr2"])/(2*eps)),
                "dkmin_dslot":float((bp["kmin"]-bm["kmin"])/(2*eps)),
            }
        per_slot[slot]={"eps":float(eps),"per_L":rows}

    ranking={}
    for L in LS:
        vals=[]
        for slot,r in per_slot.items():
            q=r["per_L"][str(L)]
            if q.get("available"):
                vals.append({
                    "slot":slot,
                    "dcr2_dslot":q["dcr2_dslot"],
                    "abs_dcr2_dslot":abs(q["dcr2_dslot"]),
                    "dkmin_dslot":q["dkmin_dslot"],
                })
        vals.sort(key=lambda z:z["abs_dcr2_dslot"],reverse=True)
        ranking[str(L)]=vals[:20]

    # cross-L robust rank by geometric/mean normalized rank position
    scores={}
    for L,rows in ranking.items():
        for rank,row in enumerate(rows):
            scores.setdefault(row["slot"],[]).append(rank+1)
    consensus=[
        {"slot":s,"mean_top20_rank":float(np.mean(ranks)),"appears_in_L":len(ranks)}
        for s,ranks in scores.items()
    ]
    consensus.sort(key=lambda z:(-z["appears_in_L"],z["mean_top20_rank"]))

    report={
        "scope":"diagnostic 41-slot sensitivity at inner radial zero; not action-consistent repair",
        "center_u_target":U0,
        "baseline":base,
        "per_slot":per_slot,
        "top20_by_L":ranking,
        "consensus":consensus[:20],
        "next_target":"TRACE_TOP_ROBUST_SLOTS_TO_ARCHIVAL_ACTION_PRIMITIVES",
    }
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps({
        "baseline":base,
        "consensus":consensus[:12],
    },indent=2))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
