#!/usr/bin/env python3
"""Targeted inner-tail background-null Hessian continuation v2.

The action-level finite predictor-corrector stalled although first-order
controllability exists.  This solver uses the exact additive background-null
f2(phi,X,F,Y) Hessian freedom H t=0, which leaves the background/value/first
jets unchanged while modifying the perturbation operator through the exact
Appendix-A response.

Unlike the earlier 3-shape global-tail trial, v2 uses localized one-sided
C-infinity bases across the complete inner native tail, including the endpoint
u=0.61.  It directly maximizes the minimum generalized radial characteristic
c_r^2 for all required L with K>0 as a hard guard.

Diagnostic only; no production promotion.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))

from ssz_p5.config import DEFAULT_L  # noqa:E402
from ssz_p5.numerics import module  # noqa:E402
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central  # noqa:E402
from ssz_p5.production.holonomic_hessian_y import TRANSVERSE,complete_hessian,response_from_hessian,chain_residual  # noqa:E402

OUT=ROOT/"data/generated/spectral/INNER_HESSIAN_RADIAL_CONTINUATION_V2.json"
STREAM_OUT=ROOT/"data/generated/spectral/INNER_HESSIAN_RADIAL_CONTINUATION_V2_DIRECT41.csv"

CENTERS=(0.61015,0.6112,0.6125,0.6140,0.6155,0.6170,0.6185,0.61955)
WIDTH=0.00125
LS=tuple(int(x) for x in DEFAULT_L)
EPS=5e-4
TRUST=1.5
MAX_IT=12
INNER_HI=0.62


def bump(u,c):
    z=(np.asarray(u,float)-c)/WIDTH
    q=np.zeros_like(z)
    m=np.abs(z)<1
    q[m]=np.exp(-1/(1-z[m]**2))/np.exp(-1)
    q=np.where(u<=INNER_HI,q,0.0)
    return q


def columns(action):
    u=action.u.to_numpy(float)
    cols=[];meta=[]
    for c in CENTERS:
        sh=bump(u,c)
        for j,name in enumerate(TRANSVERSE):
            q=np.zeros((len(action),6),float)
            q[:,j]=sh
            H,_=complete_hessian(action,q)
            cres,cn=chain_residual(action,H)
            resp=response_from_hessian(action,H)
            cols.append({slot:resp[slot].to_numpy(float) for slot in ("v5","c3","e3","v1","v4","c2")})
            meta.append({
                "center":c,"transverse":name,
                "max_chain_residual":float(np.max(np.abs(cres))),
                "max_normalized_chain_residual":float(np.max(cn)),
            })
    return cols,meta


def add(base,cols,p):
    d=base.copy()
    for slot in ("v5","c3","e3","v1","v4","c2"):
        x=d[slot].to_numpy(float).copy()
        for a,col in zip(p,cols):
            if a:
                x+=float(a)*col[slot]
        d[slot]=x
    return d


def arrays(stream,red,L):
    a=red.canonical_audit(stream,int(L))
    K=np.asarray(a["K"],float);G=np.asarray(a["G"],float)
    Ks=(K+K.transpose(0,2,1))/2
    Gs=(G+G.transpose(0,2,1))/2
    ke=np.linalg.eigvalsh(Ks)[:,0]
    cr=np.full(len(stream),np.nan)
    for i in np.flatnonzero(ke>0):
        w,U=np.linalg.eigh(Ks[i])
        inv=U@np.diag(1/np.sqrt(w))@U.T
        C=inv@Gs[i]@inv
        cr[i]=np.linalg.eigvalsh((C+C.T)/2)[0]
    return ke,cr


def sample_rows(u,n=56):
    ids=np.flatnonzero(u<=INNER_HI)
    pick=np.unique(np.linspace(0,len(ids)-1,min(n,len(ids))).round().astype(int))
    # include closest rows to all centers and both endpoints
    ext=[ids[0],ids[-1]]
    ext += [int(np.argmin(np.abs(u-c))) for c in CENTERS]
    return np.unique(np.concatenate([ids[pick],np.asarray(ext,int)]))


def objective(stream,red,u,rows):
    cr=[];ctags=[];kv=[];ktags=[]
    for L in LS:
        k,c=arrays(stream,red,L)
        for i in rows:
            cr.append(c[i] if np.isfinite(c[i]) else -1e6); ctags.append((L,int(i)))
            kv.append(k[i]); ktags.append((L,int(i)))
    return np.asarray(cr),ctags,np.asarray(kv),ktags


def worst(vals,tags,u,key):
    j=int(np.argmin(vals));L,i=tags[j]
    return {key:float(vals[j]),"L":int(L),"u":float(u[i]),"row":int(i)}


def main():
    b=build_onshell_central(ROOT)
    action=b.action.sort_values("x").reset_index(drop=True)
    base=b.direct41.sort_values("x").reset_index(drop=True)
    u=base.u.to_numpy(float)
    red=module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    cols,meta=columns(action)
    rows=sample_rows(u)
    p=np.zeros(len(cols))
    cur=base.copy()
    cr,ct,kv,kt=objective(cur,red,u,rows)
    init_cr=worst(cr,ct,u,"min_cr2"); init_k=worst(kv,kt,u,"min_K")
    hist=[]

    for it in range(MAX_IT):
        Jc=np.zeros((len(cr),len(cols)))
        Jk=np.zeros((len(kv),len(cols)))
        for j in range(len(cols)):
            pp=p.copy();pp[j]+=EPS
            cc,cct,kk,kkt=objective(add(base,cols,pp),red,u,rows)
            Jc[:,j]=(cc-cr)/EPS
            Jk[:,j]=(kk-kv)/EPS

        obj=np.zeros(len(cols)+1);obj[-1]=-1
        Aub=[];bub=[]
        for q in range(len(cr)):
            row=np.zeros(len(cols)+1);row[:-1]=-Jc[q];row[-1]=1
            Aub.append(row);bub.append(float(cr[q]))
        for q in range(len(kv)):
            row=np.zeros(len(cols)+1);row[:-1]=-Jk[q]
            Aub.append(row);bub.append(float(kv[q]-1e-8))
        sol=linprog(obj,A_ub=np.asarray(Aub),b_ub=np.asarray(bub),
                    bounds=[(-TRUST,TRUST)]*len(cols)+[(None,None)],method="highs")
        if not sol.success:
            hist.append({"iteration":it,"status":"LP_FAIL","message":sol.message});break

        dp=sol.x[:-1]
        old=float(np.min(cr))
        best=None;accepted=False
        for frac in (1,.5,.25,.125,.0625,.03125,.015625):
            pt=p+frac*dp
            st=add(base,cols,pt)
            cc,cct,kk,kkt=objective(st,red,u,rows)
            wi=worst(cc,cct,u,"min_cr2");wk=worst(kk,kkt,u,"min_K")
            rec={"iteration":it,"fraction":frac,"predicted_margin":float(sol.x[-1]),
                 "worst_cr2":wi,"worst_K":wk,"step_max_abs":float(np.max(np.abs(frac*dp)))}
            score=wi["min_cr2"]
            if best is None or score>best[0]:
                best=(score,pt,st,cc,cct,kk,kkt,rec)
            if score>old+max(1e-6,1e-4*max(1,abs(old))) and wk["min_K"]>1e-8:
                p,cur,cr,ct,kv,kt=pt,st,cc,cct,kk,kkt
                rec["status"]="ACCEPT";hist.append(rec);accepted=True;break
        if not accepted:
            if best and best[0]>old and best[-1]["worst_K"]["min_K"]>1e-8:
                _,p,cur,cr,ct,kv,kt,rec=best
                rec["status"]="ACCEPT_BEST";hist.append(rec);accepted=True
            else:
                if best:
                    best[-1]["status"]="STALL";hist.append(best[-1])
                break

        # exact full inner verify for early stop
        inner=np.flatnonzero(u<=INNER_HI)
        fcr=[];fct=[];fk=[];fkt=[]
        for L in LS:
            k,c=arrays(cur,red,L)
            for i in inner:
                fcr.append(c[i] if np.isfinite(c[i]) else -1e6);fct.append((L,int(i)))
                fk.append(k[i]);fkt.append((L,int(i)))
        if min(fcr)>0 and min(fk)>0:
            break

    inner=np.flatnonzero(u<=INNER_HI)
    prod=np.flatnonzero((u>0.62)&(u<0.70))
    full={}
    all_cr=True;all_k=True;prod_ok=True
    for L in LS:
        k,c=arrays(cur,red,L)
        rec={
            "inner_min_K":float(np.min(k[inner])),
            "inner_min_cr2":float(np.nanmin(c[inner])),
            "inner_negative_cr2_rows":int(np.sum(c[inner]<0)),
            "production_min_K":float(np.min(k[prod])),
        }
        full[str(L)]=rec
        all_k &= rec["inner_min_K"]>0
        all_cr &= rec["inner_min_cr2"]>0
        prod_ok &= rec["production_min_K"]>0

    final_cr=min((v["inner_min_cr2"],L) for L,v in full.items())
    if all_cr and all_k and prod_ok:
        status="INNER_HESSIAN_RADIAL_CONTINUATION_PASS"
        nxt="COMBINE_WITH_OUTER_REPAIR_AND_RUN_FULL_NATIVE_PRECERTIFICATE"
    elif final_cr[0]>init_cr["min_cr2"]:
        status="INNER_HESSIAN_RADIAL_CONTINUATION_IMPROVED_NOT_CLOSED"
        nxt="ADD_MORE_LOCAL_HESSIAN_KNOTS_OR_COUPLED_ACTION_CONTROLS"
    else:
        status="INNER_HESSIAN_RADIAL_CONTINUATION_STALLED"
        nxt="EXPAND_BEYOND_F2_BACKGROUND_NULL_HESSIAN"

    STREAM_OUT.parent.mkdir(parents=True,exist_ok=True)
    cur.to_csv(STREAM_OUT,index=False)
    report={
        "status":status,"next_target":nxt,
        "scope":"localized background-null f2 Hessian continuation on full inner native tail",
        "initial_worst_cr2":init_cr,"initial_worst_K":init_k,
        "centers":list(CENTERS),"width":WIDTH,
        "control_meta":meta,
        "history":hist,
        "parameters":[float(x) for x in p],
        "parameter_max_abs":float(np.max(np.abs(p))),
        "full_verification":full,
        "all_inner_K_pass":bool(all_k),
        "all_inner_radial_pass":bool(all_cr),
        "all_production_K_pass":bool(prod_ok),
        "stream_csv":str(STREAM_OUT.relative_to(ROOT)),
    }
    OUT.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps({
        "status":status,
        "initial_cr2":init_cr,
        "final_min_cr2":final_cr,
        "all_inner_K":all_k,
        "all_inner_radial":all_cr,
        "production":prod_ok,
        "next":nxt,
    },indent=2,allow_nan=False))
    return 0

if __name__=="__main__":
    raise SystemExit(main())
