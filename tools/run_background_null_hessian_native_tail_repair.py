#!/usr/bin/env python3
"""Endpoint-aware action-consistent background-null Hessian repair.

The previous two-envelope trial accidentally forced the control envelope to
vanish at the exact native endpoints u=0.61 and u=u_max, precisely where the
worst inner c_r^2 and outer K values occur.  This corrected trial uses
one-sided smooth tail bases which vanish with zero derivatives at the
production seams (u=0.62 and u=0.70) but remain active at the native endpoints.

Every deformation is an exact additive background-null f2(phi,X,F,Y) Hessian
satisfying H@(phi',X',F',Y')=0 via holonomic_hessian_y.  Only the exact
Appendix-A response (v5,c3,e3,v1,v4,c2) is added to the current 41-stream.
The background and first action jets are therefore untouched.

Outer K and inner c_r^2 are optimized independently because the controls have
disjoint support.  The final candidate is then verified on every native row and
all required L.  Nothing is promoted automatically.
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

OUT=ROOT/"data/generated/spectral/BACKGROUND_NULL_HESSIAN_NATIVE_TAIL_REPAIR.json"
STREAM_OUT=ROOT/"data/generated/spectral/BACKGROUND_NULL_HESSIAN_NATIVE_TAIL_REPAIR_DIRECT41.csv"

LS_OPT=(6,42,110,1000)
LS_VERIFY=tuple(int(x) for x in DEFAULT_L)
EPS=2e-4
TRUST=0.8
MAX_IT=8
SHAPES=3


def smooth01(t):
    t=np.clip(np.asarray(t,float),0.0,1.0)
    return t**3*(10.0-15.0*t+6.0*t*t)


def tail_basis(u,region):
    u=np.asarray(u,float)
    if region=="inner":
        t=np.clip((0.62-u)/0.01,0.0,1.0)
        support=u<=0.62
    elif region=="outer":
        umax=float(np.max(u))
        t=np.clip((u-0.70)/max(umax-0.70,1e-15),0.0,1.0)
        support=u>=0.70
    else:
        raise ValueError(region)
    s=smooth01(t)
    B=[
        s,
        s*(4.0*t*(1.0-t)),
        s*(16.0*t*t*(1.0-t)*(1.0-t)),
    ]
    return [np.where(support,b,0.0) for b in B]


def response_columns(action,region):
    basis=tail_basis(action.u.to_numpy(float),region)
    cols=[]; meta=[]
    for ib,b in enumerate(basis):
        for j,name in enumerate(TRANSVERSE):
            q=np.zeros((len(action),6),float)
            q[:,j]=b
            H,_=complete_hessian(action,q)
            cres,cnorm=chain_residual(action,H)
            resp=response_from_hessian(action,H)
            cols.append({slot:resp[slot].to_numpy(float) for slot in ("v5","c3","e3","v1","v4","c2")})
            meta.append({
                "region":region,"shape":ib,"transverse":name,
                "endpoint_value":float(b[0] if region=="outer" else b[-1]),
                "max_abs_chain_residual":float(np.max(np.abs(cres))),
                "max_normalized_chain_residual":float(np.max(cnorm)),
            })
    return cols,meta


def add(base,cols,p):
    d=base.copy()
    for slot in ("v5","c3","e3","v1","v4","c2"):
        x=d[slot].to_numpy(float).copy()
        for a,c in zip(p,cols):
            if a: x+=float(a)*c[slot]
        d[slot]=x
    return d


def audit_arrays(stream,red,L):
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


def sampled_rows(u,region,n=17):
    ids=np.flatnonzero(u<=0.62) if region=="inner" else np.flatnonzero(u>=0.70)
    pick=np.unique(np.linspace(0,len(ids)-1,min(n,len(ids))).round().astype(int))
    # always endpoints and seam-nearest rows
    return ids[pick]


def objective(stream,red,u,region):
    ids=sampled_rows(u,region)
    vals=[];tags=[]
    for L in LS_OPT:
        k,c=audit_arrays(stream,red,L)
        if region=="outer":
            for i in ids:
                vals.append(k[i]);tags.append(("K",L,int(i)))
        else:
            for i in ids:
                vals.append(c[i]);tags.append(("cr2",L,int(i)))
                vals.append(k[i]);tags.append(("Kguard",L,int(i)))
    return np.asarray(vals,float),tags


def worst(vals,tags,u):
    out={}
    for kind in sorted(set(t[0] for t in tags)):
        ids=[i for i,t in enumerate(tags) if t[0]==kind]
        q=ids[int(np.nanargmin(vals[ids]))]
        _,L,row=tags[q]
        out[kind]={"value":float(vals[q]),"L":int(L),"u":float(u[row]),"row":int(row)}
    return out


def optimize_region(base,cols,red,u,region):
    p=np.zeros(len(cols))
    st=add(base,cols,p)
    vals,tags=objective(st,red,u,region)
    hist=[]
    target_kind="K" if region=="outer" else "cr2"
    for it in range(MAX_IT):
        finite=np.where(np.isfinite(vals),vals,-1e6)
        J=np.zeros((len(vals),len(cols)))
        for j in range(len(cols)):
            pp=p.copy();pp[j]+=EPS
            vv,_=objective(add(base,cols,pp),red,u,region)
            J[:,j]=(np.where(np.isfinite(vv),vv,-1e6)-finite)/EPS

        c=np.zeros(len(cols)+1);c[-1]=-1
        Aub=[];bub=[]
        for q,tag in enumerate(tags):
            row=np.zeros(len(cols)+1)
            if tag[0]==target_kind:
                row[:-1]=-J[q];row[-1]=1.0
                Aub.append(row);bub.append(finite[q])
            elif tag[0]=="Kguard":
                row[:-1]=-J[q]
                Aub.append(row);bub.append(finite[q]-1e-8)
        sol=linprog(c,A_ub=np.asarray(Aub),b_ub=np.asarray(bub),
                    bounds=[(-TRUST,TRUST)]*len(cols)+[(None,None)],method="highs")
        if not sol.success:
            hist.append({"iteration":it,"status":"LP_FAIL","message":sol.message});break
        dp=sol.x[:-1]
        old=worst(vals,tags,u)
        oldscore=old[target_kind]["value"]
        best=None;accepted=False
        for frac in (1.0,.5,.25,.125,.0625,.03125):
            pt=p+frac*dp
            cand=add(base,cols,pt)
            vv,tt=objective(cand,red,u,region)
            if not np.all(np.isfinite(vv)): continue
            ww=worst(vv,tt,u)
            if region=="inner" and ww["Kguard"]["value"]<=1e-8: continue
            score=ww[target_kind]["value"]
            rec={"iteration":it,"fraction":frac,"target":score,"worst":ww,
                 "step_max_abs":float(np.max(np.abs(frac*dp)))}
            if best is None or score>best[0]:best=(score,pt,cand,vv,tt,rec)
            if score>oldscore+max(1e-7,1e-5*max(1.0,abs(oldscore))):
                p,st,vals,tags=pt,cand,vv,tt
                rec["status"]="ACCEPT";hist.append(rec);accepted=True;break
        if not accepted:
            if best and best[0]>oldscore:
                _,p,st,vals,tags,rec=best
                rec["status"]="ACCEPT_BEST";hist.append(rec);accepted=True
            else:
                hist.append({"iteration":it,"status":"STALL","worst":old});break
        if worst(vals,tags,u)[target_kind]["value"]>0:
            break
    return p,st,{"initial":worst(objective(base,red,u,region)[0],objective(base,red,u,region)[1],u),
                 "final":worst(vals,tags,u),"history":hist}


def verify(stream,red,u):
    inner=np.flatnonzero(u<=0.62);outer=np.flatnonzero(u>=0.70);prod=np.flatnonzero((u>0.62)&(u<0.70))
    rep={}
    for L in LS_VERIFY:
        k,c=audit_arrays(stream,red,L)
        rep[str(L)]={
            "production_min_K":float(np.min(k[prod])),
            "outer_min_K":float(np.min(k[outer])),
            "outer_negative_K_rows":int(np.sum(k[outer]<=0)),
            "outer_min_cr2_where_K_positive":float(np.nanmin(c[outer])) if np.any(np.isfinite(c[outer])) else None,
            "outer_negative_cr2_rows_where_K_positive":int(np.sum(c[outer]<0)),
            "inner_min_K":float(np.min(k[inner])),
            "inner_min_cr2":float(np.nanmin(c[inner])),
            "inner_negative_cr2_rows":int(np.sum(c[inner]<0)),
        }
    return rep


def main():
    b=build_onshell_central(ROOT)
    action=b.action.sort_values("x").reset_index(drop=True)
    base=b.direct41.sort_values("x").reset_index(drop=True)
    u=base.u.to_numpy(float)
    red=module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")

    outer_cols,outer_meta=response_columns(action,"outer")
    inner_cols,inner_meta=response_columns(action,"inner")

    po,outer_stream,outer_opt=optimize_region(base,outer_cols,red,u,"outer")
    pi,_,inner_opt=optimize_region(base,inner_cols,red,u,"inner")

    # combine disjoint-support exact responses
    final=add(add(base,outer_cols,po),inner_cols,pi)
    full=verify(final,red,u)

    outerK=all(v["outer_min_K"]>0 for v in full.values())
    outerG=all(v["outer_min_cr2_where_K_positive"] is not None and v["outer_min_cr2_where_K_positive"]>0 for v in full.values())
    innerK=all(v["inner_min_K"]>0 for v in full.values())
    innerG=all(v["inner_min_cr2"]>0 for v in full.values())
    prod=all(v["production_min_K"]>0 for v in full.values())

    STREAM_OUT.parent.mkdir(parents=True,exist_ok=True)
    final.to_csv(STREAM_OUT,index=False)

    if outerK and outerG and innerK and innerG and prod:
        status="BACKGROUND_NULL_HESSIAN_NATIVE_TAIL_REPAIR_PASS"
        nxt="RUN_NATIVE_PRECERTIFICATE_AND_GLOBAL_KRGM_CONSTRUCTION"
    elif (
        outer_opt["final"]["K"]["value"]>outer_opt["initial"]["K"]["value"]
        or inner_opt["final"]["cr2"]["value"]>inner_opt["initial"]["cr2"]["value"]
    ):
        status="BACKGROUND_NULL_HESSIAN_NATIVE_TAIL_REPAIR_IMPROVED_NOT_CLOSED"
        nxt="EXPAND_ENDPOINT_AWARE_HESSIAN_BASIS_AND_ADD_RADIAL_G_GUARDS"
    else:
        status="BACKGROUND_NULL_HESSIAN_NATIVE_TAIL_REPAIR_STALLED"
        nxt="USE_ADDITIONAL_BACKGROUND_NULL_ACTION_SECTORS"

    report={
        "status":status,"next_target":nxt,
        "correction_note":"supersedes prior envelope trial which vanished at both pathological native endpoints",
        "outer_control_meta":outer_meta,"inner_control_meta":inner_meta,
        "outer_optimization":outer_opt,"inner_optimization":inner_opt,
        "outer_parameters":[float(x) for x in po],"inner_parameters":[float(x) for x in pi],
        "full_verification":full,
        "all_outer_K_pass":outerK,"all_outer_radial_pass":outerG,
        "all_inner_K_pass":innerK,"all_inner_radial_pass":innerG,
        "all_production_K_pass":prod,
        "stream_csv":str(STREAM_OUT.relative_to(ROOT)),
    }
    OUT.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps(report,indent=2,allow_nan=False))
    return 0

if __name__=="__main__":raise SystemExit(main())
