#!/usr/bin/env python3
"""Finite Eq131-compatible inner radial predictor-corrector.

The inner-tail controllability audit found a common first-order action tangent
that increases the weakest radial characteristic c_r^2 for all required L
while preserving E00,E11,JA,Eq131 and not decreasing K.

This script attempts the corresponding finite continuation on every native
inner-tail row u<=0.62.  It modifies only u<0.61995 so the registered production
window is untouched, re-emits the full action after every step, and verifies
both K>0 and c_r^2>0 on the complete inner tail.

Diagnostic only: no production overwrite or promotion.
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
from ssz_p5.production.full_action_lower import complete_total_action_jets, emit_lower_slots  # noqa:E402
from ssz_p5.production.central_action import background_residuals  # noqa:E402

OUT = ROOT / "data/generated/spectral/CURRENT_MEMBER_INNER_RADIAL_FINITE_CONTINUATION.json"
ACTION_OUT = ROOT / "data/generated/spectral/CURRENT_MEMBER_INNER_RADIAL_FINITE_CONTINUATION_ACTION.csv"
STREAM_OUT = ROOT / "data/generated/spectral/CURRENT_MEMBER_INNER_RADIAL_FINITE_CONTINUATION_DIRECT41.csv"

DIRECTIONS = ("f2","f2X","f2F","f3","f3X","f4","f4X","f4XX")
CENTERS = (0.6190,0.6170,0.6150,0.6130,0.6112,0.61035)
WIDTH = 7.0e-4
LS = tuple(int(x) for x in DEFAULT_L)
FD_EPS = 1e-5
MAX_IT = 8
TRUST = 0.04
INNER_HI = 0.61995


def bump(u,u0):
    z=(np.asarray(u,float)-u0)/WIDTH
    q=np.zeros_like(z)
    m=np.abs(z)<1
    q[m]=np.exp(-1/(1-z[m]**2))/np.exp(-1)
    return q


def reemit(action):
    zk=module("ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py")
    d=action.copy().sort_values("x").reset_index(drop=True)
    r=d.x.to_numpy(float)
    f=d.f.to_numpy(float); h=d.h.to_numpy(float)
    X=d.X.to_numpy(float); A=d.A0prime.to_numpy(float); ph=d.phiprime.to_numpy(float)
    F=h*A*A/(2*f); Y=4*X*F
    Xp=zk.dr(r,X,1,9,8); Fp=zk.dr(r,F,1,9,8); Yp=zk.dr(r,Y,1,9,8)
    d["f2phi"]=(zk.dr(r,d.f2.to_numpy(float),1,9,8)
                 -d.f2X.to_numpy(float)*Xp
                 -d.f2F.to_numpy(float)*Fp
                 -d.f2Y.to_numpy(float)*Yp)/ph
    d=complete_total_action_jets(d)
    lower,d=emit_lower_slots(d)
    out=zk.emit(d,
        selected_v5=lower.v5.to_numpy(float),
        selected_c3=lower.c3.to_numpy(float),
        selected_e3=lower.e3.to_numpy(float),
        v6_phi_selector="action")
    out["a5"]=(zk.dr(r,out.a2.to_numpy(float),1,9,8)
               -zk.dr(r,out.a1.to_numpy(float),2,9,8)
               -zk.dr(r,out.A0prime.to_numpy(float)*out.v4.to_numpy(float)/2,1,9,8)
               +out.A0prime.to_numpy(float)*out.v5.to_numpy(float)/2)
    return d,out


def apply_params(base,p):
    d=base.copy()
    u=d.u.to_numpy(float)
    pp=np.asarray(p,float).reshape(len(CENTERS),len(DIRECTIONS))
    for ic,c in enumerate(CENTERS):
        sh=bump(u,c)
        # hard guard: exactly zero in/above production-side interface
        sh=np.where(u<INNER_HI,sh,0.0)
        for j,name in enumerate(DIRECTIONS):
            d[name]=d[name].to_numpy(float)+pp[ic,j]*sh
    return d


def constraints(comp,colloc):
    bg=background_residuals(comp)
    out=[]
    for i in colloc:
        out.extend((bg["E00"][i],bg["E11"][i],bg["JA"][i],bg["eq131_residual"][i]))
    return np.asarray(out,float)


def metrics(stream,red,rows):
    crvals=[]; crtags=[]; kvals=[]; ktags=[]
    for L in LS:
        a=red.canonical_audit(stream,int(L))
        K=np.asarray(a["K"],float); G=np.asarray(a["G"],float)
        for i in rows:
            Ks=(K[i]+K[i].T)/2; Gs=(G[i]+G[i].T)/2
            w,U=np.linalg.eigh(Ks)
            km=float(w[0]); kvals.append(km); ktags.append((L,i))
            if km<=0:
                cr=-1e6
            else:
                inv=U@np.diag(1/np.sqrt(w))@U.T
                C=inv@Gs@inv
                cr=float(np.linalg.eigvalsh((C+C.T)/2)[0])
            crvals.append(cr); crtags.append((L,i))
    return np.asarray(crvals),crtags,np.asarray(kvals),ktags


def worst(vals,tags,u,key):
    j=int(np.argmin(vals)); L,i=tags[j]
    return {key:float(vals[j]),"L":int(L),"u":float(u[i]),"row":int(i)}


def main():
    red=module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    build=build_onshell_central(ROOT)
    base=build.action.sort_values("x").reset_index(drop=True)
    u=base.u.to_numpy(float)
    inner=np.flatnonzero(u<=0.62)
    if len(inner)<50:
        raise RuntimeError("insufficient inner rows")
    colloc=np.array([int(np.argmin(np.abs(u-c))) for c in CENTERS],dtype=int)
    pick=np.linspace(0,len(inner)-1,min(44,len(inner))).round().astype(int)
    sample=np.unique(np.concatenate([inner[pick],colloc,[inner[0],inner[-1]]]))

    npar=len(CENTERS)*len(DIRECTIONS)
    p=np.zeros(npar)
    hist=[]

    comp,stream=reemit(apply_params(base,p))
    E=constraints(comp,colloc)
    cr,crtags,kv,ktags=metrics(stream,red,sample)
    initial_cr=worst(cr,crtags,u,"min_cr2")
    initial_k=worst(kv,ktags,u,"min_K")

    for it in range(MAX_IT):
        B=np.zeros((len(E),npar))
        Gcr=np.zeros((len(cr),npar))
        Gk=np.zeros((len(kv),npar))
        for j in range(npar):
            pp=p.copy(); pp[j]+=FD_EPS
            cp,sp=reemit(apply_params(base,pp))
            Ep=constraints(cp,colloc)
            crp,_,kp,_=metrics(sp,red,sample)
            B[:,j]=(Ep-E)/FD_EPS
            Gcr[:,j]=(crp-cr)/FD_EPS
            Gk[:,j]=(kp-kv)/FD_EPS

        rn=np.linalg.norm(B,axis=1)
        keep=rn>1e-12
        Be=B[keep]/rn[keep,None]
        ee=E[keep]/rn[keep]

        # maximize common radial margin t, while requiring linearized K>0
        obj=np.zeros(npar+1); obj[-1]=-1
        Aub=[]; bub=[]
        for q in range(len(cr)):
            row=np.zeros(npar+1); row[:npar]=-Gcr[q]; row[-1]=1
            Aub.append(row); bub.append(float(cr[q]))
        kfloor=1e-8
        for q in range(len(kv)):
            row=np.zeros(npar+1); row[:npar]=-Gk[q]
            Aub.append(row); bub.append(float(kv[q]-kfloor))
        Aeq=np.column_stack([Be,np.zeros(len(Be))]) if len(Be) else None
        beq=-ee if len(Be) else None
        sol=linprog(obj,A_ub=np.asarray(Aub),b_ub=np.asarray(bub),
                    A_eq=Aeq,b_eq=beq,
                    bounds=[(-TRUST,TRUST)]*npar+[(None,None)],method="highs")
        if not sol.success:
            hist.append({"iteration":it,"status":"LP_FAIL","message":sol.message})
            break

        dp=sol.x[:-1]
        cur=float(np.min(cr))
        accepted=False
        best=None
        for frac in (1.0,.5,.25,.125,.0625,.03125):
            pt=p+frac*dp
            ct,st=reemit(apply_params(base,pt))
            Et=constraints(ct,colloc)
            crt,crtt,kvt,kvtt=metrics(st,red,sample)
            em=float(np.max(np.abs(Et)))
            wi=worst(crt,crtt,u,"min_cr2")
            wk=worst(kvt,kvtt,u,"min_K")
            rec={"iteration":it,"fraction":frac,"predicted_margin":float(sol.x[-1]),
                 "finite_worst_cr2":wi,"finite_worst_K":wk,
                 "constraint_max_abs":em,"step_max_abs":float(np.max(np.abs(frac*dp)))}
            score=wi["min_cr2"]-10*em
            if best is None or score>best[0]:
                best=(score,pt,ct,st,Et,crt,crtt,kvt,kvtt,rec)
            if wi["min_cr2"]>cur+1e-5 and wk["min_K"]>0 and em<2e-3:
                p,comp,stream,E,cr,crtags,kv,ktags=pt,ct,st,Et,crt,crtt,kvt,kvtt
                rec["status"]="ACCEPT"; hist.append(rec); accepted=True; break
        if not accepted:
            if best and best[-1]["finite_worst_K"]["min_K"]>0 and best[-1]["constraint_max_abs"]<5e-3 and best[-1]["finite_worst_cr2"]["min_cr2"]>cur:
                _,p,comp,stream,E,cr,crtags,kv,ktags,rec=best
                rec["status"]="ACCEPT_BEST"; hist.append(rec); accepted=True
            else:
                if best:
                    best[-1]["status"]="STALL"; hist.append(best[-1])
                break

        fullcr,fullcrtags,fullk,fullktags=metrics(stream,red,inner)
        if np.min(fullcr)>0 and np.min(fullk)>0 and np.max(np.abs(E))<2e-3:
            break

    final_action=apply_params(base,p)
    comp,stream=reemit(final_action)
    fullcr,fullcrtags,fullk,fullktags=metrics(stream,red,inner)
    final_cr=worst(fullcr,fullcrtags,u,"min_cr2")
    final_k=worst(fullk,fullktags,u,"min_K")
    bg=background_residuals(comp)
    imask=u<=0.62
    bgmax={k:float(np.max(np.abs(np.asarray(v)[imask]))) for k,v in bg.items()}

    prod=(u>0.62)&(u<0.70)
    prod_shift=0.0
    for name in DIRECTIONS:
        prod_shift=max(prod_shift,float(np.max(np.abs(
            final_action[name].to_numpy(float)[prod]-base[name].to_numpy(float)[prod]))))

    if final_cr["min_cr2"]>0 and final_k["min_K"]>0 and all(bgmax[k]<2e-3 for k in ("E00","E11","JA","eq131_residual")):
        status="INNER_RADIAL_FINITE_CONTINUATION_PASS"
        nxt="COMBINE_WITH_OUTER_K_REPAIR_AND_REAUDIT_FULL_NATIVE_MEMBER"
    elif final_cr["min_cr2"]>initial_cr["min_cr2"]:
        status="INNER_RADIAL_FINITE_CONTINUATION_IMPROVED_NOT_CLOSED"
        nxt="REFINE_INNER_BASIS_AND_NONLINEAR_CORRECTOR"
    else:
        status="INNER_RADIAL_FINITE_CONTINUATION_STALLED"
        nxt="EXPAND_INNER_ACTION_BASIS_OR_ROOT_CORRECTOR"

    OUT.parent.mkdir(parents=True,exist_ok=True)
    comp.to_csv(ACTION_OUT,index=False)
    stream.to_csv(STREAM_OUT,index=False)
    report={
        "scope":"finite Eq131-compatible inner radial continuation; no promotion",
        "inner_domain":[float(u[inner].min()),float(u[inner].max())],
        "centers":list(CENTERS),"directions":list(DIRECTIONS),"width":WIDTH,
        "initial_worst_cr2":initial_cr,"initial_worst_K":initial_k,
        "iterations":hist,
        "final_worst_cr2":final_cr,"final_worst_K":final_k,
        "background_max_abs_inner":bgmax,
        "production_action_max_abs_shift":prod_shift,
        "parameter_max_abs":float(np.max(np.abs(p))),
        "nonzero_parameter_count":int(np.sum(np.abs(p)>1e-10)),
        "status":status,"next_target":nxt,
        "artifacts":{
            "action_csv":str(ACTION_OUT.relative_to(ROOT)),
            "direct41_csv":str(STREAM_OUT.relative_to(ROOT)),
        }
    }
    OUT.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps(report,indent=2,allow_nan=False))
    return 0

if __name__=="__main__":
    raise SystemExit(main())
