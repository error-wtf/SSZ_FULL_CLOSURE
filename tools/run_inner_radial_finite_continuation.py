#!/usr/bin/env python3
"""Finite predictor-corrector trial for the native inner radial c_r^2 failure.

Licensed by INNER_RADIAL_EQ131_COMPATIBLE_COMMON_TANGENT_EXISTS.

The trial modifies only the native inner tail (u<0.6198), leaving the registered
production window untouched.  It uses smooth compact action-level bumps in
f2,f2X,f2F,f3,f3X,f4,f4X,f4XX, re-completes all jets, re-emits all 41 slots,
and iteratively maximizes the worst c_r,min^2 over all native inner rows and
required L while correcting E00,E11,JA,Eq131 at collocation points.

No result is promoted or frozen as production member.
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
from ssz_p5.production.full_action_lower import complete_total_action_jets,emit_lower_slots  # noqa:E402
from ssz_p5.production.central_action import background_residuals  # noqa:E402

OUT=ROOT/"data/generated/spectral/INNER_RADIAL_FINITE_CONTINUATION_TRIAL.json"

DIRS=("f2","f2X","f2F","f3","f3X","f4","f4X","f4XX")
CENTERS=(0.6115,0.6132,0.6150,0.6167,0.6184)
WIDTH=0.00115
INNER_HI=0.6198
LS=tuple(DEFAULT_L)
FD_EPS=1e-5
TRUST=0.04
MAX_IT=8


def bump(u,u0):
    z=(np.asarray(u,float)-u0)/WIDTH
    q=np.zeros_like(z);m=np.abs(z)<1
    q[m]=np.exp(-1/(1-z[m]**2))/np.exp(-1)
    return q


def reemit(action):
    zk=module("ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py")
    d=action.copy().sort_values("x").reset_index(drop=True)
    r=d.x.to_numpy(float);f=d.f.to_numpy(float);h=d.h.to_numpy(float)
    X=d.X.to_numpy(float);A=d.A0prime.to_numpy(float);ph=d.phiprime.to_numpy(float)
    F=h*A*A/(2*f);Y=4*X*F
    Xp=zk.dr(r,X,1,9,8);Fp=zk.dr(r,F,1,9,8);Yp=zk.dr(r,Y,1,9,8)
    d["f2phi"]=(zk.dr(r,d.f2.to_numpy(float),1,9,8)-d.f2X.to_numpy(float)*Xp
                -d.f2F.to_numpy(float)*Fp-d.f2Y.to_numpy(float)*Yp)/ph
    d=complete_total_action_jets(d)
    lower,d=emit_lower_slots(d)
    out=zk.emit(d,selected_v5=lower.v5.to_numpy(float),selected_c3=lower.c3.to_numpy(float),
                selected_e3=lower.e3.to_numpy(float),v6_phi_selector="action")
    out["a5"]=(zk.dr(r,out.a2.to_numpy(float),1,9,8)-zk.dr(r,out.a1.to_numpy(float),2,9,8)
               -zk.dr(r,out.A0prime.to_numpy(float)*out.v4.to_numpy(float)/2,1,9,8)
               +out.A0prime.to_numpy(float)*out.v5.to_numpy(float)/2)
    return d,out


def apply(base,p):
    d=base.copy();u=d.u.to_numpy(float)
    pp=np.asarray(p,float).reshape(len(CENTERS),len(DIRS))
    for ic,c in enumerate(CENTERS):
        sh=bump(u,c)
        # hard guard: exactly zero above inner boundary
        sh=np.where(u<INNER_HI,sh,0.0)
        for j,n in enumerate(DIRS):
            d[n]=d[n].to_numpy(float)+pp[ic,j]*sh
    return d


def constraints(comp,colloc):
    bg=background_residuals(comp);v=[]
    for i in colloc:
        v.extend((bg["E00"][i],bg["E11"][i],bg["JA"][i],bg["eq131_residual"][i]))
    return np.asarray(v,float)


def cr2_vector(stream,red,rows):
    vals=[];tags=[];kmins=[]
    for L in LS:
        a=red.canonical_audit(stream,int(L))
        K=np.asarray(a["K"],float);G=np.asarray(a["G"],float)
        for i in rows:
            Ks=(K[i]+K[i].T)/2;Gs=(G[i]+G[i].T)/2
            w,U=np.linalg.eigh(Ks);kmins.append(float(np.min(w)))
            if np.min(w)<=0:
                vals.append(-1e6);tags.append((L,i));continue
            inv=U@np.diag(1/np.sqrt(w))@U.T
            C=inv@Gs@inv;C=(C+C.T)/2
            vals.append(float(np.linalg.eigvalsh(C)[0]));tags.append((L,i))
    return np.asarray(vals,float),tags,float(np.min(kmins))


def worst(v,tags,u):
    j=int(np.argmin(v));L,i=tags[j]
    return {"min_cr2":float(v[j]),"L":int(L),"u":float(u[i]),"row":int(i)}


def main():
    red=module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    b=build_onshell_central(ROOT)
    base=b.action.sort_values("x").reset_index(drop=True)
    u=base.u.to_numpy(float)
    inner=np.flatnonzero(u<INNER_HI)
    colloc=np.array([int(np.argmin(np.abs(u-c))) for c in CENTERS],int)
    pick=np.linspace(0,len(inner)-1,42).round().astype(int)
    sample=np.unique(np.concatenate([inner[pick],colloc,[inner[0],inner[-1]]]))

    npar=len(CENTERS)*len(DIRS);p=np.zeros(npar);hist=[]
    comp,st=reemit(apply(base,p));E=constraints(comp,colloc)
    vals,tags,kmin=cr2_vector(st,red,sample)
    init=worst(vals,tags,u)
    full0,t0,k0=cr2_vector(st,red,inner);init_full=worst(full0,t0,u)

    for it in range(MAX_IT):
        B=np.zeros((len(E),npar));GJ=np.zeros((len(vals),npar))
        for j in range(npar):
            pp=p.copy();pp[j]+=FD_EPS
            cp,sp=reemit(apply(base,pp));Ep=constraints(cp,colloc)
            vp,_,_=cr2_vector(sp,red,sample)
            B[:,j]=(Ep-E)/FD_EPS;GJ[:,j]=(vp-vals)/FD_EPS
        rn=np.linalg.norm(B,axis=1);keep=rn>1e-12
        Be=B[keep]/rn[keep,None];ee=E[keep]/rn[keep]
        obj=np.zeros(npar+1);obj[-1]=-1
        Aub=np.column_stack([-GJ,np.ones(len(vals))]);bub=vals.copy()
        Aeq=np.column_stack([Be,np.zeros(len(Be))]) if len(Be) else None
        sol=linprog(obj,A_ub=Aub,b_ub=bub,A_eq=Aeq,b_eq=-ee if len(Be) else None,
                    bounds=[(-TRUST,TRUST)]*npar+[(None,None)],method="highs")
        if not sol.success:
            hist.append({"iteration":it,"status":"LP_FAIL","message":sol.message});break
        dp=sol.x[:-1];cur=float(vals.min());accepted=False;best=None
        for frac in (1,.5,.25,.125,.0625,.03125):
            pt=p+frac*dp;ct,stt=reemit(apply(base,pt));Et=constraints(ct,colloc)
            vt,tt,km=cr2_vector(stt,red,sample);wi=worst(vt,tt,u);em=float(np.max(np.abs(Et)))
            rec={"iteration":it,"fraction":frac,"predicted_margin":float(sol.x[-1]),
                 "finite_sample_worst":wi,"constraint_max_abs":em,"min_K_sample":km}
            score=wi["min_cr2"]-10*em
            if best is None or score>best[0]:best=(score,pt,ct,stt,Et,vt,tt,rec)
            if wi["min_cr2"]>cur+1e-5 and em<2e-3 and km>0:
                p,comp,st,E,vals,tags=pt,ct,stt,Et,vt,tt
                rec["status"]="ACCEPT";hist.append(rec);accepted=True;break
        if not accepted:
            if best and best[7]["finite_sample_worst"]["min_cr2"]>cur and best[7]["constraint_max_abs"]<5e-3 and best[7]["min_K_sample"]>0:
                _,p,comp,st,E,vals,tags,rec=best;rec["status"]="ACCEPT_BEST";hist.append(rec);accepted=True
            else:
                if best: best[7]["status"]="STALL";hist.append(best[7])
                break
        vf,tf,kf=cr2_vector(st,red,inner)
        if vf.min()>0 and kf>0 and np.max(np.abs(E))<2e-3:break

    final_action=apply(base,p);comp,st=reemit(final_action)
    vf,tf,kf=cr2_vector(st,red,inner);fin=worst(vf,tf,u)
    bg=background_residuals(comp);mask=u<INNER_HI
    bgmax={k:float(np.max(np.abs(np.asarray(v)[mask]))) for k,v in bg.items()}
    prod=(u>0.62)&(u<0.70);shift=0.
    for n in DIRS:
        shift=max(shift,float(np.max(np.abs(final_action[n].to_numpy(float)[prod]-base[n].to_numpy(float)[prod]))))
    if fin["min_cr2"]>0 and kf>0 and all(bgmax[k]<2e-3 for k in ("E00","E11","JA","eq131_residual")):
        status="INNER_NATIVE_CR2_FINITE_CONTINUATION_PASS"
    elif fin["min_cr2"]>init_full["min_cr2"]:
        status="INNER_NATIVE_CR2_IMPROVED_NOT_CLOSED"
    else: status="INNER_NATIVE_CR2_CONTINUATION_STALLED"
    rep={"scope":"finite inner action continuation; no promotion","inner_domain":[float(u[inner].min()),INNER_HI],
         "initial_full":init_full,"final_full":fin,"final_min_K":kf,"background_max_abs_inner":bgmax,
         "production_action_max_abs_shift":shift,"history":hist,"status":status}
    OUT.parent.mkdir(parents=True,exist_ok=True);OUT.write_text(json.dumps(rep,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"status":status,"initial":init_full,"final":fin,"minK":kf,"bg":bgmax},indent=2))
    return 0

if __name__=="__main__":raise SystemExit(main())
