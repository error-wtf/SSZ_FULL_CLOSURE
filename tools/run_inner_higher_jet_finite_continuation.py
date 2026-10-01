#!/usr/bin/env python3
"""Finite inner continuation in the higher-action-jet sector.

Uses f3XX, f4XX and f4XXX localized on the native inner tail.  The projected
audit found a common radial-improving tangent at every tested inner point while
linearized E00/E11/JA/Eq131 constraints were imposed.

This finite predictor-corrector:
  * modifies only u<=0.62;
  * re-completes total-action phi jets and re-emits the 41-slot stream;
  * enforces E00,E11,JA,Eq131 at all control centers each iteration;
  * maximizes the weakest c_r^2 for every required L;
  * keeps K strictly positive;
  * verifies the complete native inner tail and leaves production untouched.

Diagnostic only; no promotion.
"""
from __future__ import annotations
import json,sys
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

OUT=ROOT/"data/generated/spectral/INNER_HIGHER_JET_FINITE_CONTINUATION.json"
ACTION_OUT=ROOT/"data/generated/spectral/INNER_HIGHER_JET_FINITE_CONTINUATION_ACTION.csv"
STREAM_OUT=ROOT/"data/generated/spectral/INNER_HIGHER_JET_FINITE_CONTINUATION_DIRECT41.csv"

JETS=("f3XX","f4XX","f4XXX")
CENTERS=(0.6100,0.6108,0.6120,0.6135,0.6150,0.6165,0.6180,0.61935)
WIDTH=8.5e-4
LS=tuple(int(x) for x in DEFAULT_L)
EPS=2e-5
MAX_IT=10
TRUST=8.0
INNER_HI=0.62


def bump(u,c):
    z=(np.asarray(u,float)-c)/WIDTH
    q=np.zeros_like(z);m=np.abs(z)<1
    q[m]=np.exp(-1/(1-z[m]**2))/np.exp(-1)
    return np.where(u<=INNER_HI,q,0.0)


def apply(base,p):
    d=base.copy();u=d.u.to_numpy(float)
    pp=np.asarray(p,float).reshape(len(CENTERS),len(JETS))
    for ic,c in enumerate(CENTERS):
        sh=bump(u,c)
        for j,name in enumerate(JETS):
            d[name]=d[name].to_numpy(float)+pp[ic,j]*sh
    return d


def reemit(action):
    zk=module("ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py")
    d=action.copy().sort_values("x").reset_index(drop=True)
    r=d.x.to_numpy(float);f=d.f.to_numpy(float);h=d.h.to_numpy(float)
    X=d.X.to_numpy(float);A=d.A0prime.to_numpy(float);ph=d.phiprime.to_numpy(float)
    F=h*A*A/(2*f);Y=4*X*F
    Xp=zk.dr(r,X,1,9,8);Fp=zk.dr(r,F,1,9,8);Yp=zk.dr(r,Y,1,9,8)
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


def constraints(comp,ids):
    b=background_residuals(comp)
    vals=[]
    for i in ids:
        vals.extend((b["E00"][i],b["E11"][i],b["JA"][i],b["eq131_residual"][i]))
    return np.asarray(vals,float)


def audit(stream,red,rows):
    cr=[];ct=[];kv=[];kt=[]
    for L in LS:
        a=red.canonical_audit(stream,int(L))
        K=np.asarray(a["K"],float);G=np.asarray(a["G"],float)
        for i in rows:
            Ks=(K[i]+K[i].T)/2;Gs=(G[i]+G[i].T)/2
            w,U=np.linalg.eigh(Ks);km=float(w[0])
            kv.append(km);kt.append((L,int(i)))
            if km<=0:
                cv=-1e6
            else:
                inv=U@np.diag(1/np.sqrt(w))@U.T
                C=inv@Gs@inv
                cv=float(np.linalg.eigvalsh((C+C.T)/2)[0])
            cr.append(cv);ct.append((L,int(i)))
    return np.asarray(cr),ct,np.asarray(kv),kt


def worst(vals,tags,u,key):
    j=int(np.argmin(vals));L,i=tags[j]
    return {key:float(vals[j]),"L":int(L),"u":float(u[i]),"row":int(i)}


def main():
    b=build_onshell_central(ROOT)
    base=b.action.sort_values("x").reset_index(drop=True)
    u=base.u.to_numpy(float)
    red=module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    inner=np.flatnonzero(u<=INNER_HI)
    colloc=np.array([int(np.argmin(np.abs(u-c))) for c in CENTERS],int)
    pick=np.linspace(0,len(inner)-1,min(64,len(inner))).round().astype(int)
    sample=np.unique(np.concatenate([inner[pick],colloc,[inner[0],inner[-1]]]))

    npar=len(CENTERS)*len(JETS)
    p=np.zeros(npar)
    comp,stream=reemit(apply(base,p))
    E=constraints(comp,colloc)
    cr,ct,kv,kt=audit(stream,red,sample)
    init_cr=worst(cr,ct,u,"min_cr2");init_k=worst(kv,kt,u,"min_K")
    hist=[]

    for it in range(MAX_IT):
        B=np.zeros((len(E),npar));Gc=np.zeros((len(cr),npar));Gk=np.zeros((len(kv),npar))
        for j in range(npar):
            pp=p.copy();pp[j]+=EPS
            cp,sp=reemit(apply(base,pp))
            Ep=constraints(cp,colloc)
            c2,_,k2,_=audit(sp,red,sample)
            B[:,j]=(Ep-E)/EPS
            Gc[:,j]=(c2-cr)/EPS
            Gk[:,j]=(k2-kv)/EPS

        rn=np.linalg.norm(B,axis=1);keep=rn>1e-12
        Be=B[keep]/rn[keep,None] if np.any(keep) else np.zeros((0,npar))
        ee=E[keep]/rn[keep] if np.any(keep) else np.zeros(0)

        obj=np.zeros(npar+1);obj[-1]=-1
        Aub=[];bub=[]
        for q in range(len(cr)):
            row=np.zeros(npar+1);row[:-1]=-Gc[q];row[-1]=1
            Aub.append(row);bub.append(float(cr[q]))
        for q in range(len(kv)):
            row=np.zeros(npar+1);row[:-1]=-Gk[q]
            Aub.append(row);bub.append(float(kv[q]-1e-8))
        Aeq=np.column_stack([Be,np.zeros(len(Be))]) if len(Be) else None
        beq=-ee if len(Be) else None
        sol=linprog(obj,A_ub=np.asarray(Aub),b_ub=np.asarray(bub),
                    A_eq=Aeq,b_eq=beq,
                    bounds=[(-TRUST,TRUST)]*npar+[(None,None)],method="highs")
        if not sol.success:
            hist.append({"iteration":it,"status":"LP_FAIL","message":sol.message});break

        dp=sol.x[:-1];old=float(np.min(cr));best=None;accepted=False
        for frac in (1,.5,.25,.125,.0625,.03125,.015625,.0078125):
            pt=p+frac*dp
            ct0,st0=reemit(apply(base,pt))
            Et=constraints(ct0,colloc)
            c2,c2t,k2,k2t=audit(st0,red,sample)
            wi=worst(c2,c2t,u,"min_cr2");wk=worst(k2,k2t,u,"min_K")
            em=float(np.max(np.abs(Et)))
            rec={"iteration":it,"fraction":frac,"predicted_margin":float(sol.x[-1]),
                 "worst_cr2":wi,"worst_K":wk,"constraint_max_abs":em,
                 "step_max_abs":float(np.max(np.abs(frac*dp)))}
            score=wi["min_cr2"]-2*em
            if best is None or score>best[0]:
                best=(score,pt,ct0,st0,Et,c2,c2t,k2,k2t,rec)
            if wi["min_cr2"]>old+1e-5 and wk["min_K"]>1e-8 and em<2e-3:
                p,comp,stream,E,cr,ct,kv,kt=pt,ct0,st0,Et,c2,c2t,k2,k2t
                rec["status"]="ACCEPT";hist.append(rec);accepted=True;break
        if not accepted:
            if best and best[-1]["worst_cr2"]["min_cr2"]>old and best[-1]["worst_K"]["min_K"]>1e-8 and best[-1]["constraint_max_abs"]<5e-3:
                _,p,comp,stream,E,cr,ct,kv,kt,rec=best
                rec["status"]="ACCEPT_BEST";hist.append(rec);accepted=True
            else:
                if best:
                    best[-1]["status"]="STALL";hist.append(best[-1])
                break

        fc,fct,fk,fkt=audit(stream,red,inner)
        if np.min(fc)>0 and np.min(fk)>0 and np.max(np.abs(E))<2e-3:
            break

    final_action=apply(base,p)
    comp,stream=reemit(final_action)
    fc,fct,fk,fkt=audit(stream,red,inner)
    final_cr=worst(fc,fct,u,"min_cr2");final_k=worst(fk,fkt,u,"min_K")
    bg=background_residuals(comp)
    im=u<=INNER_HI
    bgmax={k:float(np.max(np.abs(np.asarray(v)[im]))) for k,v in bg.items()}

    prod=(u>0.62)&(u<0.70)
    prod_shift=max(float(np.max(np.abs(final_action[j].to_numpy(float)[prod]-base[j].to_numpy(float)[prod]))) for j in JETS)

    bg_ok=all(bgmax[k]<2e-3 for k in ("E00","E11","JA","eq131_residual"))
    if final_cr["min_cr2"]>0 and final_k["min_K"]>0 and bg_ok:
        status="INNER_HIGHER_JET_FINITE_CONTINUATION_PASS"
        nxt="COMBINE_WITH_OUTER_REPAIR_AND_BUILD_FULL_NATIVE_CANDIDATE"
    elif final_cr["min_cr2"]>init_cr["min_cr2"]:
        status="INNER_HIGHER_JET_FINITE_CONTINUATION_IMPROVED_NOT_CLOSED"
        nxt="COMBINE_HIGHER_JETS_WITH_F2_HESSIAN_AND_NONLINEAR_ROOT_CORRECTOR"
    else:
        status="INNER_HIGHER_JET_FINITE_CONTINUATION_STALLED"
        nxt="BUILD_COUPLED_HIGHER_JET_PLUS_F2_HESSIAN_CORRECTOR"

    OUT.parent.mkdir(parents=True,exist_ok=True)
    comp.to_csv(ACTION_OUT,index=False);stream.to_csv(STREAM_OUT,index=False)
    report={
      "status":status,"next_target":nxt,
      "scope":"finite f3XX/f4XX/f4XXX inner continuation with full re-emission",
      "initial_worst_cr2":init_cr,"initial_worst_K":init_k,
      "final_worst_cr2":final_cr,"final_worst_K":final_k,
      "background_max_abs_inner":bgmax,
      "production_higher_jet_max_abs_shift":prod_shift,
      "centers":list(CENTERS),"width":WIDTH,"directions":list(JETS),
      "parameter_max_abs":float(np.max(np.abs(p))),
      "nonzero_parameter_count":int(np.sum(np.abs(p)>1e-10)),
      "iterations":hist,
      "artifacts":{
        "action_csv":str(ACTION_OUT.relative_to(ROOT)),
        "direct41_csv":str(STREAM_OUT.relative_to(ROOT)),
      }
    }
    OUT.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps(report,indent=2,allow_nan=False))
    return 0

if __name__=="__main__":raise SystemExit(main())
