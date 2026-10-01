#!/usr/bin/env python3
"""Coupled inner-tail corrector: higher action jets + background-null f2 Hessian.

Combines the two independently admissible control sectors already audited:
  A) higher total-action jets f3XX, f4XX, f4XXX, with full action re-emission;
  B) exact additive background-null f2 Hessian controls H t = 0, added through
     their Appendix-A 41-slot response.

The higher-jet-only finite continuation is allowed to run independently.  This
solver is the next nonlinear closure attempt if that sector alone is
insufficient.  It maximizes the weakest radial c_r^2 over the complete native
inner tail while preserving K>0 and actively correcting E00,E11,JA,Eq131.

No production promotion is performed.
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
from ssz_p5.production.holonomic_hessian_y import TRANSVERSE,complete_hessian,response_from_hessian  # noqa:E402

OUT=ROOT/"data/generated/spectral/INNER_COUPLED_HIGHERJET_HESSIAN_CORRECTOR.json"
ACTION_OUT=ROOT/"data/generated/spectral/INNER_COUPLED_HIGHERJET_HESSIAN_CORRECTOR_ACTION.csv"
STREAM_OUT=ROOT/"data/generated/spectral/INNER_COUPLED_HIGHERJET_HESSIAN_CORRECTOR_DIRECT41.csv"

HJETS=("f3XX","f4XX","f4XXX")
HCENTERS=(0.6101,0.6122,0.6145,0.6167,0.6187)
JCENTERS=(0.6100,0.6112,0.6126,0.6140,0.6155,0.6170,0.6184,0.61945)
HWIDTH=0.00155
JWIDTH=0.0009
LS=tuple(int(x) for x in DEFAULT_L)
INNER_HI=0.62
EPS_J=2e-5
EPS_H=5e-4
MAX_IT=8
TRUST_J=6.0
TRUST_H=1.0


def bump(u,c,w):
    z=(np.asarray(u,float)-c)/w
    q=np.zeros_like(z);m=np.abs(z)<1
    q[m]=np.exp(-1/(1-z[m]**2))/np.exp(-1)
    return np.where(u<=INNER_HI,q,0.0)


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
    st=zk.emit(d,
        selected_v5=lower.v5.to_numpy(float),
        selected_c3=lower.c3.to_numpy(float),
        selected_e3=lower.e3.to_numpy(float),
        v6_phi_selector="action")
    st["a5"]=(zk.dr(r,st.a2.to_numpy(float),1,9,8)
              -zk.dr(r,st.a1.to_numpy(float),2,9,8)
              -zk.dr(r,st.A0prime.to_numpy(float)*st.v4.to_numpy(float)/2,1,9,8)
              +st.A0prime.to_numpy(float)*st.v5.to_numpy(float)/2)
    return d,st


def hessian_columns(action):
    u=action.u.to_numpy(float)
    cols=[];names=[]
    for c in HCENTERS:
        sh=bump(u,c,HWIDTH)
        for j,name in enumerate(TRANSVERSE):
            q=np.zeros((len(action),6),float);q[:,j]=sh
            H,_=complete_hessian(action,q)
            resp=response_from_hessian(action,H)
            cols.append({slot:resp[slot].to_numpy(float) for slot in ("v5","c3","e3","v1","v4","c2")})
            names.append(f"{name}@{c:.4f}")
    return cols,names


def apply_hjets(base,p):
    d=base.copy();u=d.u.to_numpy(float)
    pp=np.asarray(p,float).reshape(len(JCENTERS),len(HJETS))
    for ic,c in enumerate(JCENTERS):
        sh=bump(u,c,JWIDTH)
        for j,name in enumerate(HJETS):
            d[name]=d[name].to_numpy(float)+pp[ic,j]*sh
    return d


def add_hessian(stream,cols,p):
    d=stream.copy()
    for slot in ("v5","c3","e3","v1","v4","c2"):
        x=d[slot].to_numpy(float).copy()
        for a,col in zip(p,cols):
            if a:x+=float(a)*col[slot]
        d[slot]=x
    return d


def constraints(comp,ids):
    b=background_residuals(comp)
    vals=[]
    for i in ids: vals.extend((b["E00"][i],b["E11"][i],b["JA"][i],b["eq131_residual"][i]))
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
            if km<=0: cv=-1e6
            else:
                inv=U@np.diag(1/np.sqrt(w))@U.T
                C=inv@Gs@inv
                cv=float(np.linalg.eigvalsh((C+C.T)/2)[0])
            cr.append(cv);ct.append((L,int(i)))
    return np.asarray(cr),ct,np.asarray(kv),kt


def worst(vals,tags,u,key):
    j=int(np.argmin(vals));L,i=tags[j]
    return {key:float(vals[j]),"L":int(L),"u":float(u[i]),"row":int(i)}


def evaluate(base,jp,hp,hcols,red,rows,colloc):
    action=apply_hjets(base,jp)
    comp,st=reemit(action)
    st=add_hessian(st,hcols,hp)
    E=constraints(comp,colloc)
    cr,ct,k,kt=audit(st,red,rows)
    return action,comp,st,E,cr,ct,k,kt


def main():
    b=build_onshell_central(ROOT)
    base=b.action.sort_values("x").reset_index(drop=True)
    u=base.u.to_numpy(float)
    red=module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    hcols,hnames=hessian_columns(base)
    jnames=[f"{n}@{c:.4f}" for c in JCENTERS for n in HJETS]
    nj=len(jnames);nh=len(hnames)

    inner=np.flatnonzero(u<=INNER_HI)
    colloc=np.array([int(np.argmin(np.abs(u-c))) for c in sorted(set(JCENTERS+HCENTERS))],int)
    pick=np.linspace(0,len(inner)-1,min(58,len(inner))).round().astype(int)
    rows=np.unique(np.concatenate([inner[pick],colloc,[inner[0],inner[-1]]]))

    jp=np.zeros(nj);hp=np.zeros(nh)
    action,comp,st,E,cr,ct,kv,kt=evaluate(base,jp,hp,hcols,red,rows,colloc)
    initial_cr=worst(cr,ct,u,"min_cr2");initial_k=worst(kv,kt,u,"min_K")
    hist=[]

    for it in range(MAX_IT):
        npar=nj+nh
        B=np.zeros((len(E),npar));Gc=np.zeros((len(cr),npar));Gk=np.zeros((len(kv),npar))

        # Higher-jet finite differences.
        for j in range(nj):
            p2=jp.copy();p2[j]+=EPS_J
            _,cp,sp,Ep,c2,_,k2,_=evaluate(base,p2,hp,hcols,red,rows,colloc)
            B[:,j]=(Ep-E)/EPS_J;Gc[:,j]=(c2-cr)/EPS_J;Gk[:,j]=(k2-kv)/EPS_J

        # Hessian response is exactly background-null; only operator derivatives.
        for h in range(nh):
            q2=hp.copy();q2[h]+=EPS_H
            _,_,sp,Ep,c2,_,k2,_=evaluate(base,jp,q2,hcols,red,rows,colloc)
            j=nj+h
            B[:,j]=(Ep-E)/EPS_H;Gc[:,j]=(c2-cr)/EPS_H;Gk[:,j]=(k2-kv)/EPS_H

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
        bounds=[(-TRUST_J,TRUST_J)]*nj+[(-TRUST_H,TRUST_H)]*nh+[(None,None)]
        sol=linprog(obj,A_ub=np.asarray(Aub),b_ub=np.asarray(bub),
                    A_eq=np.column_stack([Be,np.zeros(len(Be))]) if len(Be) else None,
                    b_eq=-ee if len(Be) else None,bounds=bounds,method="highs")
        if not sol.success:
            hist.append({"iteration":it,"status":"LP_FAIL","message":sol.message});break

        dj=sol.x[:nj];dh=sol.x[nj:-1]
        old=float(np.min(cr));best=None;accepted=False
        for frac in (1,.5,.25,.125,.0625,.03125,.015625):
            j2=jp+frac*dj;h2=hp+frac*dh
            a2,c2,s2,E2,cc,cct,kk,kkt=evaluate(base,j2,h2,hcols,red,rows,colloc)
            wi=worst(cc,cct,u,"min_cr2");wk=worst(kk,kkt,u,"min_K");em=float(np.max(np.abs(E2)))
            rec={"iteration":it,"fraction":frac,"predicted_margin":float(sol.x[-1]),
                 "worst_cr2":wi,"worst_K":wk,"constraint_max_abs":em,
                 "higher_step_max":float(np.max(np.abs(frac*dj))),
                 "hessian_step_max":float(np.max(np.abs(frac*dh)))}
            score=wi["min_cr2"]-2*em
            if best is None or score>best[0]:best=(score,j2,h2,a2,c2,s2,E2,cc,cct,kk,kkt,rec)
            if wi["min_cr2"]>old+1e-5 and wk["min_K"]>1e-8 and em<2e-3:
                _,jp,hp,action,comp,st,E,cr,ct,kv,kt,_=best
                rec["status"]="ACCEPT";hist.append(rec);accepted=True;break
        if not accepted:
            if best and best[-1]["worst_cr2"]["min_cr2"]>old and best[-1]["worst_K"]["min_K"]>1e-8 and best[-1]["constraint_max_abs"]<5e-3:
                _,jp,hp,action,comp,st,E,cr,ct,kv,kt,rec=best
                rec["status"]="ACCEPT_BEST";hist.append(rec);accepted=True
            else:
                if best:best[-1]["status"]="STALL";hist.append(best[-1])
                break

        fc,fct,fk,fkt=audit(st,red,inner)
        if np.min(fc)>0 and np.min(fk)>0 and np.max(np.abs(E))<2e-3:break

    action,comp,st,E,_,_,_,_=evaluate(base,jp,hp,hcols,red,rows,colloc)
    fc,fct,fk,fkt=audit(st,red,inner)
    final_cr=worst(fc,fct,u,"min_cr2");final_k=worst(fk,fkt,u,"min_K")
    bg=background_residuals(comp);im=u<=INNER_HI
    bgmax={k:float(np.max(np.abs(np.asarray(v)[im]))) for k,v in bg.items()}
    prod=(u>0.62)&(u<0.70)
    prod_shift=max(float(np.max(np.abs(action[n].to_numpy(float)[prod]-base[n].to_numpy(float)[prod]))) for n in HJETS)

    bg_ok=all(bgmax[k]<2e-3 for k in ("E00","E11","JA","eq131_residual"))
    if final_cr["min_cr2"]>0 and final_k["min_K"]>0 and bg_ok:
        status="INNER_COUPLED_CORRECTOR_PASS";nxt="COMBINE_WITH_OUTER_REPAIR_AND_RUN_FULL_NATIVE_PRECERTIFICATE"
    elif final_cr["min_cr2"]>initial_cr["min_cr2"]:
        status="INNER_COUPLED_CORRECTOR_IMPROVED_NOT_CLOSED";nxt="ADAPTIVE_KNOT_REFINEMENT_AROUND_REMAINING_WORST_RADIUS"
    else:
        status="INNER_COUPLED_CORRECTOR_STALLED";nxt="REVISIT_INNER_INTERFACE_OR_ADDITIONAL_ACTION_SECTOR"

    OUT.parent.mkdir(parents=True,exist_ok=True)
    comp.to_csv(ACTION_OUT,index=False);st.to_csv(STREAM_OUT,index=False)
    report={
      "status":status,"next_target":nxt,
      "initial_worst_cr2":initial_cr,"initial_worst_K":initial_k,
      "final_worst_cr2":final_cr,"final_worst_K":final_k,
      "background_max_abs_inner":bgmax,
      "production_higher_jet_max_abs_shift":prod_shift,
      "higher_parameters":{jnames[i]:float(jp[i]) for i in range(nj) if abs(jp[i])>1e-10},
      "hessian_parameters":{hnames[i]:float(hp[i]) for i in range(nh) if abs(hp[i])>1e-10},
      "iterations":hist,
      "artifacts":{"action_csv":str(ACTION_OUT.relative_to(ROOT)),"direct41_csv":str(STREAM_OUT.relative_to(ROOT))}
    }
    OUT.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"status":status,"initial":initial_cr,"final":final_cr,"K":final_k,"background":bgmax,"next":nxt},indent=2))
    return 0

if __name__=="__main__":raise SystemExit(main())
