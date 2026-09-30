#!/usr/bin/env python3
"""Finite Eq131-compatible predictor-corrector trial on the native outer tail.

This is the finite continuation licensed by the preceding projected
controllability audit.  It is deliberately conservative:

* only u>0.7003 is modified, so the registered production window is untouched;
* action-level primitives are varied with compact C-infinity bumps;
* after every finite predictor step the full action is re-completed and all
  41 slots are absolutely re-emitted;
* E00, E11, JA and Eq131 are used as nonlinear corrector constraints;
* the weakest kinetic eigenvalue across sampled tail rows and required L is the
  continuation objective;
* the tangent/Jacobian is rebuilt after every accepted step.

This is a diagnostic finite-member search.  No result is promoted or frozen as
the production member by this script.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ssz_p5.numerics import module  # noqa:E402
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central  # noqa:E402
from ssz_p5.production.full_action_lower import complete_total_action_jets, emit_lower_slots  # noqa:E402
from ssz_p5.production.central_action import background_residuals  # noqa:E402

OUT = ROOT / "data/generated/spectral/CURRENT_MEMBER_EQ131_FINITE_CONTINUATION_TRIAL.json"

DIRECTIONS = ("f2", "f2X", "f2F", "f3", "f3X", "f4", "f4X", "f4XX")
CENTERS = (0.7013, 0.7030, 0.7048, 0.7062)
WIDTH = 0.00115
LS = (6, 20, 42, 110, 420, 1000)
FD_EPS = 2e-5
MAX_IT = 7
TRUST = 0.08
TAIL_LO = 0.70035
TAIL_HI = 0.7080


def bump(u, u0):
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


def apply_params(base, p):
    d=base.copy()
    u=d.u.to_numpy(float)
    p=np.asarray(p,float).reshape(len(CENTERS),len(DIRECTIONS))
    for ic,c in enumerate(CENTERS):
        sh=bump(u,c)
        for j,name in enumerate(DIRECTIONS):
            d[name]=d[name].to_numpy(float)+p[ic,j]*sh
    return d


def constraint_vector(completed, colloc):
    bg=background_residuals(completed)
    vals=[]
    for i in colloc:
        vals.extend((bg["E00"][i],bg["E11"][i],bg["JA"][i],bg["eq131_residual"][i]))
    return np.asarray(vals,float)


def kinetic_vector(stream, red, sample):
    vals=[]
    tags=[]
    for L in LS:
        a=red.canonical_audit(stream,int(L))
        K=np.asarray(a["K"],float)
        Ks=(K+K.swapaxes(1,2))/2
        ev=np.linalg.eigvalsh(Ks)[:,0]
        for i in sample:
            vals.append(ev[i]); tags.append((L,i))
    return np.asarray(vals,float),tags


def worst_info(vals,tags,u):
    q=int(np.argmin(vals)); L,i=tags[q]
    return {"min_K":float(vals[q]),"L":int(L),"u":float(u[i]),"row":int(i)}


def main():
    red=module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    build=build_onshell_central(ROOT)
    base=build.action.sort_values("x").reset_index(drop=True)
    u=base.u.to_numpy(float)

    # Stay strictly outside production.
    active=np.flatnonzero((u>=TAIL_LO)&(u<=TAIL_HI))
    if len(active)<40:
        raise RuntimeError("insufficient tail rows")

    colloc=np.array([int(np.argmin(np.abs(u-c))) for c in CENTERS],dtype=int)
    # Uniformly subsample kinetic rows while always including collocation rows.
    pick=np.linspace(0,len(active)-1,48).round().astype(int)
    sample=np.unique(np.concatenate([active[pick],colloc]))

    npar=len(CENTERS)*len(DIRECTIONS)
    p=np.zeros(npar)
    history=[]

    comp,stream=reemit(apply_params(base,p))
    E=constraint_vector(comp,colloc)
    kvals,tags=kinetic_vector(stream,red,sample)
    initial=worst_info(kvals,tags,u)

    for it in range(MAX_IT):
        # Finite-difference nonlinear Jacobians around current accepted point.
        B=np.zeros((len(E),npar))
        G=np.zeros((len(kvals),npar))
        for j in range(npar):
            pp=p.copy(); pp[j]+=FD_EPS
            cp,sp=reemit(apply_params(base,pp))
            Ep=constraint_vector(cp,colloc)
            kp,_=kinetic_vector(sp,red,sample)
            B[:,j]=(Ep-E)/FD_EPS
            G[:,j]=(kp-kvals)/FD_EPS

        # Row-equilibrate equality corrector for numerical conditioning.
        rn=np.linalg.norm(B,axis=1)
        keep=rn>1e-12
        Be=B[keep]/rn[keep,None]
        ee=E[keep]/rn[keep]

        # x=(dp,t): maximize t subject to finite-linearized K>=t and
        # linearized background correction E+B dp=0.
        obj=np.zeros(npar+1); obj[-1]=-1
        Aub=np.column_stack([-G,np.ones(len(kvals))])
        bub=kvals.copy()
        Aeq=np.column_stack([Be,np.zeros(len(Be))]) if len(Be) else None
        beq=-ee if len(Be) else None
        bounds=[(-TRUST,TRUST)]*npar+[(None,None)]
        sol=linprog(obj,A_ub=Aub,b_ub=bub,A_eq=Aeq,b_eq=beq,bounds=bounds,method="highs")
        if not sol.success:
            history.append({"iteration":it,"status":"LP_FAIL","message":sol.message})
            break

        dp=sol.x[:-1]
        curmin=float(kvals.min())
        accepted=False
        best=None
        for frac in (1.0,0.5,0.25,0.125,0.0625):
            pt=p+frac*dp
            ct,st=reemit(apply_params(base,pt))
            Et=constraint_vector(ct,colloc)
            kt,tagt=kinetic_vector(st,red,sample)
            emax=float(np.max(np.abs(Et)))
            wi=worst_info(kt,tagt,u)
            rec={
                "iteration":it,
                "fraction":frac,
                "predicted_margin":float(sol.x[-1]),
                "finite_worst":wi,
                "constraint_max_abs":emax,
                "step_max_abs":float(np.max(np.abs(frac*dp))),
            }
            score=wi["min_K"]-10.0*emax
            if best is None or score>best[0]:
                best=(score,pt,ct,st,Et,kt,tagt,rec)
            # conservative acceptance: improve K and keep constraints small.
            if wi["min_K"]>curmin+max(1e-6,1e-4*max(1.0,abs(curmin))) and emax<2e-3:
                p=pt; comp=ct; stream=st; E=Et; kvals=kt; tags=tagt
                rec["status"]="ACCEPT"
                history.append(rec); accepted=True
                break
        if not accepted:
            if best and best[7]["constraint_max_abs"]<5e-3 and best[7]["finite_worst"]["min_K"]>curmin:
                _,p,comp,stream,E,kvals,tags,rec=best
                rec["status"]="ACCEPT_BEST"
                history.append(rec); accepted=True
            else:
                if best:
                    best[7]["status"]="STALL"
                    history.append(best[7])
                break
        if float(kvals.min())>0 and float(np.max(np.abs(E)))<2e-3:
            break

    # Full-row verification over active native tail, not just sampled rows.
    full_vals,full_tags=kinetic_vector(stream,red,active)
    final=worst_info(full_vals,full_tags,u)
    bg=background_residuals(comp)
    region=(u>=TAIL_LO)&(u<=TAIL_HI)
    bgmax={k:float(np.max(np.abs(np.asarray(v)[region]))) for k,v in bg.items()}

    prod=(u>0.62)&(u<0.70)
    # verify production stream was untouched at action level
    final_action=apply_params(base,p)
    prod_shift=0.0
    for name in DIRECTIONS:
        prod_shift=max(prod_shift,float(np.max(np.abs(final_action[name].to_numpy(float)[prod]-base[name].to_numpy(float)[prod]))))

    if final["min_K"]>0 and bgmax["eq131_residual"]<2e-3 and bgmax["E00"]<2e-3 and bgmax["E11"]<2e-3 and bgmax["JA"]<2e-3:
        status="FINITE_EQ131_K_CONTINUATION_PASS"
        nxt="RUN_FULL_G_RADIAL_AND_INTERFACE_PROVENANCE_ON_FINITE_TRIAL"
    elif final["min_K"]>initial["min_K"]:
        status="FINITE_EQ131_K_CONTINUATION_IMPROVES_BUT_NOT_CLOSED"
        nxt="REFINE_BASIS_TRUST_AND_NONLINEAR_CORRECTOR"
    else:
        status="FINITE_EQ131_K_CONTINUATION_STALLED"
        nxt="EXPAND_ACTION_BASIS_OR_USE_NONLINEAR_ROOT_CORRECTOR"

    report={
        "scope":"finite predictor-corrector trial; no promotion/no production overwrite",
        "tail_domain":[TAIL_LO,TAIL_HI],
        "directions":list(DIRECTIONS),
        "centers":list(CENTERS),
        "width":WIDTH,
        "iterations":history,
        "initial_worst":initial,
        "final_worst_full_active_rows":final,
        "background_max_abs_active_region":bgmax,
        "production_action_max_abs_shift":prod_shift,
        "parameter_max_abs":float(np.max(np.abs(p))),
        "nonzero_parameter_count":int(np.sum(np.abs(p)>1e-10)),
        "status":status,
        "next_target":nxt,
    }
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps(report,indent=2,allow_nan=False))
    return 0

if __name__=="__main__":
    raise SystemExit(main())
