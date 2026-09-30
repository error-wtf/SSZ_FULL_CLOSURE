#!/usr/bin/env python3
"""Inner-tail radial-characteristic controllability audit.

The current member has positive K on u<=0.62 but a nearly L-independent
negative minimum radial characteristic.  This audit asks whether the same
background/on-shell action freedom used in the outer-tail analysis contains
local directions that preserve E00,E11,JA,Eq131 to first order while increasing
the weakest generalized radial eigenvalue c_r^2 for all required L.

No finite repair or production promotion is performed here.
"""

from __future__ import annotations
import json, sys
from pathlib import Path
import itertools
import numpy as np
from scipy.optimize import linprog

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))

from ssz_p5.numerics import module  # noqa:E402
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central  # noqa:E402
from ssz_p5.production.full_action_lower import complete_total_action_jets,emit_lower_slots  # noqa:E402
from ssz_p5.production.central_action import background_residuals  # noqa:E402

OUT=ROOT/"data/generated/spectral/CURRENT_MEMBER_INNER_RADIAL_PROJECTED_CONTROLLABILITY.json"
DIRECTIONS=("f2","f2X","f2F","f3","f3X","f4","f4X","f4XX")
LS=(6,12,20,42,110,420,1000)
POINTS=(0.6190,0.6160,0.6130,0.6106)
EPS=1e-5
WIDTH=6e-4

def bump(u,u0):
    z=(np.asarray(u,float)-u0)/WIDTH
    q=np.zeros_like(z); m=np.abs(z)<1
    q[m]=np.exp(-1/(1-z[m]**2))/np.exp(-1)
    return q

def reemit(action):
    zk=module("ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py")
    d=action.copy().sort_values("x").reset_index(drop=True)
    r=d.x.to_numpy(float); f=d.f.to_numpy(float); h=d.h.to_numpy(float)
    X=d.X.to_numpy(float); A=d.A0prime.to_numpy(float); ph=d.phiprime.to_numpy(float)
    F=h*A*A/(2*f); Y=4*X*F
    Xp=zk.dr(r,X,1,9,8); Fp=zk.dr(r,F,1,9,8); Yp=zk.dr(r,Y,1,9,8)
    d["f2phi"]=(zk.dr(r,d.f2.to_numpy(float),1,9,8)-d.f2X.to_numpy(float)*Xp-d.f2F.to_numpy(float)*Fp-d.f2Y.to_numpy(float)*Yp)/ph
    d=complete_total_action_jets(d)
    lower,d=emit_lower_slots(d)
    out=zk.emit(d,selected_v5=lower.v5.to_numpy(float),selected_c3=lower.c3.to_numpy(float),selected_e3=lower.e3.to_numpy(float),v6_phi_selector="action")
    out["a5"]=zk.dr(r,out.a2.to_numpy(float),1,9,8)-zk.dr(r,out.a1.to_numpy(float),2,9,8)-zk.dr(r,out.A0prime.to_numpy(float)*out.v4.to_numpy(float)/2,1,9,8)+out.A0prime.to_numpy(float)*out.v5.to_numpy(float)/2
    return d,out

def bgvec(comp,i):
    b=background_residuals(comp)
    return np.array([b["E00"][i],b["E11"][i],b["JA"][i],b["eq131_residual"][i]],float)

def metrics(stream,red,L,i):
    a=red.canonical_audit(stream,int(L))
    K=np.asarray(a["K"],float); G=np.asarray(a["G"],float)
    Ks=(K[i]+K[i].T)/2; Gs=(G[i]+G[i].T)/2
    w,U=np.linalg.eigh(Ks)
    kmin=float(w[0])
    if kmin<=0:
        return kmin,np.nan
    inv=U@np.diag(1/np.sqrt(w))@U.T
    C=inv@Gs@inv
    cr=float(np.linalg.eigvalsh((C+C.T)/2)[0])
    return kmin,cr

def nullspace(B):
    U,S,Vt=np.linalg.svd(B,full_matrices=True)
    tol=max(B.shape)*np.finfo(float).eps*(S[0] if len(S) else 1.0)*100
    rank=int(np.sum(S>tol))
    return Vt[rank:].T,rank,[float(x) for x in S]

def common_lp(B,gcr,gk):
    n=B.shape[1]
    rn=np.linalg.norm(B,axis=1); keep=rn>1e-14
    Be=B[keep]/rn[keep,None]
    c=np.zeros(n+1); c[-1]=-1
    Aub=[]; bub=[]
    # Require every L radial derivative >= t.
    for L in LS:
        row=np.zeros(n+1); row[:n]=-gcr[L]; row[-1]=1
        Aub.append(row); bub.append(0.)
        # Also do not decrease kmin to first order.
        row2=np.zeros(n+1); row2[:n]=-gk[L]
        Aub.append(row2); bub.append(0.)
    Aeq=np.column_stack([Be,np.zeros(len(Be))]) if len(Be) else None
    sol=linprog(c,A_ub=np.asarray(Aub),b_ub=np.asarray(bub),A_eq=Aeq,b_eq=np.zeros(len(Be)) if len(Be) else None,
                bounds=[(-1,1)]*n+[(None,None)],method="highs")
    if not sol.success:
        return {"success":False,"message":sol.message}
    dp=sol.x[:n]
    return {
        "success":True,
        "common_cr2_gain":float(sol.x[-1]),
        "all_positive":bool(sol.x[-1]>0),
        "cr2_gains":{str(L):float(gcr[L]@dp) for L in LS},
        "K_gains":{str(L):float(gk[L]@dp) for L in LS},
        "direction":{DIRECTIONS[j]:float(dp[j]) for j in range(n) if abs(dp[j])>1e-8},
    }

def point(base,u0):
    loc=base[(base.u>=u0-0.0025)&(base.u<=u0+0.0025)].copy().sort_values("x").reset_index(drop=True)
    comp0,st0=reemit(loc); u=loc.u.to_numpy(float); i=int(np.argmin(np.abs(u-u0)))
    shape=bump(u,u[i]); E0=bgvec(comp0,i)
    red=module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    base_m={L:metrics(st0,red,L,i) for L in LS}
    B=np.zeros((4,len(DIRECTIONS)))
    gcr={L:np.zeros(len(DIRECTIONS)) for L in LS}
    gk={L:np.zeros(len(DIRECTIONS)) for L in LS}
    for j,name in enumerate(DIRECTIONS):
        d=loc.copy(); d[name]=d[name].to_numpy(float)+EPS*shape
        comp,st=reemit(d)
        B[:,j]=(bgvec(comp,i)-E0)/EPS
        for L in LS:
            km,cr=metrics(st,red,L,i)
            gk[L][j]=(km-base_m[L][0])/EPS
            gcr[L][j]=(cr-base_m[L][1])/EPS
    N,rank,S=nullspace(B)
    per={}
    for L in LS:
        pg=N@(N.T@gcr[L]) if N.size else np.zeros(len(DIRECTIONS))
        per[str(L)]={
            "kmin":base_m[L][0],"cr2":base_m[L][1],
            "cr2_gradient_norm":float(np.linalg.norm(gcr[L])),
            "projected_cr2_gradient_norm":float(np.linalg.norm(pg)),
            "structurally_controllable":bool(np.linalg.norm(pg)>1e-10*max(1.,np.linalg.norm(gcr[L]))),
        }
    return {
        "u_requested":u0,"u_grid":float(u[i]),"background":dict(zip(("E00","E11","JA","Eq131"),map(float,E0))),
        "rank_B":rank,"nullity":len(DIRECTIONS)-rank,"singular_values_B":S,
        "per_L":per,"simultaneous_all_L":common_lp(B,gcr,gk),
    }

def main():
    base=build_onshell_central(ROOT).action.sort_values("x").reset_index(drop=True)
    pts=[point(base,u0) for u0 in POINTS]
    all_ctrl=all(all(q["structurally_controllable"] for q in p["per_L"].values()) for p in pts)
    all_common=all(p["simultaneous_all_L"].get("all_positive",False) for p in pts)
    if all_common:
        diagnosis="INNER_RADIAL_EQ131_COMPATIBLE_COMMON_TANGENT_EXISTS"
        nxt="BUILD_FINITE_INNER_RADIAL_PREDICTOR_CORRECTOR"
    elif all_ctrl:
        diagnosis="INNER_RADIAL_PER_L_TANGENTS_EXIST_NO_COMMON_DIRECTION"
        nxt="EXPAND_INNER_ACTION_BASIS_OR_MULTIPOINT_CONTINUATION"
    else:
        diagnosis="INNER_RADIAL_EQ131_CONSTRAINED_CONTROLLABILITY_LOST"
        nxt="EXPAND_ACTION_SPACE_OR_REVISIT_INNER_BRANCH"
    out={
        "scope":"local first-order inner-tail radial controllability; no promotion",
        "points":pts,
        "all_points_per_L_structurally_controllable":all_ctrl,
        "all_points_common_all_L_positive_direction":all_common,
        "diagnosis":diagnosis,"next_target":nxt,
    }
    OUT.parent.mkdir(parents=True,exist_ok=True); OUT.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
    print(json.dumps(out,indent=2,allow_nan=False)); return 0

if __name__=="__main__": raise SystemExit(main())
