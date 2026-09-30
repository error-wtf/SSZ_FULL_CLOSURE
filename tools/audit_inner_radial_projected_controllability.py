#!/usr/bin/env python3
"""Inner-tail Eq131-compatible radial-characteristic controllability audit.

The inner c_r^2 crossing at u~=0.619021 is already present in the archival
pre-A2 source and is dominated by the psi characteristic eigenvector.  The
raw 41-slot sensitivity is strongest in fixed Einstein/background slots
(a1,a3,a4,b1,b4,b5,d1), so the relevant physical question is not whether
those slots can be edited, but whether the available action-jet nullspace can
compensate them while preserving the background equations.

At several inner-tail radii this audit computes action-level finite-difference
responses for directions
  f2,f2X,f2F,f3,f3X,f4,f4X,f4XX,
projects against
  E00=E11=JA=Eq131=0,
and asks whether a common constrained tangent can increase the minimum radial
characteristic c_r^2 for all tested L simultaneously.

No finite member is promoted.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))

from ssz_p5.numerics import module  # noqa:E402
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central  # noqa:E402
from ssz_p5.production.full_action_lower import complete_total_action_jets,emit_lower_slots  # noqa:E402
from ssz_p5.production.central_action import background_residuals  # noqa:E402

OUT=ROOT/"data/generated/spectral/INNER_RADIAL_PROJECTED_CONTROLLABILITY.json"
DIRECTIONS=("f2","f2X","f2F","f3","f3X","f4","f4X","f4XX")
LS=(6,20,42,110,420,1000)
POINTS=(0.6130,0.6160,0.6185)
WIDTH=7e-4
EPS=1e-5


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
    d["f2phi"]=(zk.dr(r,d.f2.to_numpy(float),1,9,8)
                -d.f2X.to_numpy(float)*Xp-d.f2F.to_numpy(float)*Fp
                -d.f2Y.to_numpy(float)*Yp)/ph
    d=complete_total_action_jets(d)
    lower,d=emit_lower_slots(d)
    out=zk.emit(d,selected_v5=lower.v5.to_numpy(float),selected_c3=lower.c3.to_numpy(float),
                selected_e3=lower.e3.to_numpy(float),v6_phi_selector="action")
    out["a5"]=(zk.dr(r,out.a2.to_numpy(float),1,9,8)-zk.dr(r,out.a1.to_numpy(float),2,9,8)
               -zk.dr(r,out.A0prime.to_numpy(float)*out.v4.to_numpy(float)/2,1,9,8)
               +out.A0prime.to_numpy(float)*out.v5.to_numpy(float)/2)
    return d,out


def constraint(comp,i):
    b=background_residuals(comp)
    return np.array([b["E00"][i],b["E11"][i],b["JA"][i],b["eq131_residual"][i]],float)


def cr2(stream,red,L,i):
    a=red.canonical_audit(stream,int(L))
    K=np.asarray(a["K"],float)[i];G=np.asarray(a["G"],float)[i]
    K=(K+K.T)/2;G=(G+G.T)/2
    w,U=np.linalg.eigh(K)
    if np.min(w)<=0:return None
    inv=U@np.diag(1/np.sqrt(w))@U.T
    C=inv@G@inv;C=(C+C.T)/2
    return float(np.linalg.eigvalsh(C)[0])


def nullspace(B):
    U,S,Vt=np.linalg.svd(B,full_matrices=True)
    tol=max(B.shape)*np.finfo(float).eps*(S[0] if len(S) else 1)*100
    rank=int(np.sum(S>tol))
    return Vt[rank:].T,rank,[float(z) for z in S]


def common_lp(B,grads):
    n=B.shape[1]
    rn=np.linalg.norm(B,axis=1);keep=rn>1e-14
    Be=B[keep]/rn[keep,None]
    c=np.zeros(n+1);c[-1]=-1
    Aub=[];bub=[]
    for L in LS:
        row=np.zeros(n+1);row[:n]=-grads[L];row[-1]=1
        Aub.append(row);bub.append(0)
    Aeq=np.column_stack([Be,np.zeros(len(Be))]) if len(Be) else None
    sol=linprog(c,A_ub=np.asarray(Aub),b_ub=np.asarray(bub),
                A_eq=Aeq,b_eq=np.zeros(len(Be)) if len(Be) else None,
                bounds=[(-1,1)]*n+[(None,None)],method="highs")
    if not sol.success:return {"success":False,"message":sol.message}
    dp=sol.x[:n];g={str(L):float(grads[L]@dp) for L in LS}
    return {
        "success":True,"common_gain":float(sol.x[-1]),
        "all_positive":bool(min(g.values())>0),"gains":g,
        "direction":{DIRECTIONS[j]:float(dp[j]) for j in range(n) if abs(dp[j])>1e-8},
    }


def point(base,u0):
    loc=base[(base.u>=u0-0.003)&(base.u<=u0+0.003)].copy().sort_values("x").reset_index(drop=True)
    comp0,s0=reemit(loc);u=loc.u.to_numpy(float);i=int(np.argmin(np.abs(u-u0)))
    red=module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    E0=constraint(comp0,i);basec={L:cr2(s0,red,L,i) for L in LS}
    B=np.zeros((4,len(DIRECTIONS)));gr={L:np.zeros(len(DIRECTIONS)) for L in LS}
    sh=bump(u,u[i])
    for j,name in enumerate(DIRECTIONS):
        d=loc.copy();d[name]=d[name].to_numpy(float)+EPS*sh
        cp,sp=reemit(d);B[:,j]=(constraint(cp,i)-E0)/EPS
        for L in LS:
            q=cr2(sp,red,L,i)
            gr[L][j]=(q-basec[L])/EPS
    N,rank,S=nullspace(B)
    per={}
    for L in LS:
        pg=N@(N.T@gr[L]) if N.size else np.zeros(len(DIRECTIONS))
        per[str(L)]={
            "cr2":basec[L],
            "gradient_norm":float(np.linalg.norm(gr[L])),
            "projected_norm":float(np.linalg.norm(pg)),
            "projected_fraction":float(np.linalg.norm(pg)/max(np.linalg.norm(gr[L]),1e-300)),
            "controllable":bool(np.linalg.norm(pg)>1e-10*max(1,np.linalg.norm(gr[L]))),
        }
    return {
        "u_requested":u0,"u_grid":float(u[i]),"constraints":E0.tolist(),
        "rank_B":rank,"nullity":len(DIRECTIONS)-rank,"singular_values":S,
        "per_L":per,"simultaneous":common_lp(B,gr),
    }


def main():
    base=build_onshell_central(ROOT).action.sort_values("x").reset_index(drop=True)
    pts=[point(base,u0) for u0 in POINTS]
    all_common=all(p["simultaneous"].get("all_positive",False) for p in pts)
    all_ind=all(all(v["controllable"] for v in p["per_L"].values()) for p in pts)
    if all_common:
        diagnosis="INNER_RADIAL_EQ131_COMPATIBLE_COMMON_TANGENT_EXISTS"
        nxt="BUILD_FINITE_INNER_RADIAL_PREDICTOR_CORRECTOR"
    elif all_ind:
        diagnosis="INNER_RADIAL_PER_L_TANGENTS_EXIST_NO_COMMON_DIRECTION"
        nxt="EXPAND_INNER_ACTION_BASIS_OR_MULTIPOINT_CONTINUATION"
    else:
        diagnosis="INNER_RADIAL_ACTION_NULLSPACE_LOSES_CONTROLLABILITY"
        nxt="TEST_GEOMETRY_FIXED_EINSTEIN_CHANNEL_COMPATIBILITY"
    rep={"scope":"inner radial action-level projected controllability; no finite promotion",
         "points":pts,"diagnosis":diagnosis,"next_target":nxt}
    OUT.parent.mkdir(parents=True,exist_ok=True);OUT.write_text(json.dumps(rep,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"diagnosis":diagnosis,"next_target":nxt,
                      "common":[p["simultaneous"] for p in pts]},indent=2))
    return 0


if __name__=="__main__":raise SystemExit(main())
