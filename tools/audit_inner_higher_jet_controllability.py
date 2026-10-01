#!/usr/bin/env python3
"""Inner-tail higher-action-jet controllability audit.

After the background-null f2 Hessian basis stalled at finite amplitude, test the
next legitimate same-background freedoms: f3XX, f4XX and f4XXX.  These are
higher action jets; f3,f3X,f4,f4X and the background fields remain fixed.

For localized bumps centered in the native inner tail, this audit:
  * re-completes all mixed phi jets from the total-action chain rules;
  * re-emits the full Appendix-A 41-slot stream;
  * verifies E00,E11,JA,Eq131 remain unchanged to numerical precision;
  * computes projected first-order gains of the weakest radial c_r^2 for all L;
  * searches for one common direction that increases c_r^2 for every required L
    while not decreasing K.

No finite repair or promotion is performed.
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

OUT=ROOT/"data/generated/spectral/INNER_HIGHER_JET_PROJECTED_CONTROLLABILITY.json"
JETS=("f3XX","f4XX","f4XXX")
POINTS=(0.6190,0.6160,0.6130,0.6106,0.6100)
WIDTH=6.5e-4
EPS=2e-5
LS=tuple(int(x) for x in DEFAULT_L)


def bump(u,c):
    z=(np.asarray(u,float)-c)/WIDTH
    q=np.zeros_like(z);m=np.abs(z)<1
    q[m]=np.exp(-1/(1-z[m]**2))/np.exp(-1)
    return np.where(u<=0.62,q,0.0)


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


def metrics(stream,red,L,i):
    a=red.canonical_audit(stream,int(L))
    K=np.asarray(a["K"],float);G=np.asarray(a["G"],float)
    Ks=(K[i]+K[i].T)/2;Gs=(G[i]+G[i].T)/2
    w,U=np.linalg.eigh(Ks);km=float(w[0])
    if km<=0:return km,np.nan
    inv=U@np.diag(1/np.sqrt(w))@U.T
    C=inv@Gs@inv
    return km,float(np.linalg.eigvalsh((C+C.T)/2)[0])


def bgvec(comp,i):
    b=background_residuals(comp)
    return np.array([b["E00"][i],b["E11"][i],b["JA"][i],b["eq131_residual"][i]],float)


def point(base,u0):
    loc=base[(base.u>=u0-0.0022)&(base.u<=u0+0.0022)].copy().sort_values("x").reset_index(drop=True)
    u=loc.u.to_numpy(float);i=int(np.argmin(np.abs(u-u0)))
    comp0,st0=reemit(loc);E0=bgvec(comp0,i)
    red=module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    base_m={L:metrics(st0,red,L,i) for L in LS}
    dirs=[]
    for c in POINTS:
        if abs(c-u0)<=0.0025:
            for name in JETS:dirs.append((c,name))
    Gcr={L:np.zeros(len(dirs)) for L in LS}
    Gk={L:np.zeros(len(dirs)) for L in LS}
    B=np.zeros((4,len(dirs)))
    for j,(c,name) in enumerate(dirs):
        d=loc.copy();d[name]=d[name].to_numpy(float)+EPS*bump(u,c)
        comp,st=reemit(d)
        B[:,j]=(bgvec(comp,i)-E0)/EPS
        for L in LS:
            km,cr=metrics(st,red,L,i)
            Gk[L][j]=(km-base_m[L][0])/EPS
            Gcr[L][j]=(cr-base_m[L][1])/EPS

    # Treat background changes as equalities, though expected nearly zero.
    rn=np.linalg.norm(B,axis=1);keep=rn>1e-12
    Be=B[keep]/rn[keep,None] if np.any(keep) else np.zeros((0,len(dirs)))

    cobj=np.zeros(len(dirs)+1);cobj[-1]=-1
    Aub=[];bub=[]
    for L in LS:
        row=np.zeros(len(dirs)+1);row[:-1]=-Gcr[L];row[-1]=1
        Aub.append(row);bub.append(0.)
        row=np.zeros(len(dirs)+1);row[:-1]=-Gk[L]
        Aub.append(row);bub.append(0.)
    Aeq=np.column_stack([Be,np.zeros(len(Be))]) if len(Be) else None
    sol=linprog(cobj,A_ub=np.asarray(Aub),b_ub=np.asarray(bub),
                A_eq=Aeq,b_eq=np.zeros(len(Be)) if len(Be) else None,
                bounds=[(-1,1)]*len(dirs)+[(None,None)],method="highs")
    lp={"success":False}
    if sol.success:
        dp=sol.x[:-1]
        lp={
          "success":True,
          "common_cr2_gain":float(sol.x[-1]),
          "all_positive":bool(sol.x[-1]>0),
          "cr2_gains":{str(L):float(Gcr[L]@dp) for L in LS},
          "K_gains":{str(L):float(Gk[L]@dp) for L in LS},
          "direction":{f"{dirs[j][1]}@{dirs[j][0]:.4f}":float(dp[j]) for j in range(len(dp)) if abs(dp[j])>1e-8},
        }
    return {
      "u_requested":u0,"u_grid":float(u[i]),
      "directions":[f"{n}@{c:.4f}" for c,n in dirs],
      "background_base":dict(zip(("E00","E11","JA","Eq131"),map(float,E0))),
      "background_jacobian_max_abs":float(np.max(np.abs(B))) if B.size else 0.0,
      "base_per_L":{str(L):{"K":base_m[L][0],"cr2":base_m[L][1]} for L in LS},
      "simultaneous_all_L":lp,
    }


def main():
    base=build_onshell_central(ROOT).action.sort_values("x").reset_index(drop=True)
    pts=[point(base,p) for p in POINTS]
    all_common=all(p["simultaneous_all_L"].get("all_positive",False) for p in pts)
    if all_common:
        status="INNER_HIGHER_JETS_COMMON_TANGENT_EXISTS"
        nxt="BUILD_FINITE_HIGHER_JET_INNER_CONTINUATION"
    elif any(p["simultaneous_all_L"].get("all_positive",False) for p in pts):
        status="INNER_HIGHER_JETS_PARTIAL_CONTROLLABILITY"
        nxt="COMBINE_HIGHER_JETS_WITH_F2_HESSIAN_IN_MULTIPOINT_SOLVER"
    else:
        status="INNER_HIGHER_JETS_NO_COMMON_TANGENT"
        nxt="EXPAND_TO_ADDITIONAL_ACTION_SECTORS_OR_INTERFACE_BRANCH"
    out={"status":status,"next_target":nxt,"points":pts,
         "scope":"f3XX/f4XX/f4XXX same-background higher-jet controllability; no finite promotion"}
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
    print(json.dumps(out,indent=2,allow_nan=False))
    return 0

if __name__=="__main__":raise SystemExit(main())
