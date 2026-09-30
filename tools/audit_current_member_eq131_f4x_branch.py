#!/usr/bin/env python3
"""Test whether freeing f4X can satisfy Eq.131 while retaining the archived v6/f3 branch.

The current Eq.131 A2 ODE drives the fixed-f4/f4X branch into the new high-L
kinetic ghost.  The boundary value is not the cause: the current and source-f3
anchors have identical v6 at u~=0.61.

Here we hold the archived/raw v6 profile fixed (hence its derivative), keep f4
fixed, and solve Eq.131 pointwise for f4X.  We then reconstruct f3, f2F, f2,
f2X from the same JA/E00/E11 equations, complete the action jets, emit all
41 slots, and test K/G for every required L.

This is a diagnostic branch only.  It is not promoted.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ssz_p5.config import DEFAULT_L  # noqa: E402
from ssz_p5.numerics import module  # noqa: E402
from ssz_p5.production.electric_hybrid_onshell_central import _load_inputs, build_onshell_central  # noqa: E402
from ssz_p5.production.full_action_lower import complete_total_action_jets, emit_lower_slots  # noqa: E402
from ssz_p5.production.central_action import background_residuals  # noqa: E402

OUT = ROOT / "data/generated/spectral/CURRENT_MEMBER_EQ131_F4X_BRANCH_AUDIT.json"


def emit_action(d: pd.DataFrame) -> pd.DataFrame:
    zk = module("ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py")
    d = complete_total_action_jets(d.copy().sort_values("x").reset_index(drop=True))
    lower, d = emit_lower_slots(d)
    out = zk.emit(
        d,
        selected_v5=lower.v5.to_numpy(float),
        selected_c3=lower.c3.to_numpy(float),
        selected_e3=lower.e3.to_numpy(float),
        v6_phi_selector="action",
    )
    r = out.x.to_numpy(float)
    out["a5"] = (
        zk.dr(r, out.a2.to_numpy(float), 1, 9, 8)
        - zk.dr(r, out.a1.to_numpy(float), 2, 9, 8)
        - zk.dr(r, out.A0prime.to_numpy(float) * out.v4.to_numpy(float) / 2.0, 1, 9, 8)
        + out.A0prime.to_numpy(float) * out.v5.to_numpy(float) / 2.0
    )
    return out


def audit(stream: pd.DataFrame, red, L: int) -> dict:
    d = stream.sort_values("x").reset_index(drop=True)
    u = d.u.to_numpy(float)
    a = red.canonical_audit(d, int(L))
    K = np.asarray(a["K"], float)
    G = np.asarray(a["G"], float)
    Ks = (K + K.swapaxes(1,2))/2
    Gs = (G + G.swapaxes(1,2))/2
    ke = np.linalg.eigvalsh(Ks)[:,0]
    prod = (u>0.62)&(u<0.70)
    tail = (u>=0.70)&(u<0.71)
    cr = np.full(len(d), np.nan)
    for i in np.flatnonzero(ke>0):
        w,U=np.linalg.eigh(Ks[i])
        inv=U@np.diag(1/np.sqrt(w))@U.T
        C=inv@Gs[i]@inv
        cr[i]=np.linalg.eigvalsh((C+C.T)/2)[0]
    return {
        "production_min_K": float(np.min(ke[prod])),
        "tail_min_K": float(np.min(ke[tail])),
        "tail_negative_K_rows": int(np.sum(ke[tail]<=0)),
        "tail_min_cr2_where_defined": float(np.nanmin(cr[tail])) if np.any(np.isfinite(cr[tail])) else None,
        "tail_negative_cr2_rows_where_defined": int(np.sum(cr[tail] < 0)),
    }


def main() -> int:
    zk=module("ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py")
    red=module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    prof,raw,selected,_rank=_load_inputs(ROOT)

    r=raw.x.to_numpy(float); u=raw.u.to_numpy(float)
    f=raw.f.to_numpy(float); h=raw.h.to_numpy(float); ph=raw.phiprime.to_numpy(float)
    ap=raw.A0prime.to_numpy(float); X=raw.X.to_numpy(float)
    fp=zk.dr(r,f,1,9,8); fpp=zk.dr(r,f,2,9,8); hp=zk.dr(r,h,1,9,8); app=zk.dr(r,ap,1,9,8)
    v6=raw.v6.to_numpy(float); v6p=zk.dr(r,v6,1,9,8)
    f4=prof.f4.to_numpy(float); f4Xold=prof.N4.to_numpy(float)
    f3X=prof.f3X_integrated.to_numpy(float); f4XX=selected.f4XX_recovered.to_numpy(float)
    a4=raw.a4.to_numpy(float)

    ids=np.flatnonzero((u>=0.61)&(u<0.71))

    def vals(i:int, F4X:float):
        rr=r[i]; ff=f[i]; hh=h[i]; pp=ph[i]; aa=ap[i]
        FP=fp[i]; FPP=fpp[i]; HP=hp[i]; APP=app[i]; A4=a4[i]
        F4=f4[i]; V6=v6[i]
        C=2*hh**1.5*aa/(rr*np.sqrt(ff))
        F3=(V6/C+4*F4-hh*pp**2*F4X)/(rr*pp)
        F2F=-(4*rr*hh*pp*F3+8*(1-hh)*F4+2*hh**2*pp**2*F4X)/rr**2
        V10=-np.sqrt(ff*hh)/(2*rr)*(rr*F2F+2*hh*pp*F3+(hh*FP/ff)*(rr*pp*F3-4*F4+hh*pp**2*F4X))
        A7=(1-(4*hh*aa**2/ff)*F4)/(4*rr**2*np.sqrt(ff*hh))
        num=A4*(2*ff**2*(rr*HP+2*hh)+hh*rr**2*FP**2-ff*rr*(rr*FP*HP+2*hh*(rr*FPP+FP)))
        num-=ff*hh*rr**2*(aa*(4*V10*aa+V6*FP)+ff*V6*APP+8*A7*ff**2)
        rhs=num/(ff**2*hh*rr**2*aa)
        return F3,F2F,v6p[i]-rhs

    f4Xnew=np.full(len(raw),np.nan); f3=np.full(len(raw),np.nan); f2F=np.full(len(raw),np.nan)
    tiny=0
    for i in ids:
        q0=vals(i,0.0); q1=vals(i,1.0); slope=q1[-1]-q0[-1]
        if abs(slope)<1e-12:
            tiny+=1
            continue
        xx=-q0[-1]/slope
        q=vals(i,xx)
        f4Xnew[i]=xx; f3[i]=q[0]; f2F[i]=q[1]

    if tiny == len(ids):
        result = {
            "scope": "Eq131 f4X sensitivity audit; diagnostic only",
            "rows": int(len(ids)),
            "zero_sensitivity_rows": int(tiny),
            "diagnosis": "EQ131_IS_POINTWISE_INDEPENDENT_OF_F4X_ON_THIS_JA_REDUCED_BRANCH",
            "next_target": "TEST_F4_BRANCH_WITH_MODERN_FULL_ACTION_REEMISSION",
            "interpretation": (
                "With v6 fixed and f3/f2F eliminated through the exact Appendix-A/JA relations, "
                "the Eq131 residual has numerically zero linear sensitivity to f4X on every row. "
                "Therefore f4X cannot be used as the missing branch freedom for Eq131 here."
            ),
        }
        OUT.parent.mkdir(parents=True, exist_ok=True)
        OUT.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
        print(json.dumps(result, indent=2, allow_nan=False))
        return 0
    if tiny or not np.isfinite(f4Xnew[ids]).all():
        raise RuntimeError(f"partial f4X sensitivity failure on {tiny} rows")

    f2old=prof.f2.to_numpy(float); f2Xold=prof.f2X.to_numpy(float)
    A=r**2*f; kap=h*ph**2
    E00=r*f*hp-(f*(1-h)+r**2*(f*f2old-h*ap**2*f2F)-2*r*h**2*ph*ap**2*f3+h*ap**2*(4*(h-1)*f4-h**2*ph**2*f4Xnew))
    E11=r*h*fp-(f*(1-h)+r**2*(f*f2old+f*h*ph**2*f2Xold-h*ap**2*f2F)-2*r*h**2*ph*ap**2*(3*f3-h*ph**2*f3X)+h*ap**2*(4*(3*h-1)*f4-h*(9*h-4)*ph**2*f4Xnew+h**3*ph**4*f4XX))
    f2=f2old+E00/A
    f2X=f2Xold+(E11-E00)/(A*kap)

    base=build_onshell_central(ROOT).action.copy().sort_values("x").reset_index(drop=True)
    for target,delta in (("f2XX","principal_delta_f2XX"),("f2XF","principal_delta_f2XF"),("f2FF","principal_delta_f2FF")):
        base[target]=base[target].to_numpy(float)-base[delta].to_numpy(float)

    # base rows correspond to ids sorted in x; map source arrays by x interpolation.
    rb=base.x.to_numpy(float)
    def ip(arr): return np.interp(rb,r[np.argsort(r)],np.asarray(arr)[np.argsort(r)])
    base["f4X"]=ip(f4Xnew)
    base["f3"]=ip(f3); base["f2F"]=ip(f2F); base["f2"]=ip(f2); base["f2X"]=ip(f2X)

    # Refresh f2phi before total completion.
    fb=base.f.to_numpy(float); hb=base.h.to_numpy(float); ab=base.A0prime.to_numpy(float)
    Xb=base.X.to_numpy(float); phb=base.phiprime.to_numpy(float)
    Fbg=hb*ab**2/(2*fb); Ybg=4*Xb*Fbg
    Xp=zk.dr(rb,Xb,1,9,8); Fp=zk.dr(rb,Fbg,1,9,8); Yp=zk.dr(rb,Ybg,1,9,8)
    base["f2phi"]=(zk.dr(rb,base.f2.to_numpy(float),1,9,8)-base.f2X.to_numpy(float)*Xp-base.f2F.to_numpy(float)*Fp-base.f2Y.to_numpy(float)*Yp)/phb

    stream=emit_action(base)
    bg=background_residuals(complete_total_action_jets(base.copy()))
    ub=base.u.to_numpy(float); region=(ub>=0.61)&(ub<0.71); tail=(ub>=0.70)&(ub<0.71)
    perL={str(L):audit(stream,red,int(L)) for L in DEFAULT_L}

    dfx=base.f4X.to_numpy(float)-np.interp(rb,r[np.argsort(r)],f4Xold[np.argsort(r)])
    result={
        "scope":"pointwise Eq131 solve for f4X at fixed archived v6/f4; diagnostic only",
        "f4X_delta":{"median_abs":float(np.median(np.abs(dfx))),"max_abs":float(np.max(np.abs(dfx))),"min_new":float(np.min(base.f4X)),"max_new":float(np.max(base.f4X))},
        "background_max_abs_region":{k:float(np.max(np.abs(np.asarray(v)[region]))) for k,v in bg.items()},
        "eq131_max_abs_tail":float(np.max(np.abs(np.asarray(bg["eq131_residual"])[tail]))),
        "per_L":perL,
        "all_L_tail_K_pass":all(perL[str(L)]["tail_min_K"]>0 for L in DEFAULT_L),
        "high_L_tail_K_pass":all(perL[str(L)]["tail_min_K"]>0 for L in (110,420,1000)),
    }
    if result["all_L_tail_K_pass"]:
        diagnosis="F4X_BRANCH_SATISFIES_EQ131_AND_ALL_L_TAIL_K"
        nxt="TEST_FULL_RADIAL_G_AND_INTERFACE_PROVENANCE"
    elif result["high_L_tail_K_pass"]:
        diagnosis="F4X_BRANCH_REMOVES_NEW_HIGH_L_GHOST_BUT_OLD_LOW_L_GHOST_REMAINS"
        nxt="SEPARATE_LOW_L_CENTRAL_INSTABILITY_FROM_NEW_EQ131_BRANCH"
    else:
        diagnosis="F4X_BRANCH_DOES_NOT_REMOVE_NEW_HIGH_L_GHOST"
        nxt="TEST_EQ131_FORMULA_OR_OTHER_FIXED_ACTION_SECTORS"
    result["diagnosis"]=diagnosis; result["next_target"]=nxt
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps(result,indent=2,allow_nan=False))
    return 0

if __name__=="__main__":
    raise SystemExit(main())
