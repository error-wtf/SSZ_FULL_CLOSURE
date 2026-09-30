#!/usr/bin/env python3
"""Eq.131 boundary-condition branch scan for the Electric A2 solve.

The current A2 construction solves the first-order Eq.131 ODE for v6 with one
boundary value at u~=0.61 taken from the archived raw stream.  The previous
audit showed that enforcing Eq.131 creates the new high-L tail ghost, whereas
a source-f3 algebraic branch keeps L>=110 healthy but violates Eq.131.

This audit varies ONLY the Eq.131 boundary value.  The alternative endpoint is
not fitted: it is the v6 implied at the same anchor by the archived/source f3
through the exact Appendix-A relation v6=4 alpha2.  We scan a fixed interpolation
parameter lambda between the current boundary and that source-f3 boundary,
re-solve Eq.131, reconstruct f3/f2F/f2/f2X from the same equations, re-emit the
full action, and test K for all required L.

Diagnostic only: no selected lambda is promoted as a production member.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ssz_p5.config import DEFAULT_L  # noqa: E402
from ssz_p5.numerics import module  # noqa: E402
from ssz_p5.production.central_action import central_action_inputs, background_residuals  # noqa: E402
from ssz_p5.production.electric_hybrid_onshell_central import (  # noqa: E402
    _load_inputs,
    build_onshell_central,
)
from ssz_p5.production.full_action_lower import complete_total_action_jets, emit_lower_slots  # noqa: E402

OUT = ROOT / "data/generated/spectral/CURRENT_MEMBER_EQ131_BOUNDARY_BRANCH_SCAN.json"
LAMBDAS = (-1.0, -0.5, 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0)


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


def kmatrix_metrics(stream: pd.DataFrame, reducer, L: int) -> dict:
    d = stream.sort_values("x").reset_index(drop=True)
    u = d.u.to_numpy(float)
    a = reducer.canonical_audit(d, int(L))
    K = np.asarray(a["K"], float)
    Ks = (K + K.swapaxes(1, 2)) / 2
    k = np.linalg.eigvalsh(Ks)[:, 0]
    prod = (u > 0.62) & (u < 0.70)
    tail = (u >= 0.70) & (u < 0.71)
    return {
        "production_min_K": float(np.min(k[prod])),
        "tail_min_K": float(np.min(k[tail])),
        "tail_negative_rows": int(np.sum(k[tail] <= 0)),
    }


def main() -> int:
    zk = module("ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py")
    reducer = module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    prof, raw, selected, _rank = _load_inputs(ROOT)

    r = raw.x.to_numpy(float)
    u = raw.u.to_numpy(float)
    f = raw.f.to_numpy(float)
    h = raw.h.to_numpy(float)
    ph = raw.phiprime.to_numpy(float)
    ap = raw.A0prime.to_numpy(float)
    X = raw.X.to_numpy(float)
    fp = zk.dr(r, f, 1, 9, 8)
    fpp = zk.dr(r, f, 2, 9, 8)
    hp = zk.dr(r, h, 1, 9, 8)
    app = zk.dr(r, ap, 1, 9, 8)

    f4 = prof.f4.to_numpy(float)
    f4x = prof.N4.to_numpy(float)
    f3x = prof.f3X_integrated.to_numpy(float)
    f4xx = selected.f4XX_recovered.to_numpy(float)
    a4 = raw.a4.to_numpy(float)

    band = (u >= 0.61) & (u <= 0.715)
    idx = np.flatnonzero(band)
    order = idx[np.argsort(r[idx])]
    rs = r[order]
    fields = {
        name: PchipInterpolator(rs, arr[order], extrapolate=False)
        for name, arr in (
            ("f", f), ("h", h), ("ph", ph), ("ap", ap),
            ("fp", fp), ("fpp", fpp), ("hp", hp), ("app", app),
            ("f4", f4), ("f4x", f4x), ("a4", a4),
        )
    }

    def aux(rr: float, vv: float):
        ff=float(fields["f"](rr)); hh=float(fields["h"](rr)); pp=float(fields["ph"](rr))
        aa=float(fields["ap"](rr)); F4=float(fields["f4"](rr)); F4x=float(fields["f4x"](rr))
        c = 2.0*hh**1.5*aa/(rr*np.sqrt(ff))
        F3 = (vv/c + 4.0*F4 - hh*pp**2*F4x)/(rr*pp)
        F2F = -(4.0*rr*hh*pp*F3 + 8.0*(1.0-hh)*F4 + 2.0*hh**2*pp**2*F4x)/rr**2
        FP=float(fields["fp"](rr))
        v10=-np.sqrt(ff*hh)/(2.0*rr)*(rr*F2F + 2.0*hh*pp*F3 + (hh*FP/ff)*(rr*pp*F3 - 4.0*F4 + hh*pp**2*F4x))
        a7=(1.0-(4.0*hh*aa**2/ff)*F4)/(4.0*rr**2*np.sqrt(ff*hh))
        return F3,F2F,v10,a7

    def rhs(rr, y):
        vv=float(y[0]); ff=float(fields["f"](rr)); hh=float(fields["h"](rr)); aa=float(fields["ap"](rr))
        FP=float(fields["fp"](rr)); FPP=float(fields["fpp"](rr)); HP=float(fields["hp"](rr)); APP=float(fields["app"](rr))
        A4=float(fields["a4"](rr)); _,_,v10,a7=aux(rr,vv)
        num=A4*(2*ff**2*(rr*HP+2*hh)+hh*rr**2*FP**2-ff*rr*(rr*FP*HP+2*hh*(rr*FPP+FP)))
        num-=ff*hh*rr**2*(aa*(4*v10*aa+vv*FP)+ff*vv*APP+8*a7*ff**2)
        den=ff**2*hh*rr**2*aa
        return [num/den]

    i0=int(np.argmin(np.abs(u-0.61)))
    r0=float(r[i0]); r1=float(rs.min())
    v0_current=float(raw.v6.iloc[i0])

    source=central_action_inputs(ROOT).sort_values("x")
    f3_source_anchor=float(np.interp(r0,source.x.to_numpy(float),source.f3.to_numpy(float)))
    ff=f[i0]; hh=h[i0]; pp=ph[i0]; aa=ap[i0]; rr=r[i0]; F4=f4[i0]; F4x=f4x[i0]
    c=2.0*hh**1.5*aa/(rr*np.sqrt(ff))
    v0_source_f3=c*(rr*pp*f3_source_anchor - 4.0*F4 + hh*pp**2*F4x)

    # Start from current pre-Hessian action to retain all fixed action jets.
    base=build_onshell_central(ROOT).action.copy().sort_values("x").reset_index(drop=True)
    for target,delta in (("f2XX","principal_delta_f2XX"),("f2XF","principal_delta_f2XF"),("f2FF","principal_delta_f2FF")):
        base[target]=base[target].to_numpy(float)-base[delta].to_numpy(float)

    ridx=np.flatnonzero((u>=0.61)&(u<0.71))
    ridx=ridx[np.argsort(r[ridx])]
    rrgrid=r[ridx]

    scans={}
    for lam in LAMBDAS:
        v0=v0_current + lam*(v0_source_f3-v0_current)
        sol=solve_ivp(rhs,(r0,r1),[v0],rtol=2e-11,atol=2e-12,dense_output=True,max_step=2e-4)
        if not sol.success:
            scans[str(lam)]={"solve_success":False,"message":sol.message}
            continue

        v6=np.full(len(raw),np.nan); v6[idx]=sol.sol(r[idx])[0]
        f3=np.full(len(raw),np.nan); f2F=np.full(len(raw),np.nan)
        for j in idx:
            f3[j],f2F[j],_,_=aux(float(r[j]),float(v6[j]))

        f2old=prof.f2.to_numpy(float); f2xold=prof.f2X.to_numpy(float)
        A=r**2*f; kappa=h*ph**2
        E00=r*f*hp-(f*(1-h)+r**2*(f*f2old-h*ap**2*f2F)-2*r*h**2*ph*ap**2*f3+h*ap**2*(4*(h-1)*f4-h**2*ph**2*f4x))
        E11=r*h*fp-(f*(1-h)+r**2*(f*f2old+f*h*ph**2*f2xold-h*ap**2*f2F)-2*r*h**2*ph*ap**2*(3*f3-h*ph**2*f3x)+h*ap**2*(4*(3*h-1)*f4-h*(9*h-4)*ph**2*f4x+h**3*ph**4*f4xx))
        f2new=f2old+E00/A
        f2xnew=f2xold+(E11-E00)/(A*kappa)

        d=base.copy()
        # base is already the same ridx sorted by x.
        d["f3"]=f3[ridx]; d["f2F"]=f2F[ridx]; d["f2"]=f2new[ridx]; d["f2X"]=f2xnew[ridx]; d["v6_A2_resolved"]=v6[ridx]

        Fbg=d.h.to_numpy(float)*d.A0prime.to_numpy(float)**2/(2*d.f.to_numpy(float))
        Xp=zk.dr(rrgrid,d.X.to_numpy(float),1,9,8); Fp=zk.dr(rrgrid,Fbg,1,9,8)
        d["f2phi"]=(zk.dr(rrgrid,d.f2.to_numpy(float),1,9,8)-d.f2X.to_numpy(float)*Xp-d.f2F.to_numpy(float)*Fp)/d.phiprime.to_numpy(float)

        stream=emit_action(d)
        bg=background_residuals(complete_total_action_jets(d.copy()))
        du=d.u.to_numpy(float); region=(du>=0.61)&(du<0.71); tail=(du>=0.70)&(du<0.71)
        perL={str(L):kmatrix_metrics(stream,reducer,int(L)) for L in DEFAULT_L}
        scans[str(lam)]={
            "solve_success":True,
            "v0":float(v0),
            "per_L":perL,
            "all_L_tail_K_pass":all(perL[str(L)]["tail_min_K"]>0 for L in DEFAULT_L),
            "high_L_tail_K_pass":all(perL[str(L)]["tail_min_K"]>0 for L in (110,420,1000)),
            "eq131_max_abs_region":float(np.max(np.abs(np.asarray(bg["eq131_residual"])[region]))),
            "eq131_max_abs_tail":float(np.max(np.abs(np.asarray(bg["eq131_residual"])[tail]))),
            "E00_max_abs_region":float(np.max(np.abs(np.asarray(bg["E00"])[region]))),
            "E11_max_abs_region":float(np.max(np.abs(np.asarray(bg["E11"])[region]))),
            "JA_max_abs_region":float(np.max(np.abs(np.asarray(bg["JA"])[region]))),
        }

    healthy=[float(k) for k,v in scans.items() if v.get("solve_success") and v.get("high_L_tail_K_pass")]
    current=scans["0.0"]
    source_boundary=scans["1.0"]
    if healthy:
        diagnosis="EQ131_HAS_BOUNDARY_BRANCHES_WITH_HEALTHY_NEW_HIGH_L_TAIL"
        next_target="CHECK_ACTION_PROVENANCE_INTERFACE_MATCHING_FOR_HEALTHY_EQ131_BRANCHES"
    else:
        diagnosis="EQ131_BOUNDARY_SCAN_DID_NOT_FIND_HEALTHY_NEW_HIGH_L_TAIL"
        next_target="TEST_EQ131_FORMULA_BRANCH_OR_FIXED_F4_SECTOR"
    payload={
        "scope":"fixed boundary-condition branch scan; diagnostic only; no optimization/no promotion",
        "v0_current_raw":v0_current,
        "v0_source_f3_implied":float(v0_source_f3),
        "lambda_definition":"v0(lambda)=v0_current + lambda*(v0_source_f3_implied-v0_current)",
        "lambdas":list(LAMBDAS),
        "scans":scans,
        "healthy_highL_lambdas":healthy,
        "current_lambda0":current,
        "source_boundary_lambda1":source_boundary,
        "diagnosis":diagnosis,
        "next_target":next_target,
        "guard":"A healthy diagnostic branch is not a valid member until same-action provenance, interface matching, all-L K/G, and downstream gates pass.",
    }
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(payload,indent=2,allow_nan=False)+"\n")
    print(json.dumps({
        "diagnosis":diagnosis,
        "next_target":next_target,
        "v0_current_raw":v0_current,
        "v0_source_f3_implied":float(v0_source_f3),
        "healthy_highL_lambdas":healthy,
        "lambda0_highL":current.get("high_L_tail_K_pass"),
        "lambda1_highL":source_boundary.get("high_L_tail_K_pass"),
        "lambda0_L1000":current.get("per_L",{}).get("1000"),
        "lambda1_L1000":source_boundary.get("per_L",{}).get("1000"),
    },indent=2,allow_nan=False))
    return 0

if __name__=="__main__":
    raise SystemExit(main())
