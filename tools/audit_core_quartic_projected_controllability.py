#!/usr/bin/env python3
"""Projected controllability of the quartic core kinetic ghost with legitimate
background-null Horndeski action controls.

Base state:
  direct quartic G5=0 primitive emission on 0.715<=u<50 using the action-derived
  F,G,H,a1,c4 profiles and the documented selected c2 diagnostic primitive.

Controls (all are background-null by construction):
  1) Delta G3 = 1/2 q(phi)[X-X_b(phi)]^2
     -> exact c2,c3,e3 response from horndeski_tubular.g3_tubular_response_matrices.
  2) Delta G4 = 1/2 q(phi)[X-X_b(phi)]^2 plus the compensating Delta G3,
     Delta G2, Delta G2X section from g4xx_background_null_section.
     Its induced a1,c4 response is re-emitted through mh_general_primitives.
     Net c2 is held fixed, corresponding to the documented independent G2XX
     normal cancellation wherever dc2/dG2XX is nonzero.

This audit asks a narrow falsifiable question:
Does the legitimate background-null action tangent contain a common direction
that raises the weakest K eigenvalue for all required L at the observed core
failure region around u~3.322844?

A positive answer is only a LOCAL first-order controllability result and licenses
a finite predictor-corrector.  It is not a repaired core member.
"""
from __future__ import annotations
import json,sys
from pathlib import Path
import numpy as np,pandas as pd
from scipy.optimize import linprog

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))

from ssz_p5.coefficients.mh_action_primitives import quartic_g5zero_primitives
from ssz_p5.coefficients.mh_general_primitives import emit_from_primitives
from ssz_p5.production.horndeski_tubular import (
    g3_tubular_response_matrices,g4xx_background_null_section
)
from ssz_p5.numerics import module

ACT=ROOT/"data/production/ssz_p5_F1b_G4XX_transverse_core_candidate_2026-09-14.csv"
REF=ROOT/"data/authoritative/ssz_p5_F2_core_punctured_horndeski_unreduced_39of41_2026-09-14.csv"
OUT=ROOT/"data/generated/qnm_global_diagnostic/CORE_QUARTIC_PROJECTED_CONTROLLABILITY.json"

LS=(6,12,20,42,110,420,1000)
CENTERS=(2.45,2.85,3.15,3.322844323058915,3.50,3.80,4.20)
WIDTH=0.34
EPS_G3=1e-3
EPS_G4=1e-5
TARGETS=(2.85,3.15,3.322844323058915,3.50,3.80)


def bump(u,c):
    z=(np.asarray(u,float)-c)/WIDTH
    q=np.zeros_like(z)
    m=np.abs(z)<1
    q[m]=np.exp(-1/(1-z[m]**2))/np.exp(-1)
    return q


def build_base():
    act=pd.read_csv(ACT).sort_values("r_over_rs").reset_index(drop=True)
    ref=pd.read_csv(REF).sort_values("x").reset_index(drop=True)
    q=quartic_g5zero_primitives(act)
    if len(q)!=len(ref) or np.max(np.abs(q.x.to_numpy(float)-ref.x.to_numpy(float)))>1e-12:
        raise RuntimeError("action/reference grid mismatch")
    inp=pd.DataFrame({
      "u":act.u.to_numpy(float),
      "x":act.r_over_rs.to_numpy(float),
      "phi":act.phi.to_numpy(float),
      "f":act.A_f.to_numpy(float),
      "h":act.B_h.to_numpy(float),
      "X":act.X.to_numpy(float),
      "a1":q.a1_action.to_numpy(float),
      "c2":ref.c2.to_numpy(float),
      "c4":q.c4_action.to_numpy(float),
      "F_tensor":q.F_tensor_action.to_numpy(float),
      "G_tensor":q.G_tensor_action.to_numpy(float),
      "H_tensor":q.H_tensor_action.to_numpy(float),
    })
    emitted=emit_from_primitives(inp,regularize_photon_root=True)
    mask=(emitted.u.to_numpy(float)>=0.715)&(emitted.u.to_numpy(float)<50.0)
    inp=inp.loc[mask].sort_values("x").reset_index(drop=True)
    emitted=emitted.loc[mask].sort_values("x").reset_index(drop=True)
    return inp,emitted


def Kmins_for_L(stream,red,L,indices):
    a=red.canonical_audit(stream,float(L))
    K=np.asarray(a["K"],float)
    vals=np.linalg.eigvalsh((K+K.transpose(0,2,1))/2)[:,0]
    return {int(i):float(vals[int(i)]) for i in indices}


def apply_g3(base,resp,q):
    d=base.copy()
    for slot in ("c2","c3","e3"):
        delta=np.asarray(resp[slot]@q).ravel()
        d[slot]=d[slot].to_numpy(float)+delta
    return d


def apply_g4(inp,base,qprof):
    sec=g4xx_background_null_section(inp,qprof)
    changed=inp.copy()
    changed["a1"]=changed.a1.to_numpy(float)+sec.delta_a1.to_numpy(float)
    changed["c4"]=changed.c4.to_numpy(float)+sec.delta_c4.to_numpy(float)
    # c2 held fixed: use independent background-null G2XX normal cancellation.
    out=emit_from_primitives(changed,regularize_photon_root=True)
    return out,sec


def main():
    inp,base=build_base()
    u=base.u.to_numpy(float)
    red=module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    resp=g3_tubular_response_matrices(inp)
    target_idx=[int(np.argmin(np.abs(u-t))) for t in TARGETS]

    controls=[]
    for c in CENTERS:
        controls.append(("G3",c))
        controls.append(("G4XX",c))

    # Baseline values.
    baseK={}
    for L in LS:
        kk=Kmins_for_L(base,red,L,target_idx)
        for i,v in kk.items(): baseK[(L,i)]=v

    # Jacobian of weakest K wrt legitimate control amplitudes.
    J={(L,i):np.zeros(len(controls)) for L in LS for i in target_idx}
    bg_checks=[]
    dc2_coeff=quartic_g5zero_primitives(pd.read_csv(ACT).sort_values("r_over_rs").reset_index(drop=True))
    # align mask to quartic domain for G2XX cancellation conditioning
    qfull=dc2_coeff.dc2_dG2XX_action.to_numpy(float)
    ufull=pd.read_csv(ACT).sort_values("r_over_rs").reset_index(drop=True).u.to_numpy(float)
    mask=(ufull>=0.715)&(ufull<50.0)
    dc2=qfull[mask]

    for j,(kind,c) in enumerate(controls):
        shape=bump(u,c)
        if kind=="G3":
            pert=apply_g3(base,resp,EPS_G3*shape)
            nullmax=0.0
            eps=EPS_G3
            conditioning=None
        else:
            # Explicitly refuse support where the c2-cancelling G2XX normal is singular.
            support=np.abs(shape)>1e-10
            if np.any(support & (np.abs(dc2)<1e-10)):
                bg_checks.append({
                  "control":f"{kind}@{c}","usable":False,
                  "reason":"dc2/dG2XX singular on support"
                })
                for L in LS:
                    for i in target_idx:J[(L,i)][j]=0.0
                continue
            pert,sec=apply_g4(inp,base,EPS_G4*shape)
            nullmax=float(max(np.max(np.abs(sec.delta_E00)),np.max(np.abs(sec.delta_E11)),np.max(np.abs(sec.delta_E22))))
            eps=EPS_G4
            conditioning=float(np.min(np.abs(dc2[support]))) if np.any(support) else None
        bg_checks.append({
          "control":f"{kind}@{c}","usable":True,
          "background_null_max_abs":nullmax,
          "min_abs_dc2_dG2XX_on_support":conditioning
        })
        for L in LS:
            kk=Kmins_for_L(pert,red,L,target_idx)
            for i in target_idx:
                J[(L,i)][j]=(kk[int(i)]-baseK[(L,i)])/eps

    usable=np.array([b.get("usable",False) for b in bg_checks],bool)
    # Note bg_checks order follows controls one-for-one.
    if len(usable)!=len(controls):
        raise RuntimeError("control bookkeeping mismatch")

    # LP: maximize common first-order K gain t over all (L,target), bounded controls.
    n=len(controls)
    obj=np.zeros(n+1);obj[-1]=-1
    Aub=[];bub=[]
    for L in LS:
        for i in target_idx:
            row=np.zeros(n+1);row[:n]=-J[(L,i)];row[-1]=1
            Aub.append(row);bub.append(0.0)
    bounds=[(-1,1) if usable[j] else (0,0) for j in range(n)]+[(None,None)]
    sol=linprog(obj,A_ub=np.asarray(Aub),b_ub=np.asarray(bub),bounds=bounds,method="highs")

    point_results=[]
    if sol.success:
        p=sol.x[:-1]
        for i in target_idx:
            point_results.append({
              "u":float(u[i]),
              "baseline":{str(L):baseK[(L,i)] for L in LS},
              "linear_gain":{str(L):float(J[(L,i)]@p) for L in LS},
              "linear_predicted":{str(L):float(baseK[(L,i)]+J[(L,i)]@p) for L in LS},
            })
        direction={
          f"{controls[j][0]}@{controls[j][1]:.6f}":float(p[j])
          for j in range(n) if abs(p[j])>1e-8
        }
        gain=float(sol.x[-1])
    else:
        direction={};gain=None

    report={
      "status":(
        "CORE_QUARTIC_COMMON_K_IMPROVING_TANGENT_EXISTS"
        if sol.success and gain is not None and gain>0 else
        "CORE_QUARTIC_NO_COMMON_K_IMPROVING_TANGENT_FOUND"
      ),
      "scope":"local first-order action-level controllability around u~3.322844; no promotion",
      "centers":list(CENTERS),"width":WIDTH,"targets":list(TARGETS),
      "controls":bg_checks,
      "lp_success":bool(sol.success),
      "common_gain":gain,
      "direction":direction,
      "points":point_results,
      "next_target":(
        "BUILD_FINITE_CORE_ACTION_PREDICTOR_CORRECTOR"
        if sol.success and gain is not None and gain>0 else
        "EXPAND_CORE_ACTION_CONTROL_BASIS_OR_DERIVE_FULL_C2_ACTION_BRIDGE"
      ),
      "guard":"Base c2 remains the selected diagnostic primitive. Positive tangent proves only that legitimate background-null Horndeski normal directions can improve its finite-L kinetic obstruction locally."
    }
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps(report,indent=2,allow_nan=False))
    return 0

if __name__=="__main__":raise SystemExit(main())
