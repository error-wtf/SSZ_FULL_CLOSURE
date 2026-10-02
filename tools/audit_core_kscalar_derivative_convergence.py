#!/usr/bin/env python3
"""K_scalar derivative-convergence audit for the punctured quartic core.

The archived core table stores a positive K_scalar, but recomputing the same
formula from its own stored primitives with JET9D8 gives a negative pocket near
x~0.30.  This audit tests whether the discrepancy is numerical differentiation
or a provenance/formula mismatch.

It compares several local-polynomial derivative policies on the exact archived
primitive profiles, reconstructs Y' both numerically and from archived P1 when
available, and checks the stored algebraic identity K_scalar == 2 P1 - F.
"""
from __future__ import annotations
import json,sys
from pathlib import Path
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/"src"))
from ssz_p5.jets.jet9d8 import derivative

REF=ROOT/"data/authoritative/ssz_p5_F2_core_punctured_horndeski_unreduced_39of41_2026-09-14.csv"
OUT=ROOT/"data/generated/qnm_global_diagnostic/CORE_KSCALAR_DERIVATIVE_CONVERGENCE.json"
POLICIES=((5,4),(7,6),(9,8),(11,8),(13,10),(15,12))

def main():
 d=pd.read_csv(REF).sort_values("x").reset_index(drop=True)
 d=d[(d.x>=0.02)&(d.u>=0.715)].reset_index(drop=True)
 r=d.x.to_numpy(float);f=d.f.to_numpy(float);h=d.h.to_numpy(float)
 ph=(d.phi_r if "phi_r" in d else d.phiprime).to_numpy(float)
 a1=d.a1.to_numpy(float);F=d.F_tensor.to_numpy(float);H=d.H_tensor.to_numpy(float)
 a4=.5*np.sqrt(f*h)*H
 mu=2*(ph*a1+2*r*a4)/np.sqrt(f*h)
 Y=f*r**4*H**4/(mu**2*h)
 pref=h*mu/(2*f*r**2*H**2)
 storedK=d.K_scalar.to_numpy(float)
 rep={
   "stored":{"min_K":float(storedK.min()),"negative_rows":int(np.sum(storedK<=0))},
   "stored_identity":{},
   "policies":{}
 }
 if "P1" in d:
  P1=d.P1.to_numpy(float)
  resid=storedK-(2*P1-F)
  rep["stored_identity"]={
    "P1_present":True,
    "max_abs_K_minus_2P1plusF":float(np.max(np.abs(resid))),
    "max_scaled":float(np.max(np.abs(resid)/np.maximum(1,np.abs(storedK))))
  }
  targetYp=np.divide(P1,pref,out=np.full_like(P1,np.nan),where=np.abs(pref)>1e-300)
 else:
  P1=None;targetYp=None;rep["stored_identity"]={"P1_present":False}

 for w,deg in POLICIES:
  if w>len(r):continue
  yp=derivative(r,Y,1,window=w,degree=deg)
  p1=pref*yp;k=2*p1-F;j=int(np.argmin(k))
  z={
   "window":w,"degree":deg,"min_K":float(k[j]),"x_min":float(r[j]),"u_min":float(d.u.iloc[j]),
   "negative_rows":int(np.sum(k<=0)),
   "max_abs_K_vs_stored":float(np.max(np.abs(k-storedK))),
   "p99_abs_K_vs_stored":float(np.quantile(np.abs(k-storedK),.99))
  }
  if P1 is not None:
   m=np.isfinite(targetYp)
   z["Yprime_vs_storedP1"]={
     "max_abs":float(np.max(np.abs(yp[m]-targetYp[m]))),
     "p99_abs":float(np.quantile(np.abs(yp[m]-targetYp[m]),.99)),
     "relative_L2":float(np.linalg.norm(yp[m]-targetYp[m])/max(np.linalg.norm(targetYp[m]),1e-300))
   }
  rep["policies"][f"{w}x{deg}"]=z

 # Is the sign failure stable across all derivative policies?
 mins=[z["min_K"] for z in rep["policies"].values()]
 if all(v<0 for v in mins):
  diag="NEGATIVE_POCKET_STENCIL_ROBUST"
 elif all(v>0 for v in mins):
  diag="POSITIVE_ACROSS_STENCILS"
 else:
  diag="SIGN_IS_DERIVATIVE_POLICY_SENSITIVE"
 rep["status"]="CORE_KSCALAR_DERIVATIVE_CONVERGENCE_COMPLETE"
 rep["diagnosis"]=diag
 rep["next_target"]="TRACE_ARCHIVED_P1_PROVENANCE_OR_ANALYTIC_ACTION_DERIVATIVE" if diag=="NEGATIVE_POCKET_STENCIL_ROBUST" else "RESOLVE_DERIVATIVE_POLICY_BEFORE_CORE_SPECTRAL_CLAIM"
 OUT.parent.mkdir(parents=True,exist_ok=True);OUT.write_text(json.dumps(rep,indent=2,allow_nan=False)+"\n")
 print(json.dumps(rep,indent=2,allow_nan=False));return 0
if __name__=="__main__":raise SystemExit(main())
