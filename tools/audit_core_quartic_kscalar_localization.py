#!/usr/bin/env python3
"""Localize the quartic-core K_scalar discrepancy.

Fresh action-derived quartic primitives match the archived core primitives very
closely pointwise, yet recomputed K_scalar becomes negative.  Since K_scalar
contains a radial derivative of Y=f r^4 H^4/(mu^2 h), tiny primitive differences
can be derivative-amplified in the deep core.

This audit separates:
  * action-vs-archive primitive value mismatch,
  * archive K_scalar column vs K_scalar recomputed from archived primitives,
  * action K_scalar vs recomputed archive K_scalar,
  * one-at-a-time substitutions F,H,a1 to identify the sensitive carrier.

No coefficient tuning is performed.
"""
from __future__ import annotations
import json,sys
from pathlib import Path
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/"src"))
from ssz_p5.coefficients.mh_action_primitives import quartic_g5zero_primitives
from ssz_p5.coefficients.mh_general_primitives import emit_from_primitives

ACT=ROOT/"data/production/ssz_p5_F1b_G4XX_transverse_core_candidate_2026-09-14.csv"
REF=ROOT/"data/authoritative/ssz_p5_F2_core_punctured_horndeski_unreduced_39of41_2026-09-14.csv"
OUT=ROOT/"data/generated/qnm_global_diagnostic/CORE_QUARTIC_KSCALAR_LOCALIZATION.json"

def main():
 a=pd.read_csv(ACT).sort_values("r_over_rs").reset_index(drop=True)
 r=pd.read_csv(REF).sort_values("x").reset_index(drop=True)
 q=quartic_g5zero_primitives(a)
 if len(q)!=len(r) or np.max(np.abs(q.x-r.x))>1e-12: raise RuntimeError("grid mismatch")
 x=q.x.to_numpy(float);mask=(x>=0.02)&(q.u.to_numpy(float)>=0.715)
 q=q.loc[mask].reset_index(drop=True);a=a.loc[mask].reset_index(drop=True);r=r.loc[mask].reset_index(drop=True)
 c2=r.c2.to_numpy(float)

 def inp(F,H,a1,G=None,c4=None):
  return pd.DataFrame({
   "x":q.x,"u":q.u,"f":a.A_f,"h":a.B_h,"X":a.X,"phi":a.phi,"A0prime":0.0,
   "a1":a1,"c2":c2,"c4":q.c4_action if c4 is None else c4,
   "F_tensor":F,"G_tensor":q.G_tensor_action if G is None else G,"H_tensor":H
  })
 def calc(F,H,a1,G=None,c4=None):
  z=emit_from_primitives(inp(F,H,a1,G,c4),regularize_photon_root=True)
  k=z.K_scalar.to_numpy(float)
  j=int(np.argmin(k))
  return {"min":float(k[j]),"x":float(z.x.iloc[j]),"u":float(z.u.iloc[j]),"negative_rows":int(np.sum(k<=0)),"array":k}

 Fa=q.F_tensor_action.to_numpy(float);Ha=q.H_tensor_action.to_numpy(float);Aa=q.a1_action.to_numpy(float)
 Fr=r.F_tensor.to_numpy(float);Hr=r.H_tensor.to_numpy(float);Ar=r.a1.to_numpy(float)
 Gr=r.G_tensor.to_numpy(float);C4r=r.c4.to_numpy(float)

 cases={}
 for name,F,H,A in [
  ("action_all",Fa,Ha,Aa),
  ("archive_FHa1",Fr,Hr,Ar),
  ("archive_F_only",Fr,Ha,Aa),
  ("archive_H_only",Fa,Hr,Aa),
  ("archive_a1_only",Fa,Ha,Ar),
  ("archive_FH",Fr,Hr,Aa),
  ("archive_Ha1",Fa,Hr,Ar),
  ("archive_Fa1",Fr,Ha,Ar),
 ]:
  z=calc(F,H,A,Gr if name=="archive_FHa1" else None,C4r if name=="archive_FHa1" else None)
  arr=z.pop("array");cases[name]=z

 archive_col={}
 if "K_scalar" in r:
  ka=r.K_scalar.to_numpy(float);j=int(np.argmin(ka))
  archive_col={"present":True,"min":float(ka[j]),"x":float(r.x.iloc[j]),"u":float(r.u.iloc[j]),"negative_rows":int(np.sum(ka<=0))}
 else: archive_col={"present":False}

 mismatch={}
 for name,aa,bb in [("F",Fa,Fr),("H",Ha,Hr),("a1",Aa,Ar),("G",q.G_tensor_action.to_numpy(float),Gr),("c4",q.c4_action.to_numpy(float),C4r)]:
  e=aa-bb;j=int(np.argmax(np.abs(e)))
  mismatch[name]={"max_abs":float(np.max(np.abs(e))),"max_scaled":float(np.max(np.abs(e)/np.maximum(1,np.abs(bb)))),"x_at_max":float(q.x.iloc[j])}

 # Radial bands show whether sign loss is purely deep-core derivative conditioning.
 bands={}
 for lo in (0.02,0.03,0.05,0.1,0.2,0.4,0.8,1.0):
  m=q.x.to_numpy(float)>=lo
  z=emit_from_primitives(inp(Fa,Ha,Aa),regularize_photon_root=True).K_scalar.to_numpy(float)[m]
  bands[str(lo)]={"min_action_K_scalar":float(np.min(z)),"negative_rows":int(np.sum(z<=0))}

 rep={"status":"CORE_QUARTIC_KSCALAR_LOCALIZATION_COMPLETE","primitive_mismatch":mismatch,"archive_K_scalar_column":archive_col,"cases":cases,"action_Kscalar_by_xmin":bands}
 # Diagnose.
 arch=cases["archive_FHa1"]["min"];act=cases["action_all"]["min"]
 if arch>0 and act<=0:
  diag="ACTION_TO_PRIMITIVE_DERIVATIVE_SENSITIVITY_CAUSES_SIGN_LOSS"
 elif arch<=0:
  diag="KSCALAR_RECOMPUTATION_FAILS_EVEN_ON_ARCHIVED_PRIMITIVES"
 else:
  diag="KSCALAR_DISCREPANCY_OTHER"
 rep["diagnosis"]=diag
 rep["next_target"]="AUDIT_KSCALAR_DERIVATIVE_CONVERGENCE_AND_ANALYTIC_ACTION_DERIVATIVE" if diag!="KSCALAR_DISCREPANCY_OTHER" else "INSPECT_MEMBER_PROVENANCE"
 OUT.parent.mkdir(parents=True,exist_ok=True);OUT.write_text(json.dumps(rep,indent=2,allow_nan=False)+"\n")
 print(json.dumps(rep,indent=2,allow_nan=False));return 0
if __name__=="__main__":raise SystemExit(main())
