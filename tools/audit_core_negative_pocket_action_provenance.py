#!/usr/bin/env python3
"""Trace the quartic-core negative K_scalar pocket back to action-source quality.

The fresh/recomputed K_scalar pocket is centered near x~0.3015.  This audit
inspects the authoritative G4XX action reconstruction exactly in that band:
background solve residuals/condition number, XiH reconstruction, action jets,
and mismatch to the archived F2 primitive table.  It also records whether the
old stored positive K_scalar is an internally independent imported oracle.
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
OUT=ROOT/"data/generated/qnm_global_diagnostic/CORE_NEGATIVE_POCKET_ACTION_PROVENANCE.json"

def main():
 a=pd.read_csv(ACT).sort_values("r_over_rs").reset_index(drop=True)
 r=pd.read_csv(REF).sort_values("x").reset_index(drop=True)
 q=quartic_g5zero_primitives(a)
 # emit on full quartic/native grid using archived c2 to find exact negative rows
 inp=pd.DataFrame({"x":q.x,"u":q.u,"f":a.A_f,"h":a.B_h,"X":a.X,"phi":a.phi,"A0prime":0.,
  "a1":q.a1_action,"c2":r.c2,"c4":q.c4_action,
  "F_tensor":q.F_tensor_action,"G_tensor":q.G_tensor_action,"H_tensor":q.H_tensor_action})
 z=emit_from_primitives(inp,regularize_photon_root=True)
 k=z.K_scalar.to_numpy(float);neg=k<0;ids=np.flatnonzero(neg)
 j=int(np.argmin(k));x0=float(z.x.iloc[j])
 # Include +/- 8 neighboring rows and summary of entire negative pocket.
 lo=max(0,j-8);hi=min(len(a),j+9)
 cols=[c for c in ("background_rank","background_scaled_condition","background_residual_abs","background_residual_rel","XiH_residual","XiH_target","XiH_reconstructed","G4_background","G4X","G4XX","G4phi","G4phiX","G3","G3X","G2","G2X") if c in a]
 rows=[]
 for i in range(lo,hi):
  rec={"row":i,"x":float(a.r_over_rs.iloc[i]),"u":float(a.u.iloc[i]),"K_recomputed":float(k[i]),"K_archived":float(r.K_scalar.iloc[i])}
  rec.update({c:float(a[c].iloc[i]) for c in cols})
  rows.append(rec)
 def stat(c,mask):
  arr=a[c].to_numpy(float)[mask]
  return {"min":float(np.min(arr)),"max":float(np.max(arr)),"median":float(np.median(arr))}
 pocket={}
 if len(ids):
  m=neg
  for c in ("background_scaled_condition","background_residual_abs","background_residual_rel","XiH_residual","G4XX","G4X"):
   if c in a:pocket[c]=stat(c,m)
  pocket["x_range"]=[float(a.r_over_rs.iloc[ids[0]]),float(a.r_over_rs.iloc[ids[-1]])]
  pocket["rows"]=len(ids)
 # Compare source-quality metrics inside vs outside pocket.
 outside=~neg
 contrast={}
 for c in ("background_scaled_condition","background_residual_abs","background_residual_rel","XiH_residual"):
  if c in a:
   contrast[c]={"negative_pocket":stat(c,neg),"outside":stat(c,outside)}
 rep={"status":"CORE_NEGATIVE_POCKET_ACTION_PROVENANCE_COMPLETE","worst":{"row":j,"x":x0,"u":float(a.u.iloc[j]),"K_recomputed":float(k[j]),"K_archived":float(r.K_scalar.iloc[j])},"pocket":pocket,"source_quality_contrast":contrast,"local_rows":rows}
 # rough diagnosis
 badres=False
 if "background_residual_rel" in a and len(ids):
  p=np.nanmedian(a.background_residual_rel.to_numpy(float)[neg]);o=np.nanmedian(a.background_residual_rel.to_numpy(float)[outside])
  badres=bool(p>max(1e-10,100*o))
 rep["diagnosis"]="ACTION_BACKGROUND_SOLVE_DEGRADES_IN_NEGATIVE_POCKET" if badres else "NEGATIVE_POCKET_NOT_EXPLAINED_BY_BACKGROUND_RESIDUAL_QUALITY"
 rep["next_target"]="TRACE_STORED_KSCALAR_ORACLE_GENERATION_AND_FORMULA_PROVENANCE"
 OUT.parent.mkdir(parents=True,exist_ok=True);OUT.write_text(json.dumps(rep,indent=2,allow_nan=False)+"\n")
 print(json.dumps({"status":rep["status"],"worst":rep["worst"],"pocket":rep["pocket"],"contrast":rep["source_quality_contrast"],"diagnosis":rep["diagnosis"]},indent=2,allow_nan=False));return 0
if __name__=="__main__":raise SystemExit(main())
