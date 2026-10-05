#!/usr/bin/env python3
"""Regenerate the quartic G5=0 punctured core directly from action jets.

The deep-core 3-field stream is a documented implementation failure, while the
exact MH invariants remain positive.  The available action bridge already
reconstructs F,G,H,a1,c4 from the authoritative G3/G4/G4X/G4XX action table.

This audit takes the locked selected c2 transverse control (the allowed
background-null lower representative), emits all 41 slots anew with the current
MH primitive emitter and JET9D8, and tests the quartic core x>=0.02 where G5=0.
The x<0.02 full-G5 subcore is deliberately excluded.

It then compares:
  * exact primitive invariant K_scalar,
  * 3-field canonical K/cr2,
  * unreduced generalized-psi local descriptor growth,
against the old prestaged core stream on the identical quartic domain.
"""
from __future__ import annotations
import json,sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.linalg import eig

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from ssz_p5.coefficients.mh_action_primitives import quartic_g5zero_primitives  # noqa:E402
from ssz_p5.coefficients.mh_general_primitives import emit_from_primitives  # noqa:E402
from ssz_p5.numerics import module  # noqa:E402
from ssz_p5.reducer.unreduced_descriptor import generalized_psi_descriptor  # noqa:E402

ACT=ROOT/"data/production/ssz_p5_F1b_G4XX_transverse_core_candidate_2026-09-14.csv"
REF=ROOT/"data/prestaged/direct41/core_selected_41of41.csv"
OUTDIR=ROOT/"data/generated/qnm_global_diagnostic"
CSV=OUTDIR/"CORE_QUARTIC_ACTION_REGENERATED_41.csv"
REPORT=OUTDIR/"CORE_QUARTIC_ACTION_REGEN_AUDIT.json"
LS=(6,12,20,42,110,420,1000)

def kg(stream,L):
 red=module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
 a=red.canonical_audit(stream,float(L));K=np.asarray(a["K"],float);G=np.asarray(a["G"],float)
 Ks=(K+K.transpose(0,2,1))/2;ke=np.linalg.eigvalsh(Ks)[:,0];cr=np.full(len(ke),np.nan)
 for i in np.flatnonzero(ke>0):
  w,U=np.linalg.eigh(Ks[i]);inv=U@np.diag(1/np.sqrt(w))@U.T
  C=inv@((G[i]+G[i].T)/2)@inv;cr[i]=np.linalg.eigvalsh((C+C.T)/2)[0]
 return {"min_K":float(ke.min()),"negK":int(np.sum(ke<=0)),"min_cr2":float(np.nanmin(cr)) if np.any(np.isfinite(cr)) else None,"negcr":int(np.sum(cr<0))}

def descriptor_kr0(stream,L,ns=32):
 P=generalized_psi_descriptor(stream,float(L))
 ids=np.unique(np.linspace(4,len(stream)-5,min(ns,len(stream)-8)).round().astype(int))
 mx=-np.inf;worst=None
 for i in ids:
  sample=next(iter(P.values()));nf=sample.shape[1]
  A=[np.zeros((nf,nf),complex) for _ in range(3)]
  for (pt,pr),B in P.items():
   if pr==0:A[pt]+=B[i]
  n=nf;identity_mat=np.eye(n,dtype=complex);Z=np.zeros_like(identity_mat)
  z,V=eig(np.block([[-A[1],-A[0]],[identity_mat,Z]]),np.block([[A[2],Z],[Z,identity_mat]]),right=True,check_finite=False)
  good=np.isfinite(z.real)&np.isfinite(z.imag)&(np.abs(z)<1e8)
  z=z[good];V=V[:,good]
  for j,lam in enumerate(z):
   y=V[n:,j];n2=np.linalg.norm(A[2]@y);n1=abs(lam)*np.linalg.norm(A[1]@y);n0=np.linalg.norm(A[0]@y)
   part=(abs(lam)**2*n2)/max(n0+n1+(abs(lam)**2)*n2,1e-300)
   if part>1e-6 and lam.real>mx:
    mx=float(lam.real);worst={"row":int(i),"x":float(stream.x.iloc[i]),"u":float(stream.u.iloc[i]),"lambda_re":float(lam.real),"lambda_im":float(lam.imag),"kinetic_participation":float(part)}
 return {"max_Re_lambda_kinetic_active":float(mx),"worst":worst}

def main():
 act=pd.read_csv(ACT).sort_values("r_over_rs").reset_index(drop=True)
 ref=pd.read_csv(REF).sort_values("x").reset_index(drop=True)
 prim=quartic_g5zero_primitives(act)
 # Restrict to core production and G5=0 quartic domain x>=0.02.
 mask=(prim.x.to_numpy(float)>=0.02)&(prim.u.to_numpy(float)>=0.715)
 prim=prim.loc[mask].reset_index(drop=True)
 aa=act.sort_values("r_over_rs").reset_index(drop=True).loc[mask].reset_index(drop=True)
 # Same action/reference native grid; selected c2 is the locked transverse control.
 rx=ref.x.to_numpy(float)
 c2=np.interp(prim.x.to_numpy(float),rx,ref.c2.to_numpy(float))
 inp=pd.DataFrame({
  "x":prim.x,"u":prim.u,
  "f":aa.A_f,"h":aa.B_h,"X":aa.X,
  "phi":aa.phi,"A0prime":0.0,
  "a1":prim.a1_action,"c2":c2,"c4":prim.c4_action,
  "F_tensor":prim.F_tensor_action,"G_tensor":prim.G_tensor_action,"H_tensor":prim.H_tensor_action,
 })
 out=emit_from_primitives(inp,regularize_photon_root=True)
 out["production_region"]="punctured_H_core_quartic_action_regen"
 OUTDIR.mkdir(parents=True,exist_ok=True);out.to_csv(CSV,index=False)

 old=ref[(ref.x>=0.02)&(ref.u>=0.715)].sort_values("x").reset_index(drop=True)
 # interpolate old only for comparison if edge serialization differs
 results={}
 for L in LS:
  results[str(L)]={"regenerated":kg(out,L),"old_reference":kg(old,L)}
 # Exact action-derived invariant emitted with same primitives.
 report={
  "status":"CORE_QUARTIC_ACTION_REGEN_COMPLETE",
  "scope":"G5=0 quartic core x>=0.02 only; full-G5 x<0.02 excluded",
  "rows":len(out),"x_range":[float(out.x.min()),float(out.x.max())],"u_range":[float(out.u.min()),float(out.u.max())],
  "primitive_invariants":{
    "F_min":float(out.F_tensor.min()),"G_min":float(out.G_tensor.min()),"H_min":float(out.H_tensor.min()),"K_scalar_min":float(out.K_scalar.min())
  },
  "per_L":results,
  "descriptor_kr0_L6":descriptor_kr0(out,6),
  "descriptor_kr0_L1000":descriptor_kr0(out,1000),
 }
 allK=all(v["regenerated"]["min_K"]>0 for v in results.values())
 allR=all(v["regenerated"]["min_cr2"] is not None and v["regenerated"]["min_cr2"]>0 for v in results.values())
 report["all_regenerated_K_pass"]=allK;report["all_regenerated_radial_pass"]=allR
 if allK and allR:
  report["verdict"]="QUARTIC_CORE_DIRECT_REGEN_HEALTHY"
  report["next_target"]="IMPLEMENT_FULL_G5_SUBCORE_PRIMITIVE_BRIDGE_AND_REGENERATE_X_LT_0P02"
 else:
  report["verdict"]="QUARTIC_CORE_DIRECT_REGEN_STILL_FAILS"
  report["next_target"]="LOCALIZE_PRIMITIVE_OR_LOWER_CONTROL_CAUSE_BEFORE_FULL_G5_SUBCORE"
 REPORT.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
 print(json.dumps(report,indent=2,allow_nan=False))
 return 0
if __name__=="__main__":raise SystemExit(main())
