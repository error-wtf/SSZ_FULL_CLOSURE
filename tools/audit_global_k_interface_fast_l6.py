#!/usr/bin/env python3
"""Fast L=6 global/interface localization for the rebuilt regional stream."""
from __future__ import annotations
import importlib.util,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/"src"))
from ssz_p5.numerics import module
sp=importlib.util.spec_from_file_location("gq",ROOT/"tools/run_global_coupled_qnm_diagnostic.py")
m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)
OUT=ROOT/"data/generated/qnm_global_diagnostic/GLOBAL_K_INTERFACE_FAST_L6.json"

def metrics(a):
 K=np.asarray(a["K"],float);G=np.asarray(a["G"],float);Ks=(K+K.transpose(0,2,1))/2
 ke=np.linalg.eigvalsh(Ks)[:,0];cr=np.full(len(ke),np.nan)
 for i in np.flatnonzero(ke>0):
  w,U=np.linalg.eigh(Ks[i]);inv=U@np.diag(1/np.sqrt(w))@U.T
  C=inv@((G[i]+G[i].T)/2)@inv;cr[i]=np.linalg.eigvalsh((C+C.T)/2)[0]
 return ke,cr
def sm(k,c):
 return {"min_K":float(k.min()),"negK":int(np.sum(k<=0)),"min_cr2":float(np.nanmin(c)) if np.any(np.isfinite(c)) else None,"negcr":int(np.sum(c<0))}
def main():
 g,_=m.assemble();red=module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
 reg=g.production_region.astype(str).to_numpy();ifs=np.flatnonzero(reg[1:]!=reg[:-1])+1
 a=red.canonical_audit(g,6);k,c=metrics(a);wk=int(np.argmin(k))
 stitched={r:sm(k[reg==r],c[reg==r]) for r in dict.fromkeys(reg)}
 iso={}
 for r in dict.fromkeys(reg):
  q=g[reg==r].reset_index(drop=True)
  try:
   kk,cc=metrics(red.canonical_audit(q,6));iso[r]=sm(kk,cc)
  except Exception as e: iso[r]={"error":repr(e)}
 exclusions={}
 for n in (4,8,12,20,40):
  keep=np.ones(len(g),bool)
  for j in ifs:keep[max(0,j-n):min(len(g),j+n)]=False
  exclusions[str(n)]=sm(k[keep],c[keep])
 rep={"status":"FAST_L6_DONE","worst":{"row":wk,"u":float(g.u.iloc[wk]),"x":float(g.x.iloc[wk]),"region":reg[wk],"K":float(k[wk]),"distance_to_interface":int(np.min(np.abs(ifs-wk)))}, "stitched_by_region":stitched,"isolated_regions":iso,"exclusions":exclusions}
 OUT.parent.mkdir(parents=True,exist_ok=True);OUT.write_text(json.dumps(rep,indent=2,allow_nan=False)+"\n");print(json.dumps(rep,indent=2,allow_nan=False));return 0
if __name__=="__main__":raise SystemExit(main())
