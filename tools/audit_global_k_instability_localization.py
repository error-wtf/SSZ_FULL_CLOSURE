#!/usr/bin/env python3
"""Localize global K/R instability by region vs stitching/interface effects.

Requires the rebuilt complete outer H+SVT stream. It compares:
  1) each selected regional stream audited standalone on its own rows;
  2) the stitched global stream audited on the same rows;
  3) narrow windows around every production interface.

If a region is healthy standalone but unhealthy only after stitching, the
instability is derivative/interface-induced. If it is already unhealthy
standalone, it is intrinsic to that selected regional member.
"""
from __future__ import annotations
import json,sys
from pathlib import Path
import numpy as np,pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))

from ssz_p5.numerics import module
from ssz_p5.production.regions import OUTER_START,CENTRAL_START,INNER_START,CORE_START
from ssz_p5.production.regional_coefficients import select_lower,central_selected
from ssz_p5.production.sources import SOURCE_REGISTRY

OUTDIR=ROOT/"data/generated/qnm_global_diagnostic"
OUT=OUTDIR/"GLOBAL_K_INSTABILITY_LOCALIZATION.json"
GLOBAL=OUTDIR/"GLOBAL_REGIONAL_REFERENCE_41.csv"
OUTER=OUTDIR/"OUTER_COMPLETE_H_SVT_41.csv"
LS=(6,12,20,42)


def norm(path,region):
    d=pd.read_csv(path).copy()
    if "u" not in d and "x" in d:d["u"]=1/d.x.to_numpy(float)
    if "x" not in d and "u" in d:d["x"]=1/d.u.to_numpy(float)
    if "phiprime" not in d and "phi_r" in d:d["phiprime"]=d.phi_r
    if "A0prime" not in d:d["A0prime"]=0.
    for c in ("v5","c3","e3"):
        if c not in d or not np.isfinite(pd.to_numeric(d[c],errors="coerce")).all(): d[c]=0.
    d=select_lower(d)
    d["production_region"]=region
    return d.sort_values("x").reset_index(drop=True)


def health(d,red,L):
    a=red.canonical_audit(d,float(L))
    K=np.asarray(a["K"],float);G=np.asarray(a["G"],float)
    ks=np.linalg.eigvalsh((K+K.transpose(0,2,1))/2)[:,0]
    cr=np.full(len(d),np.nan)
    for i in np.flatnonzero(ks>0):
        w,U=np.linalg.eigh((K[i]+K[i].T)/2)
        inv=U@np.diag(1/np.sqrt(w))@U.T
        C=inv@((G[i]+G[i].T)/2)@inv
        cr[i]=np.linalg.eigvalsh((C+C.T)/2)[0]
    ik=int(np.argmin(ks))
    ic=int(np.nanargmin(cr)) if np.any(np.isfinite(cr)) else None
    return {
      "min_K":float(ks[ik]),"u_at_min_K":float(d.u.iloc[ik]),
      "negative_K_rows":int(np.sum(ks<=0)),
      "min_cr2":float(cr[ic]) if ic is not None else None,
      "u_at_min_cr2":float(d.u.iloc[ic]) if ic is not None else None,
      "negative_cr2_rows":int(np.sum(cr<0)),
    }


def deriv_jump(g,u0,col):
    x=g.x.to_numpy(float); y=g[col].to_numpy(float); u=g.u.to_numpy(float)
    j=int(np.argmin(abs(u-u0)))
    # one-sided linear slopes over nearby valid points
    lo=np.arange(max(0,j-5),j)
    hi=np.arange(j+1,min(len(g),j+6))
    if len(lo)<2 or len(hi)<2:return None
    sl=np.polyfit(x[lo],y[lo],1)[0]; sr=np.polyfit(x[hi],y[hi],1)[0]
    return {"left":float(sl),"right":float(sr),"jump":float(sr-sl)}


def main():
    assert GLOBAL.exists() and OUTER.exists()
    reg=SOURCE_REGISTRY
    regions={
      "weak_exterior_H":norm(ROOT/reg["weak_exterior_H"]["coeff_reference"],"weak_exterior_H"),
      "outer_same_action_H_SVT":norm(OUTER,"outer_same_action_H_SVT"),
      "central_exact_SVT":central_selected(ROOT),
      "inner_same_action_SVT_H":norm(ROOT/reg["inner_same_action_SVT_H"]["coeff_reference"],"inner_same_action_SVT_H"),
      "punctured_H_core":norm(ROOT/reg["punctured_H_core"]["coeff_reference"],"punctured_H_core"),
    }
    bounds={
      "weak_exterior_H":(None,OUTER_START),
      "outer_same_action_H_SVT":(OUTER_START,CENTRAL_START),
      "central_exact_SVT":(CENTRAL_START,INNER_START),
      "inner_same_action_SVT_H":(INNER_START,CORE_START),
      "punctured_H_core":(CORE_START,None),
    }
    for name,d in list(regions.items()):
        lo,hi=bounds[name];m=np.ones(len(d),bool)
        if lo is not None:m&=d.u.to_numpy(float)>=lo
        if hi is not None:m&=d.u.to_numpy(float)<hi
        regions[name]=d.loc[m].sort_values("x").reset_index(drop=True)

    g=pd.read_csv(GLOBAL).sort_values("x").reset_index(drop=True)
    red=module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    result={}
    for name,d in regions.items():
        gm=g[g.production_region.astype(str)==name].copy().sort_values("x").reset_index(drop=True)
        per={}
        for L in LS:
            per[str(L)]={
              "standalone":health(d,red,L),
              "stitched_same_rows":health(gm,red,L) if len(gm)>=12 else None,
            }
        result[name]={"rows_standalone":len(d),"rows_global":len(gm),"per_L":per}

    interfaces={}
    for tag,u0 in [("weak_outer",OUTER_START),("outer_central",CENTRAL_START),
                   ("central_inner",INNER_START),("inner_core",CORE_START)]:
        # +/- ~12 nearest global rows
        ids=np.argsort(abs(g.u.to_numpy(float)-u0))[:28]
        w=g.iloc[np.sort(ids)].sort_values("x").reset_index(drop=True)
        interfaces[tag]={
          "u0":u0,
          "regions":sorted(set(map(str,w.production_region))),
          "per_L":{str(L):health(w,red,L) for L in LS},
          "slot_derivative_jumps":{
            c:deriv_jump(g,u0,c) for c in ("v1","v9","v10","a1","a2","c2","c4","e1","e2")
            if c in g
          }
        }

    # classification
    intrinsic=[];stitch=[]
    for name,r in result.items():
        for L,p in r["per_L"].items():
            s=p["standalone"];q=p["stitched_same_rows"]
            if s["min_K"]<0 or (s["min_cr2"] is not None and s["min_cr2"]<0):
                intrinsic.append({"region":name,"L":L,"standalone":s})
            elif q and (q["min_K"]<0 or (q["min_cr2"] is not None and q["min_cr2"]<0)):
                stitch.append({"region":name,"L":L,"stitched":q})

    report={
      "status":"GLOBAL_K_INSTABILITY_LOCALIZED",
      "regions":result,"interfaces":interfaces,
      "intrinsic_failures":intrinsic,
      "stitch_only_failures":stitch,
      "diagnosis":(
        "INTRINSIC_REGIONAL_FAILURES_PRESENT" if intrinsic else
        "STITCH_INTERFACE_DERIVATIVE_FAILURE_DOMINANT" if stitch else
        "NO_FAILURE_IN_TESTED_DECOMPOSITION"
      )
    }
    OUT.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps({
      "status":report["status"],"diagnosis":report["diagnosis"],
      "intrinsic":[{"region":x["region"],"L":x["L"],"minK":x["standalone"]["min_K"],"mincr2":x["standalone"]["min_cr2"]} for x in intrinsic],
      "stitch_only":[{"region":x["region"],"L":x["L"],"minK":x["stitched"]["min_K"],"mincr2":x["stitched"]["min_cr2"]} for x in stitch],
      "interface_minima":{k:{L:v["per_L"][L] for L in v["per_L"]} for k,v in interfaces.items()}
    },indent=2))
    return 0

if __name__=="__main__":raise SystemExit(main())
