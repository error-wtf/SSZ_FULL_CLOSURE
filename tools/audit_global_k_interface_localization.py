#!/usr/bin/env python3
"""Localize the huge global finite-l K/cr2 failure by region and interface.

The global coupled QNM run after rebuilding the complete outer H+SVT stream
produced O(1e6) negative kinetic eigenvalues.  This audit distinguishes:
  1. genuine within-region operator pathology;
  2. derivative contamination caused by differentiating a stitched stream
     across coefficient/value discontinuities at regional interfaces.

For each L it computes:
  * the canonical audit on the stitched global stream;
  * minima/counts by production_region;
  * the same minima after removing progressively wider row bands around every
    interface;
  * independent per-region canonical audits, so JET9D8 never crosses an
    interface;
  * direct slot jumps at each interface;
  * nearest-interface distance of the worst K and c_r^2 rows.

Diagnostic only; no production mutation.
"""
from __future__ import annotations
import json,sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from ssz_p5.numerics import module  # noqa:E402
import importlib.util
_spec=importlib.util.spec_from_file_location("global_qnm_diag", ROOT/"tools/run_global_coupled_qnm_diagnostic.py")
_mod=importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)
assemble=_mod.assemble

OUT=ROOT/"data/generated/qnm_global_diagnostic/GLOBAL_K_INTERFACE_LOCALIZATION.json"
CSV=ROOT/"data/generated/qnm_global_diagnostic/GLOBAL_K_INTERFACE_LOCALIZATION.csv"
LS=(6,12,20,42,110,420,1000)
MARGINS=(0,2,4,8,12,20)
SLOTS=(
 "a1","a2","a3","a4","a5","a6","a7","a8","a9",
 "b1","b2","b3","b4","b5",
 "c1","c2","c3","c4","c5","c6",
 "d1","d2","d3","d4",
 "e1","e2","e3","e4",
 "v1","v2","v3","v4","v5","v6","v7","v8","v9","v10","v11","v12","v13"
)

def metrics(a):
    K=np.asarray(a["K"],float);G=np.asarray(a["G"],float)
    Ks=(K+K.transpose(0,2,1))/2
    ke=np.linalg.eigvalsh(Ks)[:,0]
    cr=np.full(len(K),np.nan)
    for i in np.flatnonzero(ke>0):
        w,U=np.linalg.eigh(Ks[i])
        inv=U@np.diag(1/np.sqrt(w))@U.T
        C=inv@((G[i]+G[i].T)/2)@inv
        cr[i]=np.linalg.eigvalsh((C+C.T)/2)[0]
    return ke,cr

def summarize(ke,cr,idx=None):
    if idx is None: idx=np.arange(len(ke))
    idx=np.asarray(idx,int)
    qk=ke[idx];qc=cr[idx]
    return {
      "rows":int(len(idx)),
      "min_K":float(np.min(qk)) if len(idx) else None,
      "negative_K_rows":int(np.sum(qk<=0)),
      "min_cr2_where_defined":float(np.nanmin(qc)) if np.any(np.isfinite(qc)) else None,
      "negative_cr2_rows_where_defined":int(np.sum(qc<0)),
    }

def interface_positions(g):
    reg=g.production_region.astype(str).to_numpy()
    return np.flatnonzero(reg[1:]!=reg[:-1])+1

def nearest_iface(i,ifs):
    if len(ifs)==0:return None,None
    j=int(ifs[np.argmin(np.abs(ifs-i))])
    return j,int(abs(j-i))

def main():
    g,_=assemble()
    red=module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    x=g.x.to_numpy(float);u=g.u.to_numpy(float)
    ifs=interface_positions(g)
    regions=list(dict.fromkeys(g.production_region.astype(str).tolist()))

    iface=[]
    for j in ifs:
        left=g.iloc[j-1];right=g.iloc[j]
        jumps={}
        for s in SLOTS:
            lv=float(left[s]);rv=float(right[s])
            scale=max(abs(lv),abs(rv),1e-12)
            jumps[s]={"abs":abs(rv-lv),"relative":abs(rv-lv)/scale}
        top=sorted(jumps.items(),key=lambda kv:kv[1]["relative"],reverse=True)[:12]
        iface.append({
          "left_row":int(j-1),"right_row":int(j),
          "left_region":str(left.production_region),"right_region":str(right.production_region),
          "x_left":float(left.x),"x_right":float(right.x),
          "u_left":float(left.u),"u_right":float(right.u),
          "top_relative_slot_jumps":[{"slot":k,**v} for k,v in top],
          "max_relative_slot_jump":float(max(v["relative"] for v in jumps.values())),
          "max_abs_slot_jump":float(max(v["abs"] for v in jumps.values())),
        })

    rows=[]
    perL={}
    for L in LS:
        a=red.canonical_audit(g,float(L));ke,cr=metrics(a)
        wk=int(np.argmin(ke))
        wc=int(np.nanargmin(cr)) if np.any(np.isfinite(cr)) else None
        jk,dk=nearest_iface(wk,ifs)
        jc,dc=nearest_iface(wc,ifs) if wc is not None else (None,None)
        rec={
          "global":summarize(ke,cr),
          "worst_K":{
            "row":wk,"x":float(x[wk]),"u":float(u[wk]),
            "region":str(g.production_region.iloc[wk]),"value":float(ke[wk]),
            "nearest_interface_row":jk,"distance_rows":dk,
          },
          "worst_cr2":None if wc is None else {
            "row":wc,"x":float(x[wc]),"u":float(u[wc]),
            "region":str(g.production_region.iloc[wc]),"value":float(cr[wc]),
            "nearest_interface_row":jc,"distance_rows":dc,
          },
          "stitched_by_region":{},
          "interface_exclusion":{},
          "isolated_region_audits":{},
        }
        for region in regions:
            ids=np.flatnonzero(g.production_region.astype(str).to_numpy()==region)
            rec["stitched_by_region"][region]=summarize(ke,cr,ids)
            rows.append({"L":L,"mode":"stitched_region","region":region,**rec["stitched_by_region"][region]})

        allids=np.arange(len(g))
        for m in MARGINS:
            keep=np.ones(len(g),bool)
            for j in ifs:
                keep[max(0,j-m):min(len(g),j+m)]=False
            ids=allids[keep]
            rec["interface_exclusion"][str(m)]=summarize(ke,cr,ids)
            rows.append({"L":L,"mode":f"exclude_{m}","region":"ALL",**rec["interface_exclusion"][str(m)]})

        for region in regions:
            q=g[g.production_region.astype(str)==region].copy().reset_index(drop=True)
            if len(q)<25:
                rec["isolated_region_audits"][region]={"status":"TOO_SHORT","rows":len(q)}
                continue
            try:
                aa=red.canonical_audit(q,float(L));kk,cc=metrics(aa)
                ss=summarize(kk,cc)
                ss["status"]="OK"
            except Exception as e:
                ss={"status":"ERROR","error":repr(e),"rows":len(q)}
            rec["isolated_region_audits"][region]=ss
            if ss.get("status")=="OK":
                rows.append({"L":L,"mode":"isolated_region","region":region,**{k:v for k,v in ss.items() if k!="status"}})
        perL[str(L)]=rec

    # Decide whether O(1e6) pathology is interface-driven.
    stitched_worst=min(perL[str(L)]["global"]["min_K"] for L in LS)
    iso_mins=[]
    for L in LS:
        for r,v in perL[str(L)]["isolated_region_audits"].items():
            if v.get("status")=="OK": iso_mins.append(v["min_K"])
    isolated_worst=min(iso_mins) if iso_mins else None
    if isolated_worst is not None and abs(stitched_worst)>100*max(abs(isolated_worst),1):
        diagnosis="GLOBAL_DERIVATIVE_INTERFACE_AMPLIFICATION_DOMINATES"
        next_target="BUILD_PIECEWISE_JET_AWARE_GLOBAL_REDUCTION_WITH_EXPLICIT_INTERFACE_MATCHING"
    else:
        diagnosis="LARGE_INSTABILITY_PERSISTS_WITHIN_REGIONS"
        next_target="LOCALIZE_ACTION_OPERATOR_CAUSE_IN_WORST_ISOLATED_REGION"

    report={
      "status":"GLOBAL_K_INTERFACE_LOCALIZATION_COMPLETE",
      "diagnosis":diagnosis,
      "next_target":next_target,
      "rows":len(g),"regions":regions,
      "interfaces":iface,
      "global_worst_min_K":float(stitched_worst),
      "isolated_region_worst_min_K":None if isolated_worst is None else float(isolated_worst),
      "per_L":perL,
    }
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    pd.DataFrame(rows).to_csv(CSV,index=False)
    print(json.dumps({
      "status":report["status"],"diagnosis":diagnosis,
      "global_worst_min_K":stitched_worst,
      "isolated_region_worst_min_K":isolated_worst,
      "interfaces":[{
        "left":q["left_region"],"right":q["right_region"],
        "max_relative_slot_jump":q["max_relative_slot_jump"]
      } for q in iface],
      "per_L":{L:{
        "worst_K":v["worst_K"],
        "exclude20":v["interface_exclusion"]["20"],
        "isolated":v["isolated_region_audits"]
      } for L,v in perL.items()}
    },indent=2,allow_nan=False))
    return 0

if __name__=="__main__": raise SystemExit(main())
