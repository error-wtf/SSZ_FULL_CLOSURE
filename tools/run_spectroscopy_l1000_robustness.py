#!/usr/bin/env python3
"""Targeted robustness audit for the finite-window spectroscopy signal.

The broad finite-window suite identified a Weisz-like concentration candidate
only at L=1000.  This audit stress-tests that signal against:
  * native, stride-2 and stride-3 grids;
  * 8, 10 and 12 radial bins;
  * registered, inner-trim, outer-trim and both-trim boxes;
  * L=420 as a nearby high-L negative/control comparison.

The statistic is mode-label independent:
  top2 = fraction of local binned K-residue carried by the two strongest modes,
  N_eff = 1/sum_n P_n^2.

Edge bins are excluded.  This remains finite-box spectroscopy, not a global QNM
claim; the direct-global certificate gate is intentionally not bypassed.
"""
from __future__ import annotations
import json,sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))

from ssz_p5.numerics import module
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central
from ssz_p5.qnm.native_window import downsample_native,solve_near_zero_box_spectrum,binned_residue_summary

OUT=ROOT/"data/generated/spectral/SPECTROSCOPY_L1000_ROBUSTNESS.json"
LS=(420,1000)
STRIDES=(1,2,3)
BINS=(8,10,12)
WINDOWS={
 "registered":(0.62,0.70),
 "inner_trim":(0.625,0.70),
 "outer_trim":(0.62,0.695),
 "both_trim":(0.625,0.695),
}
MODES=14
COMPARE=8

def solve(r,u,K,G,S,M,lo,hi,stride):
    m=(u>lo)&(u<hi)
    arr=tuple(np.asarray(A)[m] for A in (r,u,K,G,S,M))
    if stride>1: arr=downsample_native(*arr,stride)
    return solve_near_zero_box_spectrum(*arr,modes=MODES)

def concentration(spec,bins):
    _,_,P,top2,neff=binned_residue_summary(spec.r,spec.weights[:COMPARE],bins=bins)
    if len(top2)>2:
        top2i=top2[1:-1]; neffi=neff[1:-1]; Pi=P[:,1:-1]
    else:
        top2i=top2; neffi=neff; Pi=P
    j=int(np.argmax(top2i)); k=int(np.argmin(neffi))
    return {
      "top2_max":float(np.max(top2i)),
      "neff_min":float(np.min(neffi)),
      "top2_bin":j,
      "neff_bin":k,
      "dominant_mode_at_top2":int(np.argmax(Pi[:,j])),
      "strong":bool(np.max(top2i)>=0.75 and np.min(neffi)<=2.5),
    }

def main():
    b=build_onshell_central(ROOT)
    st=b.direct41.sort_values("x").reset_index(drop=True)
    r=st.x.to_numpy(float);u=st.u.to_numpy(float)
    red=module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    out={}
    for L in LS:
        a=red.canonical_audit(st,L)
        K,G,S,M=(np.asarray(a[x],float) for x in ("K","G","S","M"))
        cases={}
        strong_count=0;total=0
        for w,(lo,hi) in WINDOWS.items():
            for stride in STRIDES:
                sp=solve(r,u,K,G,S,M,lo,hi,stride)
                for bins in BINS:
                    rec=concentration(sp,bins)
                    rec["finite_box_tachyon"]=bool(np.any(sp.omega2<0))
                    key=f"{w}|stride{stride}|bins{bins}"
                    cases[key]=rec
                    strong_count+=int(rec["strong"]);total+=1
        vals_top=[x["top2_max"] for x in cases.values()]
        vals_n=[x["neff_min"] for x in cases.values()]
        out[str(L)]={
          "cases":cases,
          "strong_cases":strong_count,
          "total_cases":total,
          "strong_fraction":float(strong_count/total),
          "top2_max_range":[float(min(vals_top)),float(max(vals_top))],
          "neff_min_range":[float(min(vals_n)),float(max(vals_n))],
          "robust_strong_majority":bool(strong_count>=0.75*total),
          "robust_strong_all":bool(strong_count==total),
        }
    l1000=out["1000"]; l420=out["420"]
    if l1000["robust_strong_majority"] and not l420["robust_strong_majority"]:
        status="L1000_FINITE_WINDOW_CONCENTRATION_ROBUST_WITH_L420_CONTROL"
    elif l1000["robust_strong_majority"]:
        status="HIGH_L_FINITE_WINDOW_CONCENTRATION_ROBUST_NOT_L_SPECIFIC"
    else:
        status="L1000_FINITE_WINDOW_CONCENTRATION_NOT_ROBUST"
    payload={
      "status":status,
      "semantics":"finite-window binned K-residue robustness only; not global QNM residues",
      "thresholds":{"top2_min":0.75,"neff_max":2.5,"majority_fraction":0.75},
      "per_L":out,
      "global_claim":"BLOCKED_UNTIL_DIRECT_GLOBAL_KRGM_CERTIFICATE_AND_GLOBAL_OPERATOR",
    }
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(payload,indent=2,allow_nan=False)+"\n")
    print(json.dumps(payload,indent=2,allow_nan=False))
    return 0

if __name__=="__main__":raise SystemExit(main())
