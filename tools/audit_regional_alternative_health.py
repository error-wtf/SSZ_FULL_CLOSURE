#!/usr/bin/env python3
"""Audit newer regional alternatives before any global promotion.

Compares the currently stitched/prestaged regional streams against newer
action/re-emission products already present in the repo:

* central: locked Electric same-action current member from build_onshell_central
  versus the old selected exact-SVT/source representative;
* inner: generated inner_same_action_SVT_H_41of41 versus prestaged candidate;
* core: repaired holonomic EPSY 41of41 versus prestaged core_selected.

This is a diagnostic selection audit only.  A candidate is not promoted merely
because its K/G principal health is better; provenance/common-action/interface
certificates are reported separately.
"""
from __future__ import annotations
import json,sys
from pathlib import Path
import numpy as np,pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from ssz_p5.numerics import module
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central
from ssz_p5.production.regional_coefficients import select_lower,central_selected

OUT=ROOT/"data/generated/qnm_global_diagnostic/REGIONAL_ALTERNATIVE_HEALTH_AUDIT.json"
LS=(6,12,20,42,110,420,1000)

PATHS={
 "inner_prestaged":ROOT/"data/prestaged/direct41/inner_selected_candidate_41of41.csv",
 "inner_generated":ROOT/"data/generated/inner/inner_same_action_SVT_H_41of41.csv",
 "core_prestaged":ROOT/"data/prestaged/direct41/core_selected_41of41.csv",
 "core_repaired_holonomic_epsy":ROOT/"archive/full_working_snapshot/ssz_p5_CORE_REPAIRED_HOLONOMIC_EPSY_41of41_2026-09-16.csv",
 "core_mh_selected":ROOT/"data/regression/ssz_p5_F1b_core_MH_SELECTED_41of41_2026-09-15.csv",
 "core_f3_selected":ROOT/"data/regression/ssz_p5_F3_core_SELECTED_41of41_2026-09-15.csv",
}


def norm(d):
    d=d.copy()
    if "u" not in d and "x" in d:d["u"]=1/d.x.to_numpy(float)
    if "x" not in d and "u" in d:d["x"]=1/d.u.to_numpy(float)
    if "phiprime" not in d and "phi_r" in d:d["phiprime"]=d.phi_r
    if "A0prime" not in d:d["A0prime"]=0.
    for c in ("v5","c3","e3"):
        if c not in d or not np.isfinite(pd.to_numeric(d[c],errors="coerce")).all(): d[c]=0.
    return select_lower(d).sort_values("x").reset_index(drop=True)


def health(d,red,L):
    a=red.canonical_audit(d,float(L))
    K=np.asarray(a["K"],float);G=np.asarray(a["G"],float)
    ke=np.linalg.eigvalsh((K+K.transpose(0,2,1))/2)[:,0]
    cr=np.full(len(d),np.nan)
    for i in np.flatnonzero(ke>0):
        w,U=np.linalg.eigh((K[i]+K[i].T)/2)
        inv=U@np.diag(1/np.sqrt(w))@U.T
        C=inv@((G[i]+G[i].T)/2)@inv
        cr[i]=np.linalg.eigvalsh((C+C.T)/2)[0]
    ik=int(np.argmin(ke))
    return {
      "min_K":float(ke[ik]),"u_at_min_K":float(d.u.iloc[ik]),
      "negative_K_rows":int(np.sum(ke<=0)),
      "min_cr2":float(np.nanmin(cr)) if np.any(np.isfinite(cr)) else None,
      "negative_cr2_rows":int(np.sum(cr<0)),
    }


def candidate(name,d,red,lo=None,hi=None):
    d=norm(d)
    u=d.u.to_numpy(float);m=np.ones(len(d),bool)
    if lo is not None:m&=u>=lo
    if hi is not None:m&=u<hi
    d=d.loc[m].sort_values("x").reset_index(drop=True)
    return {
      "name":name,"rows":len(d),
      "u_range":[float(d.u.min()),float(d.u.max())],
      "per_L":{str(L):health(d,red,L) for L in LS}
    }


def main():
    red=module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    cur=build_onshell_central(ROOT).direct41
    old=central_selected(ROOT)
    cands=[]
    cands.append(candidate("central_old_selected",old,red,.61,.71))
    cands.append(candidate("central_current_electric",cur,red,.61,.71))
    for n in ("inner_prestaged","inner_generated"):
        cands.append(candidate(n,pd.read_csv(PATHS[n]),red,.71,.715))
    for n in ("core_prestaged","core_repaired_holonomic_epsy","core_mh_selected","core_f3_selected"):
        cands.append(candidate(n,pd.read_csv(PATHS[n]),red,.715,None))

    def score(c):
        # lexicographic: number failed L, worst negative row count, worst minK,
        # radial failures.  Diagnostic ranking only.
        pf=list(c["per_L"].values())
        fail=sum(x["min_K"]<=0 or (x["min_cr2"] is not None and x["min_cr2"]<=0) for x in pf)
        neg=sum(x["negative_K_rows"]+x["negative_cr2_rows"] for x in pf)
        mk=min(x["min_K"] for x in pf)
        mc=min(x["min_cr2"] for x in pf if x["min_cr2"] is not None)
        return (fail,neg,-mk,-mc)
    groups={
      "central":["central_old_selected","central_current_electric"],
      "inner":["inner_prestaged","inner_generated"],
      "core":["core_prestaged","core_repaired_holonomic_epsy","core_mh_selected","core_f3_selected"],
    }
    best={}
    by={c["name"]:c for c in cands}
    for g,names in groups.items():
        best[g]=min(names,key=lambda n:score(by[n]))

    report={
      "status":"REGIONAL_ALTERNATIVE_HEALTH_AUDIT_COMPLETE",
      "candidates":cands,
      "best_health_only":best,
      "guard":"best_health_only is NOT a promotion; common-action provenance, interface values/jets, background EOM and full-domain gates remain mandatory",
      "known_provenance_notes":{
        "central_current_electric":"current same-action locked member; healthy only registered 0.62<u<0.70; native-tail defects already established",
        "inner_generated":"INNER_DIRECT_41_CERTIFICATE.json currently FAIL: common-action full F3/F4 and interface values unresolved",
        "core_repaired_holonomic_epsy":"historical/repaired 41 stream; must not be called current action-authoritative without replay certificate"
      }
    }
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps({
      "status":report["status"],"best_health_only":best,
      "summary":{c["name"]:{
        L:{"minK":v["min_K"],"negK":v["negative_K_rows"],"mincr2":v["min_cr2"],"negcr":v["negative_cr2_rows"]}
        for L,v in c["per_L"].items()
      } for c in cands}
    },indent=2))
    return 0
if __name__=="__main__":raise SystemExit(main())
