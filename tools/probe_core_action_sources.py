#!/usr/bin/env python3
"""Inspect action-level core source schemas and derive emitter readiness.

This is a provenance/schema probe only. It determines whether the locked
action-level G4XX/G5 core files already contain enough information to regenerate
the primitive inputs required by mh_general_primitives.emit_from_primitives.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/generated/qnm_global_diagnostic/CORE_ACTION_SOURCE_SCHEMA_PROBE.json"

FILES={
 "f2_core":ROOT/"data/authoritative/ssz_p5_F2_core_punctured_horndeski_unreduced_39of41_2026-09-14.csv",
 "g4xx":ROOT/"data/production/ssz_p5_F1b_G4XX_transverse_core_candidate_2026-09-14.csv",
 "g5":ROOT/"data/production/ssz_p5_F1b_FULL_G5_subcore_onshell_candidate_2026-09-14.csv",
 "center_handover":ROOT/"data/production/ssz_p5_F1b_FINAL_Cinf_center_to_punctured_handover_2026-09-14.csv",
}

TOKENS=("G4","G4X","G4XX","G4phi","G4phiphi","G5","G5X","G5phi","a1","c2","c4",
        "F_tensor","G_tensor","H_tensor","K_scalar","f","h","phi","phi_r","phiprime","x","u")


def stats(d,c):
    x=pd.to_numeric(d[c],errors="coerce").to_numpy(float)
    fin=np.isfinite(x)
    return {
      "finite_fraction":float(np.mean(fin)),
      "min":float(np.nanmin(x)) if np.any(fin) else None,
      "max":float(np.nanmax(x)) if np.any(fin) else None,
      "max_abs":float(np.nanmax(np.abs(x))) if np.any(fin) else None,
    }


def main():
    rep={}
    dfs={}
    for name,p in FILES.items():
        d=pd.read_csv(p)
        dfs[name]=d
        wanted=[c for c in d.columns if any(t.lower() in c.lower() for t in TOKENS)]
        rep[name]={
          "path":str(p.relative_to(ROOT)),
          "rows":len(d),
          "column_count":len(d.columns),
          "columns":list(d.columns),
          "relevant_columns":wanted,
          "stats":{c:stats(d,c) for c in wanted if pd.api.types.is_numeric_dtype(d[c])},
        }

    # Alignment and direct common columns.
    pairs={}
    for a in dfs:
        for b in dfs:
            if a>=b: continue
            da,db=dfs[a],dfs[b]
            common=sorted(set(da.columns)&set(db.columns))
            rec={"common_columns":common}
            for coord in ("x","u"):
                if coord in da and coord in db and len(da)==len(db):
                    rec[f"{coord}_max_abs_delta"]=float(np.max(np.abs(
                        pd.to_numeric(da[coord],errors="coerce").to_numpy(float)-
                        pd.to_numeric(db[coord],errors="coerce").to_numpy(float))))
            pairs[f"{a}__{b}"]=rec

    required={"a1","c2","c4","F_tensor","G_tensor","H_tensor"}
    readiness={}
    for name,d in dfs.items():
        readiness[name]={
          "direct_primitive_inputs_present":sorted(required & set(d.columns)),
          "missing_primitive_inputs":sorted(required-set(d.columns)),
          "has_all_primitives":required.issubset(d.columns),
          "action_level_columns":[c for c in d.columns if c.startswith(("G2","G3","G4","G5","K","F"))]
        }

    report={
      "status":"CORE_ACTION_SOURCE_SCHEMA_PROBE_COMPLETE",
      "files":rep,"pair_alignment":pairs,"emitter_readiness":readiness,
      "required_by_mh_general_primitives":sorted(required),
    }
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps(report,indent=2,allow_nan=False))
    return 0

if __name__=="__main__": raise SystemExit(main())
