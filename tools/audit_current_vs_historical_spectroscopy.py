#!/usr/bin/env python3
"""Compare current production spectroscopy to the archived exact central 41-stream.

The archived central Zhang-Kase stream is a *different action representative*
on the same P5 background and must never be silently substituted for the
current locked member.  It is nevertheless a valuable control because the
historical closure records certify a healthy principal structure there.

This audit applies the same full-resolution reducer and the same diagnostic
local/boxed spectroscopy pipeline to both streams on the common central domain.
It asks whether lower-order negative omega^2 / spectral-selectivity patterns are
specific to the current Electric member or already present in the historical
central representative.

Diagnostic control only; no promotion and no same-action claim.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))

from ssz_p5.config import DEFAULT_L  # noqa:E402
from ssz_p5.numerics import module  # noqa:E402
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central  # noqa:E402

HIST=ROOT/"archive/full_working_snapshot/ssz_p5_CENTRAL_SELECTED_41of41_CORRECTED_A5_V12_2026-09-16.csv"
OUT=ROOT/"data/generated/spectral/CURRENT_VS_HISTORICAL_CENTRAL_SPECTROSCOPY.json"


def load_spec():
    p=ROOT/"tools/run_production_spectroscopy_suite.py"
    s=importlib.util.spec_from_file_location("specsuite",p)
    m=importlib.util.module_from_spec(s)
    assert s.loader is not None
    s.loader.exec_module(m)
    return m


def prepare_current():
    d=build_onshell_central(ROOT).direct41.sort_values("x").reset_index(drop=True)
    return d[(d.u>0.62)&(d.u<0.70)].reset_index(drop=True)


def prepare_historical():
    d=pd.read_csv(HIST).sort_values("x").reset_index(drop=True)
    return d[(d.u>0.62)&(d.u<0.70)].reset_index(drop=True)


def summarize_member(d,red,spec,L):
    A=red.canonical_audit(d,int(L))
    K=np.asarray(A["K"],float);G=np.asarray(A["G"],float)
    ke=np.linalg.eigvalsh((K+K.swapaxes(1,2))/2)[:,0]
    cr=[]
    for i in range(len(d)):
        w,U=np.linalg.eigh((K[i]+K[i].T)/2)
        if np.min(w)<=0:
            continue
        inv=U@np.diag(1/np.sqrt(w))@U.T
        C=inv@((G[i]+G[i].T)/2)@inv
        cr.append(float(np.linalg.eigvalsh((C+C.T)/2)[0]))
    box=spec.boxed_suite(d,red,int(L))
    local=spec.local_suite(d,red,int(L)) if L in (6,42,110,420,1000) else None
    return {
        "min_K":float(np.min(ke)),
        "min_cr2_where_K_positive":float(np.min(cr)) if cr else None,
        "boxed_negative_omega2":int(box["physical"]["negative_omega2_count"]),
        "boxed_omega2_min":float(box["physical"]["zero_or_negative_omega2_min"]),
        "boxed_mean_probe_Neff":float(box["comparison"]["mean_probe_effective_mode_count_physical"]),
        "boxed_mean_IPR":float(box["comparison"]["mean_mode_ipr_physical"]),
        "local":local,
    }


def main():
    spec=load_spec()
    red=module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    current=prepare_current()
    hist=prepare_historical()

    result={"current":{},"historical_exact_central":{}}
    for L in DEFAULT_L:
        result["current"][str(L)]=summarize_member(current,red,spec,int(L))
        result["historical_exact_central"][str(L)]=summarize_member(hist,red,spec,int(L))

    comparison={}
    for L in map(str,DEFAULT_L):
        c=result["current"][L];h=result["historical_exact_central"][L]
        comparison[L]={
            "delta_min_K_current_minus_hist":float(c["min_K"]-h["min_K"]),
            "delta_min_cr2_current_minus_hist":float(c["min_cr2_where_K_positive"]-h["min_cr2_where_K_positive"]),
            "boxed_negative_omega2_current":c["boxed_negative_omega2"],
            "boxed_negative_omega2_hist":h["boxed_negative_omega2"],
            "boxed_omega2_min_current":c["boxed_omega2_min"],
            "boxed_omega2_min_hist":h["boxed_omega2_min"],
            "Neff_current_over_hist":float(c["boxed_mean_probe_Neff"]/h["boxed_mean_probe_Neff"]) if h["boxed_mean_probe_Neff"] else None,
        }

    report={
        "scope":"different-action historical control on same P5 central background; no same-action inference",
        "historical_path":str(HIST.relative_to(ROOT)),
        "members":result,
        "comparison":comparison,
        "guard":"Historical exact central stream is a control, not an admissible replacement for the current locked member.",
    }
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"comparison":comparison},indent=2))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
