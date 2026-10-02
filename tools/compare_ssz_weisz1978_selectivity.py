#!/usr/bin/env python3
"""Quantitative comparison of reproduced Weisz-1978 spectra with SSZ local spectra.

This does NOT claim that the mechanisms are identical.  It tests a narrower,
observable statement shared by both calculations: a large mode space can have a
small effective number of spectrally visible modes, with strongly hierarchical
residues.  In SSZ we additionally test radius-dependent redistribution.
"""
from __future__ import annotations
import json
from pathlib import Path
import pandas as pd
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
WREP=ROOT/"data/generated/spectral/weisz1978/WEISZ1978_FK_REPRODUCTION.json"
WMOD=ROOT/"data/generated/spectral/weisz1978/WEISZ1978_FK_MODES.csv"
SREP=ROOT/"data/generated/spectral/SSZ_SPECTROSCOPY_FIXED_OBSERVABLE_V3.json"
SMOD=ROOT/"data/generated/spectral/ssz_local_resolvent/SSZ_LOCAL_RESOLVENT_MODES.csv"
OUT=ROOT/"data/generated/spectral/SSZ_WEISZ1978_QUANTITATIVE_COMPARISON.json"
CSV=ROOT/"data/generated/spectral/SSZ_WEISZ1978_SELECTIVITY_TABLE.csv"


def main():
    wr=json.loads(WREP.read_text())
    sr=json.loads(SREP.read_text())
    wm=pd.read_csv(WMOD)
    sm=pd.read_csv(SMOD)

    wrows=[]
    for (N,M,beta),g in wm.groupby(["N","M","beta"]):
        p=np.sort(g.q0_residue.to_numpy(float))[::-1]
        p=p/max(p.sum(),1e-300)
        wrows.append({
          "system":"Weisz1978_FK","case":f"{int(N)}lambda={int(M)}a beta={beta:g}",
          "top1":float(p[0]),"top2":float(p[:2].sum()),
          "top3":float(p[:3].sum()),"n_eff":float(1/np.sum(p*p)),
          "entropy":float(-np.sum(np.where(p>0,p*np.log(p),0))),
          "mode_count":len(p),
        })

    srows=[]
    for (L,u,obs),g in sm.groupby(["L","u","observable"]):
        p=g.sort_values("mode").normalized_residue.to_numpy(float)
        ps=np.sort(p)[::-1]
        srows.append({
          "system":"SSZ_local","case":f"L={int(L)} u={u:.6f} {obs}",
          "L":int(L),"u":float(u),"observable":obs,
          "top1":float(ps[0]),"top2":float(ps[:2].sum()),
          "top3":float(ps[:3].sum()),"n_eff":float(1/np.sum(p*p)),
          "entropy":float(-np.sum(np.where(p>0,p*np.log(p),0))),
          "mode_count":len(p),
        })

    table=pd.DataFrame(wrows+srows)
    CSV.parent.mkdir(parents=True,exist_ok=True)
    table.to_csv(CSV,index=False)

    w=pd.DataFrame(wrows);s=pd.DataFrame(srows)
    # Paper-inspired sparse-observability region: two modes carry >=95%.
    w_sparse=float(np.mean(w.top2>=0.95))
    s_sparse=float(np.mean(s.top2>=0.95))
    # Very sparse: effective mode count <=2.1.
    w_neff=float(np.mean(w.n_eff<=2.1))
    s_neff=float(np.mean(s.n_eff<=2.1))

    # Radius sensitivity on fixed (L,observable): spread in dominant weight.
    spans=[]
    for (L,obs),g in s.groupby(["L","observable"]):
        spans.append({
          "L":int(L),"observable":obs,
          "top1_span":float(g.top1.max()-g.top1.min()),
          "neff_span":float(g.n_eff.max()-g.n_eff.min()),
          "top2_min":float(g.top2.min()),"top2_max":float(g.top2.max()),
        })

    v3=sr["local_principal_spectroscopy_v3"]
    report={
      "status":"SSZ_WEISZ1978_QUANTITATIVE_SELECTIVITY_COMPARISON_COMPLETE",
      "weisz_reproduction_pass":bool(wr["overall_reproduction_pass"]),
      "ssz_v3_controls_pass":bool(v3["all_controls_pass"]),
      "ssz_v3_robust_residue_reordering":bool(v3["any_robust_residue_reordering"]),
      "selectivity_metrics":{
        "weisz_top1_range":[float(w.top1.min()),float(w.top1.max())],
        "weisz_top2_range":[float(w.top2.min()),float(w.top2.max())],
        "weisz_neff_range":[float(w.n_eff.min()),float(w.n_eff.max())],
        "ssz_top1_range":[float(s.top1.min()),float(s.top1.max())],
        "ssz_top2_range":[float(s.top2.min()),float(s.top2.max())],
        "ssz_neff_range":[float(s.n_eff.min()),float(s.n_eff.max())],
        "fraction_two_modes_ge_95pct":{"weisz":w_sparse,"ssz":s_sparse},
        "fraction_neff_le_2p1":{"weisz":w_neff,"ssz":s_neff},
      },
      "ssz_radius_sensitivity":spans,
      "tests":{
        "weisz_few_visible_modes_reproduced":bool(w_sparse>=0.8 and w_neff>=0.8),
        "ssz_has_weisz_like_sparse_observable_cases":bool(np.any((s.top2>=0.95)&(s.n_eff<=2.1))),
        "ssz_selectivity_changes_with_radius":bool(any(x["top1_span"]>0.05 or x["neff_span"]>0.1 for x in spans)),
        "ssz_residue_order_reversal_survives_v3_controls":bool(v3["all_controls_pass"] and v3["any_robust_residue_reordering"]),
      },
      "interpretation":{
        "supported":"Both numerical systems exhibit sparse observable spectral weight; SSZ additionally shows robust radius-dependent redistribution in the tested healthy window.",
        "not_supported":"Identity of microscopic mechanism, Anderson localization in SSZ, or global SSZ QNM/retarded-Green spectroscopy."
      },
      "table_csv":str(CSV.relative_to(ROOT)),
    }
    OUT.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps(report,indent=2,allow_nan=False))
    return 0 if all(report["tests"].values()) else 2

if __name__=="__main__": raise SystemExit(main())
