#!/usr/bin/env python3
"""Healthy local-k spectroscopy on the SSZ production window.

The boxed diagnostic contains robust negative omega^2 at low radial momentum and
its positive eigenfrequencies are not yet numerically converged.  To obtain a
cleaner spectroscopy test without artificial radial boundaries, this audit works
with the full-resolution reduced local 3x3 Hermitian pencil

    H(k,r) = k^2 G(r) + i k S(r) - M(r),
    H v = omega^2 K(r) v.

For each L and radius it finds the smallest positive k at which all three local
omega^2 branches are positive.  It then chooses one common safe k for the whole
sampled production window, tracks mode identity by K-whitened eigenvector
overlap, and measures frequency-ratio and probe-weight changes.

This is a healthy local dispersion / spectral-weight diagnostic.  It is not a
center-to-infinity QNM or Green-function residue calculation.
"""

from __future__ import annotations

import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))

from ssz_p5.numerics import module  # noqa:E402
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central  # noqa:E402

OUT=ROOT/"data/generated/spectral/HEALTHY_LOCAL_K_SPECTROSCOPY.json"
LS=(6,12,20,42,110,420,1000)
TARGET_U=np.linspace(0.628,0.692,17)
PROBES=("psi","dphi","V")


def load_spec():
    p=ROOT/"tools/run_production_spectroscopy_suite.py"
    s=importlib.util.spec_from_file_location("specsuite",p)
    m=importlib.util.module_from_spec(s)
    assert s.loader is not None
    s.loader.exec_module(m)
    return m


def min_omega2(spec,K,G,S,M,k):
    o2,*_=spec.local_spectrum(K,G,S,M,k)
    return float(np.min(o2))


def threshold(spec,K,G,S,M,kunit):
    f0=min_omega2(spec,K,G,S,M,0.0)
    if f0>0:
        return 0.0
    hi=0.25
    while hi<=256.0 and min_omega2(spec,K,G,S,M,hi*kunit)<=0:
        hi*=2
    if hi>256.0:
        return None
    lo=hi/2
    for _ in range(60):
        mid=(lo+hi)/2
        if min_omega2(spec,K,G,S,M,mid*kunit)>0:
            hi=mid
        else:
            lo=mid
    return float(hi)


def main():
    spec=load_spec()
    b=build_onshell_central(ROOT)
    d=b.direct41.sort_values("x").reset_index(drop=True)
    prod=d[(d.u>0.62)&(d.u<0.70)].reset_index(drop=True)
    u=prod.u.to_numpy(float);x=prod.x.to_numpy(float)
    rows=spec.nearest_rows(u,TARGET_U)
    kunit=math.pi/max(abs(float(x[rows[-1]]-x[rows[0]])),1e-12)
    red=module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")

    result={}
    for L in LS:
        A=red.canonical_audit(prod,int(L))
        K,G,S,M=(np.asarray(A[q]) for q in ("K","G","S","M"))
        th=[]
        for i in rows:
            th.append(threshold(spec,K[i],G[i],S[i],M[i],kunit))
        finite=[z for z in th if z is not None]
        if len(finite)!=len(th):
            common=None
        else:
            common=1.20*max(finite)

        recs=[]
        prev=None
        frequencies=[]
        probe_hist={p:[] for p in PROBES}
        if common is not None:
            kval=common*kunit
            for ii,i in enumerate(rows):
                o2,U,V,hres=spec.local_spectrum(K[i],G[i],S[i],M[i],kval)
                if prev is not None:
                    p=spec.best_overlap_permutation(prev,U)
                    o2=o2[p];U=U[:,p];V=V[:,p]
                prev=U
                if not np.all(o2>0):
                    raise RuntimeError(f"L={L} common safe k failed at u={u[i]}")
                w=spec.normalized_probe_weights(V,o2)
                om=np.sqrt(o2)
                frequencies.append(om)
                for p in PROBES:
                    probe_hist[p].append(w[p])
                recs.append({
                    "u":float(u[i]),"x":float(x[i]),
                    "omega2":[float(z) for z in o2],
                    "omega":[float(z) for z in om],
                    "probe":w,
                    "hermiticity_residual":float(hres),
                })

        ratio_stats={}
        if frequencies:
            F=np.asarray(frequencies)
            for a in range(3):
                for bb in range(a+1,3):
                    rr=F[:,a]/F[:,bb]
                    ratio_stats[f"{a}/{bb}"]={
                        "mean":float(np.mean(rr)),
                        "fractional_range":float((np.max(rr)-np.min(rr))/max(abs(np.mean(rr)),1e-300)),
                    }

        psum={}
        for p in PROBES:
            h=probe_hist[p]
            if not h:
                psum[p]=None
                continue
            dom=[r["dominant_mode"] for r in h]
            psum[p]={
                "dominant_mode_switches":int(sum(dom[j]!=dom[j-1] for j in range(1,len(dom)))),
                "unique_dominant_modes":sorted(set(dom)),
                "mean_effective_mode_count":float(np.mean([r["effective_mode_count"] for r in h])),
                "min_effective_mode_count":float(np.min([r["effective_mode_count"] for r in h])),
                "max_top_fraction":float(np.max([r["top_fraction"] for r in h])),
            }

        result[str(L)]={
            "k_unit":float(kunit),
            "threshold_factor_by_radius":[
                {"u":float(u[i]),"threshold_factor":th[j]}
                for j,i in enumerate(rows)
            ],
            "max_threshold_factor":float(max(finite)) if finite else None,
            "common_safe_factor":float(common) if common is not None else None,
            "common_safe_k":float(common*kunit) if common is not None else None,
            "records_at_common_safe_k":recs,
            "frequency_ratio_stats":ratio_stats,
            "probe_summary":psum,
        }

    any_switch=any(
        rec["probe_summary"][p] is not None and rec["probe_summary"][p]["dominant_mode_switches"]>0
        for rec in result.values() for p in PROBES
    )
    max_ratio=max(
        (rr["fractional_range"] for rec in result.values() for rr in rec["frequency_ratio_stats"].values()),
        default=0.0,
    )
    min_neff=min(
        (rec["probe_summary"][p]["min_effective_mode_count"]
         for rec in result.values() for p in PROBES if rec["probe_summary"][p] is not None),
        default=None,
    )
    max_top=max(
        (rec["probe_summary"][p]["max_top_fraction"]
         for rec in result.values() for p in PROBES if rec["probe_summary"][p] is not None),
        default=None,
    )

    report={
        "scope":"healthy all-positive local-k spectroscopy; no artificial radial boundary; not QNM",
        "result":result,
        "summary":{
            "all_L_found_finite_all_positive_k":all(rec["common_safe_factor"] is not None for rec in result.values()),
            "any_probe_dominance_switch_at_common_safe_k":any_switch,
            "max_frequency_ratio_fractional_range_at_common_safe_k":float(max_ratio),
            "minimum_probe_effective_mode_count":float(min_neff) if min_neff is not None else None,
            "maximum_single_mode_probe_weight":float(max_top) if max_top is not None else None,
        },
        "interpretation_guard":(
            "Nonconstant local frequency ratios or probe weights reject a pure common-scaling "
            "description of this local frozen-k diagnostic, but they are not QNM residues Z_n(r)."
        ),
    }
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps(report["summary"],indent=2))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
