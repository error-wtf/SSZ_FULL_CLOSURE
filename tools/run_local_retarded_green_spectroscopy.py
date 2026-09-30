#!/usr/bin/env python3
"""Retarded local Green-function spectroscopy for the healthy production window.

This turns the healthy all-positive frozen-k dispersion diagnostic into an
explicit spectral function.

At each L we find a single common radial momentum k for which all three local
branches are positive throughout 0.628<=u<=0.692, track mode identity by
K-whitened eigenvector overlap, and construct

    G^R_loc(r;k,omega)
      = [(omega+i eta)^2 K(r) - H(k,r)]^{-1},

with
    H(k,r)=k^2 G + i k S - M.

For canonical basis probe e_a the positive-frequency modal residue is

    Z_n^(a)(r) = |e_a^T v_n(r)|^2 / (2 omega_n(r))

for K-normalized eigenvectors v_n.  We report exact residue weights and check
that a direct broadened spectral function reproduces the dominant line.

This is the closest admissible spectroscopy test before a direct-global KRGM
certificate exists, but it remains a *local frozen-k* Green function rather
than the physical center-to-infinity G^R(r,r';omega).
"""

from __future__ import annotations

import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))

from ssz_p5.config import DEFAULT_L  # noqa:E402
from ssz_p5.numerics import module  # noqa:E402
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central  # noqa:E402

OUT=ROOT/"data/generated/spectral/LOCAL_RETARDED_GREEN_SPECTROSCOPY.json"
TARGET_U=np.linspace(0.628,0.692,17)
PROBES=("psi","dphi","V")


def load_spec():
    p=ROOT/"tools/run_production_spectroscopy_suite.py"
    s=importlib.util.spec_from_file_location("specsuite",p)
    m=importlib.util.module_from_spec(s)
    assert s.loader is not None
    s.loader.exec_module(m)
    return m


def min_o2(spec,K,G,S,M,k):
    return float(np.min(spec.local_spectrum(K,G,S,M,k)[0]))


def threshold(spec,K,G,S,M,kunit):
    if min_o2(spec,K,G,S,M,0)>0:
        return 0.0
    hi=0.25
    while hi<=256 and min_o2(spec,K,G,S,M,hi*kunit)<=0:
        hi*=2
    if hi>256:
        return None
    lo=hi/2
    for _ in range(56):
        mid=(lo+hi)/2
        if min_o2(spec,K,G,S,M,mid*kunit)>0: hi=mid
        else: lo=mid
    return float(hi)


def direct_rho(K,H,omega,eta,a):
    e=np.zeros(3,dtype=complex);e[a]=1
    vals=[]
    for w in omega:
        A=((w+1j*eta)**2)*K-H
        g=np.vdot(e,np.linalg.solve(A,e))
        vals.append(float(-np.imag(g)/np.pi))
    return np.asarray(vals)


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
    for L in DEFAULT_L:
        A=red.canonical_audit(prod,int(L))
        K,G,S,M=(np.asarray(A[q]) for q in ("K","G","S","M"))
        th=[threshold(spec,K[i],G[i],S[i],M[i],kunit) for i in rows]
        if any(z is None for z in th):
            raise RuntimeError(f"L={L}: no finite positive-k threshold")
        fac=1.20*max(th); kval=fac*kunit

        prev=None
        recs=[]
        for i in rows:
            o2,U,V,hres=spec.local_spectrum(K[i],G[i],S[i],M[i],kval)
            if prev is not None:
                p=spec.best_overlap_permutation(prev,U)
                o2=o2[p];U=U[:,p];V=V[:,p]
            prev=U
            if not np.all(o2>0):
                raise RuntimeError(f"L={L}: common k not positive at u={u[i]}")
            om=np.sqrt(o2)
            weights=spec.normalized_probe_weights(V,o2)

            # direct broadened spectral validation over a grid spanning all lines
            eta=max(1e-6,0.01*float(np.min(om)))
            wgrid=np.linspace(max(0.0,0.5*float(np.min(om))),1.25*float(np.max(om)),1800)
            H=kval*kval*((G[i]+G[i].T)/2)+1j*kval*((S[i]-S[i].T)/2)-((M[i]+M[i].T)/2)
            probe_spectra={}
            for a,name in enumerate(PROBES):
                rho=direct_rho(K[i],H,wgrid,eta,a)
                pk,_=find_peaks(rho)
                if len(pk):
                    order=pk[np.argsort(rho[pk])[::-1]]
                    peaks=[
                        {"omega":float(wgrid[j]),"rho":float(rho[j])}
                        for j in order[:3]
                    ]
                    dominant=float(wgrid[order[0]])
                else:
                    peaks=[];dominant=None
                probe_spectra[name]={
                    "exact_residue_weights":weights[name]["weights"],
                    "exact_dominant_mode":weights[name]["dominant_mode"],
                    "exact_top_fraction":weights[name]["top_fraction"],
                    "broadened_dominant_peak_omega":dominant,
                    "broadened_top_peaks":peaks,
                }

            recs.append({
                "u":float(u[i]),"x":float(x[i]),
                "omega":[float(z) for z in om],
                "probe_spectra":probe_spectra,
                "eta":float(eta),
                "hermiticity_residual":float(hres),
            })

        summary={}
        for name in PROBES:
            dom=[r["probe_spectra"][name]["exact_dominant_mode"] for r in recs]
            switches=[]
            for j in range(1,len(dom)):
                if dom[j]!=dom[j-1]:
                    switches.append({
                        "u_left":recs[j-1]["u"],"u_right":recs[j]["u"],
                        "mode_left":dom[j-1],"mode_right":dom[j],
                        "weights_left":recs[j-1]["probe_spectra"][name]["exact_residue_weights"],
                        "weights_right":recs[j]["probe_spectra"][name]["exact_residue_weights"],
                    })
            summary[name]={
                "dominance_switch_count":len(switches),
                "switches":switches,
                "max_top_fraction":float(max(r["probe_spectra"][name]["exact_top_fraction"] for r in recs)),
                "min_top_fraction":float(min(r["probe_spectra"][name]["exact_top_fraction"] for r in recs)),
            }

        F=np.asarray([r["omega"] for r in recs])
        ratios={}
        for a in range(3):
            for bb in range(a+1,3):
                q=F[:,a]/F[:,bb]
                ratios[f"{a}/{bb}"]={
                    "fractional_range":float((q.max()-q.min())/max(abs(q.mean()),1e-300)),
                    "min":float(q.min()),"max":float(q.max()),
                }

        result[str(L)]={
            "common_safe_factor":float(fac),
            "common_safe_k":float(kval),
            "threshold_factor_max":float(max(th)),
            "records":recs,
            "probe_summary":summary,
            "frequency_ratio_stats":ratios,
        }

    report={
        "scope":"explicit local retarded Green-function spectroscopy at all-positive frozen k; not global QNM",
        "result":result,
        "summary":{
            "total_probe_dominance_switches":int(sum(
                r["probe_summary"][p]["dominance_switch_count"]
                for r in result.values() for p in PROBES
            )),
            "any_residue_ranking_reversal":bool(any(
                r["probe_summary"][p]["dominance_switch_count"]>0
                for r in result.values() for p in PROBES
            )),
            "max_frequency_ratio_fractional_range":float(max(
                q["fractional_range"]
                for r in result.values() for q in r["frequency_ratio_stats"].values()
            )),
        },
        "guard":(
            "These residues belong to the local frozen-k resolvent. They are not "
            "the center-regular/outgoing-infinity QNM residues Z_n(r) of the "
            "full global Green function."
        ),
    }
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps(report["summary"],indent=2))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
