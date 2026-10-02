#!/usr/bin/env python3
"""Numerical local SSZ resolvent spectra on the locked current member.

This evaluates the actual current-member K,G matrices at selected radii inside
0.62<u<0.70 and computes:
  G v_n = omega_n^2 K v_n,
  Z_n(q)=|q^T v_n|^2/(2 omega_n),
  A_q(omega)=|Im q^T [G-(omega+i eta)^2 K]^-1 q|/pi.

The output contains raw mode frequencies/residues and densely sampled spectral
curves.  These are local-principal spectra, not global QNMs.
"""
from __future__ import annotations

import json,sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.linalg import eigh
from scipy.signal import find_peaks

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))

from ssz_p5.numerics import module  # noqa:E402
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central  # noqa:E402

OUTDIR=ROOT/"data/generated/spectral/ssz_local_resolvent"
MODES=OUTDIR/"SSZ_LOCAL_RESOLVENT_MODES.csv"
CURVES=OUTDIR/"SSZ_LOCAL_RESOLVENT_CURVES.csv"
REPORT=OUTDIR/"SSZ_LOCAL_RESOLVENT_REPORT.json"

LS=(6,20,42,110,420,1000)
U_TARGETS=(0.63,0.65,0.67,0.69)
ETA_FRACS=(0.01,0.005,0.002)


def sym(A): return (A+A.T)/2


def solve(K,G):
    K=sym(K);G=sym(G)
    kw=np.linalg.eigvalsh(K)
    if np.min(kw)<=0: raise ValueError("K_NOT_POSITIVE")
    lam,V=eigh(G,K,check_finite=True)
    if np.min(lam)<=0: raise ValueError("G_GENERALIZED_NOT_POSITIVE")
    return lam,np.sqrt(lam),V


def observables(n):
    out={f"field_{i}":np.eye(n)[i] for i in range(n)}
    out["equal_weight"]=np.ones(n)/np.sqrt(n)
    q=np.arange(1,n+1,dtype=float)
    out["ramp_weight"]=q/np.linalg.norm(q)
    return out


def residues(omega,V,q):
    raw=np.abs(q@V)**2/np.maximum(2*omega,1e-300)
    norm=raw/max(float(raw.sum()),1e-300)
    return raw,norm


def spectral_curve(K,G,q,omega,eta_frac):
    lo=max(0.0,float(np.min(omega))*0.65)
    hi=float(np.max(omega))*1.35
    x=np.linspace(lo,hi,5000)
    eta=max(eta_frac*float(np.median(omega)),1e-12)
    A=np.empty_like(x)
    for i,w in enumerate(x):
        z=complex(float(w),eta)
        val=q@np.linalg.solve(G.astype(complex)-z*z*K.astype(complex),q)
        A[i]=abs(np.imag(val))/np.pi
    return x,A,eta


def main():
    b=build_onshell_central(ROOT)
    d=b.direct41.sort_values("x").reset_index(drop=True)
    u=d.u.to_numpy(float)
    red=module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")

    mode_rows=[];curve_rows=[];cases=[]
    for L in LS:
        a=red.canonical_audit(d,int(L))
        Kall=np.asarray(a["K"],float)
        Gall=np.asarray(a["G"],float)
        for ut in U_TARGETS:
            i=int(np.argmin(np.abs(u-ut)))
            if not (0.62<u[i]<0.70): raise RuntimeError("target outside production")
            K=sym(Kall[i]);G=sym(Gall[i])
            lam,om,V=solve(K,G)
            obs=observables(len(om))
            for oname,q in obs.items():
                raw,norm=residues(om,V,q)
                order=np.argsort(-norm)
                for n in range(len(om)):
                    mode_rows.append({
                      "L":L,"u":float(u[i]),"observable":oname,"mode":n,
                      "omega2":float(lam[n]),"omega":float(om[n]),
                      "raw_residue":float(raw[n]),"normalized_residue":float(norm[n]),
                      "rank":int(np.where(order==n)[0][0]+1),
                    })
                ec=float(1/np.sum(norm*norm))
                ent=float(-np.sum(np.where(norm>0,norm*np.log(norm),0)))
                for ef in ETA_FRACS:
                    x,A,eta=spectral_curve(K,G,q,om,ef)
                    peaks,_=find_peaks(A)
                    pidx=peaks[np.argsort(A[peaks])[-min(len(peaks),len(om)):]] if len(peaks) else np.array([],int)
                    for xx,yy in zip(x,A):
                        curve_rows.append({
                          "L":L,"u":float(u[i]),"observable":oname,
                          "eta_frac":ef,"omega":float(xx),"spectral_density":float(yy),
                        })
                    cases.append({
                      "L":L,"u":float(u[i]),"observable":oname,"eta_frac":ef,
                      "eta":eta,"effective_mode_count":ec,"entropy":ent,
                      "top1":float(np.max(norm)),
                      "top2":float(np.sum(np.sort(norm)[-min(2,len(norm)):])),
                      "dominant_mode":int(order[0]),
                      "dominant_omega":float(om[order[0]]),
                      "peak_omegas":[float(x[k]) for k in pidx[np.argsort(x[pidx])]] if len(pidx) else [],
                      "peak_heights":[float(A[k]) for k in pidx[np.argsort(x[pidx])]] if len(pidx) else [],
                    })

    OUTDIR.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(mode_rows).to_csv(MODES,index=False)
    pd.DataFrame(curve_rows).to_csv(CURVES,index=False)

    # Robustness: residue weights themselves do not depend on eta; peak positions
    # should stabilize as eta decreases.
    grouped={}
    for c in cases:
        key=f"L{c['L']}_u{c['u']:.6f}_{c['observable']}"
        grouped.setdefault(key,[]).append(c)
    robust={}
    for key,vals in grouped.items():
        vals=sorted(vals,key=lambda x:x["eta_frac"])
        robust[key]={
          "top1":vals[0]["top1"],"top2":vals[0]["top2"],
          "effective_mode_count":vals[0]["effective_mode_count"],
          "dominant_mode":vals[0]["dominant_mode"],
          "dominant_omega":vals[0]["dominant_omega"],
          "eta_peak_counts":{str(v["eta_frac"]):len(v["peak_omegas"]) for v in vals},
        }

    report={
      "status":"SSZ_LOCAL_PRINCIPAL_RESOLVENT_SPECTRA_EVALUATED",
      "scope":"current locked member, registered healthy window only",
      "L_values":list(LS),"u_targets":list(U_TARGETS),"eta_fracs":list(ETA_FRACS),
      "case_count":len(cases),"summary":robust,
      "artifacts":{
        "modes_csv":str(MODES.relative_to(ROOT)),
        "curves_csv":str(CURVES.relative_to(ROOT)),
      },
      "guard":"These are local principal generalized-mode/resolvent spectra, not center-to-infinity QNM spectra."
    }
    REPORT.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps({
      "status":report["status"],"case_count":len(cases),
      "top1_range":[min(c["top1"] for c in cases),max(c["top1"] for c in cases)],
      "top2_range":[min(c["top2"] for c in cases),max(c["top2"] for c in cases)],
      "neff_range":[min(c["effective_mode_count"] for c in cases),max(c["effective_mode_count"] for c in cases)],
    },indent=2))
    return 0

if __name__=="__main__": raise SystemExit(main())
