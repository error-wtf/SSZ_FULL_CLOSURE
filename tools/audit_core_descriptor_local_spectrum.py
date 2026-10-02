#!/usr/bin/env python3
"""Deep-core descriptor local-dispersion audit.

The historical 3-field Schur reducer is explicitly superseded in the punctured
core because the D_h1 elimination is ill-conditioned.  This audit keeps the
8-field generalized-psi descriptor unreduced and evaluates its local polynomial
symbol for real radial wavenumbers.

For each sampled core radius, multipole L and real k_r, solve the quadratic
descriptor eigenproblem in lambda (exp(lambda t + i k_r r)). Infinite DAE
eigenvalues are discarded.  A conservative locally stable branch should have
Re(lambda) ~ 0; positive Re(lambda) indicates a local growing descriptor mode.

This is a local principal/dispersion test, not yet the global QNM boundary-value
problem.
"""
from __future__ import annotations
import json,sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.linalg import eig

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from ssz_p5.reducer.unreduced_descriptor import generalized_psi_descriptor  # noqa:E402

CORE=ROOT/"data/prestaged/direct41/core_selected_41of41.csv"
OUT=ROOT/"data/generated/qnm_global_diagnostic/CORE_DESCRIPTOR_LOCAL_SPECTRUM.json"
LS=(6,12,20,42,110,420,1000)
KRS=(0.0,0.25,0.5,1.0,2.0,4.0)
NSAMP=36

def local_poly(P,i,kr):
    # Q(lambda)=A2 l^2 + A1 l + A0 at one radius with dr -> i kr.
    sample=next(iter(P.values()));nf=sample.shape[1]
    A=[np.zeros((nf,nf),complex) for _ in range(3)]
    for (pt,pr),B in P.items():
        A[pt]+=np.asarray(B[i],complex)*(1j*kr)**pr
    return A

def roots(A0,A1,A2):
    n=A0.shape[0];I=np.eye(n,dtype=complex);Z=np.zeros_like(I)
    Lm=np.block([[-A1,-A0],[I,Z]])
    Rm=np.block([[A2,Z],[Z,I]])
    z=eig(Lm,Rm,right=False,check_finite=False)
    z=z[np.isfinite(z.real)&np.isfinite(z.imag)]
    z=z[np.abs(z)<1e8]
    return z

def main():
    d=pd.read_csv(CORE).sort_values("x").reset_index(drop=True)
    # Strict production core only.
    if "u" not in d:d["u"]=1/d.x
    d=d[d.u>=0.715].reset_index(drop=True)
    ids=np.unique(np.linspace(4,len(d)-5,min(NSAMP,max(1,len(d)-8))).round().astype(int))
    results={}
    global_max=-np.inf;global_worst=None
    for L in LS:
        P=generalized_psi_descriptor(d,float(L))
        entries=[]
        for i in ids:
            for kr in KRS:
                A0,A1,A2=local_poly(P,int(i),kr)
                z=roots(A0,A1,A2)
                if len(z)==0:
                    entries.append({"row":int(i),"u":float(d.u.iloc[i]),"x":float(d.x.iloc[i]),"kr":kr,"finite_roots":0})
                    continue
                mx=float(np.max(z.real));mn=float(np.min(z.real))
                sym=float(np.max(np.abs(np.sort(z.real)+np.sort(z.real)[::-1]))) if len(z)>1 else abs(mx+mn)
                ent={
                  "row":int(i),"u":float(d.u.iloc[i]),"x":float(d.x.iloc[i]),"kr":kr,
                  "finite_roots":int(len(z)),"max_Re_lambda":mx,"min_Re_lambda":mn,
                  "max_abs_Im_lambda":float(np.max(np.abs(z.imag))),
                  "growth_pair_symmetry_proxy":sym,
                }
                entries.append(ent)
                if mx>global_max:
                    global_max=mx;global_worst={"L":L,**ent}
        results[str(L)]={
          "max_Re_lambda":max((e.get("max_Re_lambda",-np.inf) for e in entries),default=None),
          "sample_count":len(entries),
          "entries":entries,
        }
    # Descriptor Hamiltonian/conservative local symbols commonly yield +/- pairs;
    # only call exact local pass if real parts are at numerical scale.
    tol=1e-7
    status="CORE_DESCRIPTOR_LOCAL_NO_GROWTH" if global_max<tol else "CORE_DESCRIPTOR_LOCAL_GROWTH_PRESENT"
    report={
      "status":status,
      "scope":"unreduced generalized-psi deep-core local dispersion; no D_h1 Schur inverse",
      "lambda_convention":"exp(lambda t + i k_r r)",
      "growth_tolerance":tol,
      "global_max_Re_lambda":float(global_max),
      "worst":global_worst,
      "per_L":results,
      "guard":"Local descriptor dispersion only; not a center-to-infinity QNM spectrum."
    }
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps({k:report[k] for k in ("status","global_max_Re_lambda","worst","guard")},indent=2))
    return 0

if __name__=="__main__":raise SystemExit(main())
