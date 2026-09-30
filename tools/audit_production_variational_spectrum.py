#!/usr/bin/env python3
"""Variational negative-mode witness for production-window spectroscopy.

The finite Dirichlet box spectra show negative omega^2 at every tested L, but
the positive part of the boxed spectrum is not yet well converged.  This audit
therefore avoids a matrix eigen-discretization claim and asks a sharper
question:

Does there exist an explicit smooth Dirichlet trial field Y(x)=q(x)v inside
the healthy production window whose Rayleigh quotient for the self-adjoint
reduced operator is negative?

For constant field direction v and scalar Dirichlet envelope q, the
antisymmetric S contribution vanishes identically (v^T S v = 0), so

  R[v] =
    v^T [∫(q'^2 G - q^2 M) dx] v
    / v^T [∫q^2 K dx] v.

A negative generalized eigenvalue of these integrated 3x3 matrices is an
explicit variational witness that the corresponding finite Dirichlet operator
has a negative omega^2 direction.  It is still NOT a QNM instability claim.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy import linalg

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))

from ssz_p5.config import DEFAULT_L  # noqa:E402
from ssz_p5.numerics import module  # noqa:E402
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central  # noqa:E402

OUT=ROOT/"data/generated/spectral/PRODUCTION_VARIATIONAL_NEGATIVE_MODE_WITNESS.json"
WINDOWS=((0.628,0.692),(0.632,0.688),(0.640,0.680))
HARMONICS=(1,2,3)
FIELDS=("psi","dphi","V")


def trapz_matrix(x,A):
    # numpy.trapezoid handles first axis
    return np.trapezoid(A,x,axis=0)


def witness_for(prod,A,lo,hi,n):
    u=prod.u.to_numpy(float)
    ids=np.flatnonzero((u>lo)&(u<hi))
    x=prod.x.to_numpy(float)[ids]
    K=np.asarray(A["K"],float)[ids]
    G=np.asarray(A["G"],float)[ids]
    M=np.asarray(A["M"],float)[ids]
    if x[0]>x[-1]:
        order=np.argsort(x);x=x[order];K=K[order];G=G[order];M=M[order]

    Lx=float(x[-1]-x[0])
    t=(x-x[0])/Lx
    q=np.sin(n*np.pi*t)
    qp=(n*np.pi/Lx)*np.cos(n*np.pi*t)
    Aint=trapz_matrix(x,qp[:,None,None]**2*((G+G.swapaxes(1,2))/2)
                     -q[:,None,None]**2*((M+M.swapaxes(1,2))/2))
    Bint=trapz_matrix(x,q[:,None,None]**2*((K+K.swapaxes(1,2))/2))
    vals,vecs=linalg.eigh((Aint+Aint.T)/2,(Bint+Bint.T)/2)
    j=int(np.argmin(vals))
    v=vecs[:,j]
    # Euclidean component power only for descriptive basis composition.
    p=np.abs(v/np.linalg.norm(v))**2
    return {
        "harmonic":int(n),
        "rayleigh_min":float(vals[j]),
        "all_generalized_eigenvalues":[float(z) for z in vals],
        "field_component_power":{FIELDS[i]:float(p[i]) for i in range(3)},
        "x_span":Lx,
        "rows":int(len(ids)),
    }


def main():
    b=build_onshell_central(ROOT)
    d=b.direct41.sort_values("x").reset_index(drop=True)
    prod=d[(d.u>0.62)&(d.u<0.70)].reset_index(drop=True)
    red=module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")

    result={}
    for L in DEFAULT_L:
        A=red.canonical_audit(prod,int(L))
        rows=[]
        for lo,hi in WINDOWS:
            for n in HARMONICS:
                r=witness_for(prod,A,lo,hi,n)
                r["u_window"]=[lo,hi]
                rows.append(r)
        best=min(rows,key=lambda z:z["rayleigh_min"])
        result[str(L)]={
            "trials":rows,
            "best":best,
            "negative_witness_exists":bool(best["rayleigh_min"]<0),
        }

    report={
        "scope":"explicit smooth finite-Dirichlet variational witness; not QNM",
        "result":result,
        "summary":{
            "negative_witness_for_all_L":all(v["negative_witness_exists"] for v in result.values()),
            "best_rayleigh_by_L":{L:v["best"]["rayleigh_min"] for L,v in result.items()},
            "best_field_composition_by_L":{L:v["best"]["field_component_power"] for L,v in result.items()},
        },
        "guard":(
            "A negative Rayleigh quotient certifies a negative direction of the "
            "finite self-adjoint Dirichlet diagnostic operator. It does not by "
            "itself establish a physical center-regular/outgoing QNM instability."
        ),
    }
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps(report["summary"],indent=2))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
