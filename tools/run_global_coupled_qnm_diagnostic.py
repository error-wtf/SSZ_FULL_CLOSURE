#!/usr/bin/env python3
"""Assemble a diagnostic global regional 41-slot stream and run a global QNM attempt.

IMPORTANT SCIENTIFIC SCOPE
--------------------------
This is an explicit *diagnostic* execution of the currently selected regional
architecture.  It intentionally uses the locked/prestaged regional 41-slot
representatives according to production/regions.py so that the full radial
operator can actually be exercised end-to-end now.

It is NOT a Direct-Global-KRGSM production certificate, because some regional
tables are prestaged/reference representatives rather than freshly regenerated
from one covariant action.  The report states this explicitly.

The QNM layer uses the full reduced Euler operator P(dt,dr) from the accepted
profile-aware reducer.  A finite-domain polynomial eigenproblem is built with
regular inner and outgoing Sommerfeld outer boundary rows.  Multiple outer
cutoffs/resolutions are compared.  Any roots that move are rejected.  Any
persistent Im(omega)>0 branch is reported as a diagnostic instability candidate,
not as a certified physical QNM until the direct-global action provenance and
full stability gates pass.
"""
from __future__ import annotations

import json, hashlib, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.linalg import eig

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))

from ssz_p5.config import SLOT_NAMES
from ssz_p5.numerics import module
from ssz_p5.production.regions import (
    OUTER_START,CENTRAL_START,INNER_START,CORE_START,classify_u
)
from ssz_p5.production.regional_coefficients import select_lower,central_selected
from ssz_p5.production.sources import SOURCE_REGISTRY
from ssz_p5.qnm.descriptor_pencil import semidiscrete_polynomial

OUTDIR=ROOT/"data/generated/qnm_global_diagnostic"
STREAM=OUTDIR/"GLOBAL_REGIONAL_REFERENCE_41.csv"
REPORT=OUTDIR/"GLOBAL_QNM_DIAGNOSTIC_REPORT.json"

LS=(6,12,20,42)
RESOLUTIONS=(45,65,85)
OUTER_XMAX=(18.0,25.0,35.0)


def sha256(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(1<<20),b""): h.update(b)
    return h.hexdigest()


def normalize_reference(path,region):
    d=pd.read_csv(path).copy()
    if "u" not in d.columns and "x" in d.columns:
        d["u"]=1.0/d["x"].to_numpy(float)
    if "x" not in d.columns and "u" in d.columns:
        d["x"]=1.0/d["u"].to_numpy(float)
    if "phiprime" not in d.columns and "phi_r" in d.columns:
        d["phiprime"]=d["phi_r"]
    if "A0prime" not in d.columns:
        d["A0prime"]=0.0
    # Pure-H weak exterior is the documented 39/41 table; missing lower controls
    # are zero in that pure-H representative.
    for c in ("v5","c3","e3"):
        if c not in d.columns or not np.isfinite(pd.to_numeric(d[c], errors="coerce")).all():
            d[c]=0.0
    d=select_lower(d)
    d["production_region"]=region
    return d


def assemble():
    reg=SOURCE_REGISTRY
    weak=normalize_reference(ROOT/reg["weak_exterior_H"]["coeff_reference"],"weak_exterior_H")
    outer_current=ROOT/"data/generated/qnm_global_diagnostic/OUTER_COMPLETE_H_SVT_41.csv"
    if outer_current.exists():
        outer=normalize_reference(outer_current,"outer_same_action_H_SVT")
    else:
        outer=normalize_reference(ROOT/reg["outer_same_action_H_SVT"]["coeff_reference"],"outer_same_action_H_SVT")
    central=central_selected(ROOT)
    central["production_region"]="central_exact_SVT"
    inner=normalize_reference(ROOT/reg["inner_same_action_SVT_H"]["coeff_reference"],"inner_same_action_SVT_H")
    core=normalize_reference(ROOT/reg["punctured_H_core"]["coeff_reference"],"punctured_H_core")

    pieces=[]
    for d,name,lo,hi in (
      (weak,"weak_exterior_H",None,OUTER_START),
      (outer,"outer_same_action_H_SVT",OUTER_START,CENTRAL_START),
      (central,"central_exact_SVT",CENTRAL_START,INNER_START),
      (inner,"inner_same_action_SVT_H",INNER_START,CORE_START),
      (core,"punctured_H_core",CORE_START,None),
    ):
        u=d.u.to_numpy(float)
        m=np.ones(len(d),bool)
        if lo is not None:m &= u>=lo
        if hi is not None:m &= u<hi
        q=d.loc[m].copy()
        if len(q)==0:
            raise RuntimeError(f"empty selected region {name}")
        pieces.append(q)

    g=pd.concat(pieces,ignore_index=True)
    g=g.sort_values("x").drop_duplicates("x",keep="last").reset_index(drop=True)
    if np.any(np.diff(g.x.to_numpy(float))<=0): raise RuntimeError("global x not increasing")
    slots=list(SLOT_NAMES)
    if not np.isfinite(g[slots].to_numpy(float)).all(): raise RuntimeError("nonfinite slots")
    return g,pieces


def interface_report(g):
    out={}
    for u0 in (OUTER_START,CENTRAL_START,INNER_START,CORE_START):
        x0=1/u0
        ix=np.argsort(np.abs(g.x.to_numpy(float)-x0))[:4]
        out[str(u0)]=[{
          "x":float(g.x.iloc[i]),"u":float(g.u.iloc[i]),
          "region":str(g.production_region.iloc[i])
        } for i in ix]
    return out


def qnm_pencil(P,x,indices,L):
    """Finite-domain quadratic pencil with simple regular/outgoing boundary rows.

    Lambda convention exp(lambda t), omega=-i lambda.  Outer outgoing condition
    for asymptotically unit speed: y_r + lambda y = 0.  Inner diagnostic
    regularity: y_r = 0.  Applied componentwise.
    """
    A0,A1,A2,xs=semidiscrete_polynomial(P,x,indices=indices,window=9,degree=8)
    A0=A0.toarray().astype(complex)
    A1=A1.toarray().astype(complex)
    A2=A2.toarray().astype(complex)
    nf=3;n=len(xs)

    # derivative matrix on selected nodes
    qmod=__import__("ssz_p5.qnm.descriptor_pencil",fromlist=["local_poly_differentiation_matrix"])
    D=qmod.local_poly_differentiation_matrix(xs,1,9,8).toarray()

    # Replace PDE rows at boundaries componentwise.
    for f in range(nf):
        rin=f*n
        rout=f*n+n-1
        A0[rin,:]=0;A1[rin,:]=0;A2[rin,:]=0
        A0[rin,f*n:(f+1)*n]=D[0]
        A0[rout,:]=0;A1[rout,:]=0;A2[rout,:]=0
        A0[rout,f*n:(f+1)*n]=D[-1]
        A1[rout,rout]=1.0

    N=A0.shape[0]
    I=np.eye(N,dtype=complex);Z=np.zeros_like(I)
    Lm=np.block([[-A1,-A0],[I,Z]])
    Rm=np.block([[A2,Z],[Z,I]])
    vals=eig(Lm,Rm,right=False,check_finite=False)
    vals=vals[np.isfinite(vals.real)&np.isfinite(vals.imag)]
    # remove enormous descriptor/infinite artifacts
    vals=vals[np.abs(vals)<50]
    # omega = i lambda for exp(-i omega t) if lambda=-i omega => omega=i lambda
    omega=1j*vals
    return omega,xs


def match_roots(rootsets,tol=0.12):
    """Roots from first set that have a neighbor in every other set."""
    base=rootsets[0]
    stable=[]
    for z in base:
        ds=[]
        ok=True
        for rr in rootsets[1:]:
            if len(rr)==0:ok=False;break
            d=np.min(np.abs(rr-z))
            ds.append(float(d))
            if d>tol:ok=False;break
        if ok:
            stable.append((z,max(ds) if ds else 0.0))
    stable.sort(key=lambda q:(abs(q[0].imag),abs(q[0].real)))
    return stable


def main():
    g,pieces=assemble()
    OUTDIR.mkdir(parents=True,exist_ok=True)
    g.to_csv(STREAM,index=False)

    # Exact direct-stream pivot gate before any reduction/QNM solve.  Some
    # prestaged regional files are sector/reference products rather than one
    # complete assembled same-action member; zero auxiliary pivots make that
    # distinction operationally visible and must block the coupled solver.
    pivot_detail = {}
    pivot_block = False
    for region, q in g.groupby("production_region", sort=False):
        rec = {}
        for name in ("v1","v9","v10"):
            if name in q:
                a=np.asarray(q[name],float)
                rec[name]={
                    "min_abs":float(np.min(np.abs(a))),
                    "zero_or_nearzero_rows":int(np.sum(np.abs(a)<1e-14)),
                }
                pivot_block |= rec[name]["zero_or_nearzero_rows"]>0
        if all(c in q for c in ("b1","v10","v11")):
            DeltaV=4*np.asarray(q.b1,float)*np.asarray(q.v10,float)-np.asarray(q.v11,float)**2
            rec["DeltaV"]={
                "min_abs":float(np.min(np.abs(DeltaV))),
                "zero_or_nearzero_rows":int(np.sum(np.abs(DeltaV)<1e-14)),
            }
            pivot_block |= rec["DeltaV"]["zero_or_nearzero_rows"]>0
        pivot_detail[str(region)]=rec

    if pivot_block:
        report={
          "status":"GLOBAL_COUPLED_QNM_BLOCKED_AT_DIRECT_STREAM_PIVOTS",
          "scientific_scope":"global coupled HSVT execution attempted; stopped at mandatory direct-stream constraint-pivot gate",
          "stream_sha256":sha256(STREAM),
          "rows":len(g),
          "x_range":[float(g.x.min()),float(g.x.max())],
          "u_range":[float(g.u.min()),float(g.u.max())],
          "interfaces":interface_report(g),
          "pivot_gate":pivot_detail,
          "physical_qnm_claim_allowed":False,
          "blockers":[
            "selected regional reference assembly contains zero/near-zero auxiliary pivots before reduction",
            "the offending prestaged/raw sector is not a complete direct same-action H+SVT stream",
            "QNM eigenvalues from such a singularly incomplete assembled operator would be nonphysical"
          ],
          "next_target":"REGENERATE_COMPLETE_OUTER_H_PLUS_SVT_DIRECT_41_FROM_ACTION_AND_SHARED_BASELINE_THEN_RERUN_GLOBAL_KRGSM_QNM"
        }
        REPORT.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
        print(json.dumps(report,indent=2,allow_nan=False))
        return 0

    reducer=module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    x=g.x.to_numpy(float)

    health={}
    spectra={}
    for L in LS:
        a=reducer.canonical_audit(g,float(L))
        K=np.asarray(a["K"],float);G=np.asarray(a["G"],float)
        ke=np.linalg.eigvalsh((K+K.transpose(0,2,1))/2)[:,0]
        cr=np.full(len(g),np.nan)
        for i in np.flatnonzero(ke>0):
            w,U=np.linalg.eigh((K[i]+K[i].T)/2)
            inv=U@np.diag(1/np.sqrt(w))@U.T
            C=inv@((G[i]+G[i].T)/2)@inv
            cr[i]=np.linalg.eigvalsh((C+C.T)/2)[0]
        health[str(L)]={
          "min_K":float(np.min(ke)),
          "min_cr2_where_K_positive":float(np.nanmin(cr)),
          "negative_K_rows":int(np.sum(ke<=0)),
          "negative_cr2_rows_where_defined":int(np.sum(cr<0)),
        }

        rootsets=[]
        configs=[]
        for N,xmax in zip(RESOLUTIONS,OUTER_XMAX):
            # Keep common inner edge; vary outer edge and radial count.
            ids=np.flatnonzero(x<=xmax)
            if len(ids)<N: Nuse=len(ids)
            else:Nuse=N
            choose=np.unique(np.linspace(ids[0],ids[-1],Nuse).round().astype(int))
            # Need at least 9 nodes.
            if len(choose)<9: continue
            om,xs=qnm_pencil(a["P"],x,choose,L)
            # physically interesting bounded window
            om=om[(np.abs(om.real)<5)&(np.abs(om.imag)<5)]
            rootsets.append(om)
            configs.append({"N":len(choose),"x_min":float(xs.min()),"x_max":float(xs.max()),"root_count":len(om)})
        stable=match_roots(rootsets)
        spectra[str(L)]={
          "configs":configs,
          "stable_root_count":len(stable),
          "stable_roots":[{"omega_re":float(z.real),"omega_im":float(z.imag),"max_cross_config_delta":d}
                          for z,d in stable[:40]],
          "stable_growing_count":int(sum(z.imag>1e-6 for z,d in stable)),
          "stable_damped_count":int(sum(z.imag<-1e-6 for z,d in stable)),
        }

    full_health=all(v["min_K"]>0 and v["min_cr2_where_K_positive"]>0 for v in health.values())
    any_stable=any(v["stable_root_count"]>0 for v in spectra.values())
    status=(
      "GLOBAL_QNM_DIAGNOSTIC_EXECUTED_UNCERTIFIED_MEMBER"
      if any_stable else
      "GLOBAL_QNM_DIAGNOSTIC_EXECUTED_NO_CONVERGED_ROOTS"
    )
    report={
      "status":status,
      "scientific_scope":"global regional-reference diagnostic; NOT a direct-global production certificate",
      "stream_sha256":sha256(STREAM),
      "rows":len(g),"x_range":[float(x.min()),float(x.max())],
      "u_range":[float(g.u.min()),float(g.u.max())],
      "interfaces":interface_report(g),
      "finite_l_health":health,
      "full_health_pass":full_health,
      "qnm_boundary_model":{
        "inner":"componentwise Neumann regularity on innermost numerical node",
        "outer":"componentwise Sommerfeld y_r + lambda y = 0",
        "warning":"finite-domain diagnostic only; compactified/Jost/ECS production boundary not yet implemented here"
      },
      "spectra":spectra,
      "physical_qnm_claim_allowed":False,
      "blockers":[
        "reference/prestaged regional coefficients are not one fresh direct covariant-action export",
        "current selected regional member fails full-domain K/cr2 health if any finite_l_health entry is negative",
        "finite-radius Sommerfeld diagnostic is not the accepted production resonance boundary method"
      ],
      "next_target":"IF_ROOTS_STABLE_BUILD_COMPACTIFIED_JOST_OR_ECS_ON_SAME_OPERATOR; OTHERWISE_LOCALIZE_BOUNDARY_SENSITIVITY"
    }
    REPORT.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps({
      "status":status,
      "full_health_pass":full_health,
      "health":health,
      "spectra_summary":{L:{
        "stable":v["stable_root_count"],
        "growing":v["stable_growing_count"],
        "damped":v["stable_damped_count"]
      } for L,v in spectra.items()},
      "physical_qnm_claim_allowed":False
    },indent=2))
    return 0

if __name__=="__main__": raise SystemExit(main())
