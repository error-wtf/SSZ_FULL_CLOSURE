#!/usr/bin/env python3
"""Regenerate the quartic G5=0 core 41-stream directly from action jets.

Domain: punctured core from u>=0.715 up to the onset of the deep full-G5
subcore (u<50).  This uses:
  * quartic G3/G4 action jets -> F,G,H,a1,c4 via mh_action_primitives;
  * the documented background-null G2,XX normal direction to realize c2;
  * mh_general_primitives for a fresh 41-slot emission.

No historical selected 41 coefficients are used to construct the emitted
stream.  The old F2 table is used only for the c2 target that the transverse
G2,XX normal control is designed to realize, and as a regression oracle.

This is the direct action-level repair demanded by the core implementation
audit for the quartic (G5=0) part of the core.  The deep u>=50 full-G5 subcore
is explicitly excluded and remains a separate bridge.
"""
from __future__ import annotations

import json, hashlib, sys
from pathlib import Path
import numpy as np, pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))

from ssz_p5.coefficients.mh_action_primitives import quartic_g5zero_primitives
from ssz_p5.coefficients.mh_general_primitives import emit_from_primitives
from ssz_p5.numerics import module

ACT=ROOT/"data/production/ssz_p5_F1b_G4XX_transverse_core_candidate_2026-09-14.csv"
REF=ROOT/"data/authoritative/ssz_p5_F2_core_punctured_horndeski_unreduced_39of41_2026-09-14.csv"
OUTDIR=ROOT/"data/generated/qnm_global_diagnostic"
OUT=OUTDIR/"CORE_QUARTIC_G5ZERO_DIRECT41_ACTION_REGEN.csv"
REPORT=OUTDIR/"CORE_QUARTIC_G5ZERO_ACTION_REGEN_AUDIT.json"
LS=(6,12,20,42,110,420,1000)

def sha(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):h.update(b)
    return h.hexdigest()

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

def main():
    act=pd.read_csv(ACT).sort_values("r_over_rs").reset_index(drop=True)
    ref=pd.read_csv(REF).sort_values("x").reset_index(drop=True)
    q=quartic_g5zero_primitives(act)
    if len(q)!=len(ref) or np.max(np.abs(q.x.to_numpy(float)-ref.x.to_numpy(float)))>1e-12:
        raise RuntimeError("action/reference grid mismatch")

    # Realize the chosen c2 profile with the exact background-null G2XX normal
    # coefficient.  This leaves background action values/first jets unchanged.
    dc=q.dc2_dG2XX_action.to_numpy(float)
    c2_target=ref.c2.to_numpy(float)
    domain=(act.u.to_numpy(float)>=0.715)&(act.u.to_numpy(float)<50.0)
    bad=np.abs(dc)<1e-16
    unresolved = domain & bad & (np.abs(c2_target)>1e-12)
    G2XX=np.zeros_like(dc)
    ok=domain & (~bad)
    G2XX[ok]=c2_target[ok]/dc[ok]

    inp=pd.DataFrame({
      "u":act.u.to_numpy(float),
      "x":act.r_over_rs.to_numpy(float),
      "phi":act.phi.to_numpy(float),
      "f":act.A_f.to_numpy(float),
      "h":act.B_h.to_numpy(float),
      "X":act.X.to_numpy(float),
      "a1":q.a1_action.to_numpy(float),
      "c2":c2_target,
      "c4":q.c4_action.to_numpy(float),
      "F_tensor":q.F_tensor_action.to_numpy(float),
      "G_tensor":q.G_tensor_action.to_numpy(float),
      "H_tensor":q.H_tensor_action.to_numpy(float),
    })
    emitted=emit_from_primitives(inp,regularize_photon_root=True)
    emitted["G2XX_control"]=G2XX

    mask=(emitted.u.to_numpy(float)>=0.715)&(emitted.u.to_numpy(float)<50.0)
    d=emitted.loc[mask].sort_values("x").reset_index(drop=True)
    OUTDIR.mkdir(parents=True,exist_ok=True)
    d.to_csv(OUT,index=False)

    red=module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    per={str(L):health(d,red,L) for L in LS}

    # Primitive bridge regression at the inner/core edge only; deep deviations
    # are expected because this is the action-level quartic reconstruction.
    edge=(q.u.to_numpy(float)>=0.71)&(q.u.to_numpy(float)<=0.716)
    refmaps={
      "F_tensor":q.F_tensor_action.to_numpy(float),
      "G_tensor":q.G_tensor_action.to_numpy(float),
      "H_tensor":q.H_tensor_action.to_numpy(float),
      "a1":q.a1_action.to_numpy(float),
      "c4":q.c4_action.to_numpy(float),
    }
    edge_err={}
    for c,arr in refmaps.items():
        truth=ref[c].to_numpy(float)
        edge_err[c]=float(np.max(np.abs(arr[edge]-truth[edge])/np.maximum(1.,np.abs(truth[edge]))))

    report={
      "status":"CORE_QUARTIC_G5ZERO_ACTION_REGEN_COMPLETE",
      "domain":"0.715 <= u < 50 (deep full-G5 subcore excluded)",
      "rows":len(d),"u_range":[float(d.u.min()),float(d.u.max())],
      "source_hashes":{"action_g4xx":sha(ACT),"f2_c2_target":sha(REF)},
      "edge_primitive_scaled_errors":edge_err,
      "G2XX_control":{"min":float(np.min(G2XX[mask])),"max":float(np.max(G2XX[mask])),
                      "max_abs":float(np.max(np.abs(G2XX[mask]))),
                      "reachable_fraction":float(np.mean(~unresolved[mask])),
                      "unresolved_rows":int(np.sum(unresolved[mask])),
                      "note":"c2 is supplied as the documented selected primitive for this diagnostic. Rows where dc2/dG2XX vanishes are not claimed action-normal-jet reconstructed."},
      "exact_invariant_min":{
        "F":float(np.min(d.F_tensor)),"G":float(np.min(d.G_tensor)),
        "H":float(np.min(d.H_tensor)),"K_scalar":float(np.min(d.K_scalar))
      },
      "finite_L_health":per,
      "all_finite_L_health_pass":bool(all(v["min_K"]>0 and v["min_cr2"]>0 for v in per.values())),
      "guard":"This is a quartic G5=0 primitive-emission diagnostic. F,G,H,a1,c4 are direct action-derived; c2 is the selected primitive and is action-normal-jet certified only where dc2/dG2XX is nonzero. u>=50 requires the separate full-G5 bridge."
    }
    REPORT.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps(report,indent=2,allow_nan=False))
    return 0

if __name__=="__main__":raise SystemExit(main())
