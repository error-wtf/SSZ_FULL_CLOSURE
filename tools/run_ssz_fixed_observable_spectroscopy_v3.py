#!/usr/bin/env python3
"""Degeneracy-aware fixed-observable spectroscopy v3.

This supersedes the v2 basis-control failure at very high L.  Individual
eigenvector residues are not basis-invariant inside a numerically degenerate
eigenspace because the eigensolver may rotate that subspace.  The physical
invariants are:
  (a) generalized eigenvalues,
  (b) summed residue of a degenerate cluster,
  (c) the fixed-observable resolvent q^T [G-z^2 K]^-1 q.

v3 therefore:
  * solves G v = omega^2 K v with K-normalized modes;
  * tracks nondegenerate modes continuously in radius;
  * clusters near-degenerate modes by relative omega^2 gap;
  * computes cluster residue weights for fixed original-field observables;
  * tests exact resolvent invariance under a constant invertible field change;
  * tests common K,G rescaling;
  * reports only robust rank inversions with a finite margin.

This remains local-principal spectroscopy in the registered healthy window.
It is not a center-to-infinity QNM/retarded-Green spectrum.
"""
from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ssz_p5.config import DEFAULT_L  # noqa:E402
from ssz_p5.numerics import module  # noqa:E402
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central  # noqa:E402

OUT = ROOT / "data/generated/spectral/SSZ_SPECTROSCOPY_FIXED_OBSERVABLE_V3.json"
FULL_AUDIT = ROOT / "data/generated/spectral/CURRENT_MEMBER_FULL_AVAILABLE_DOMAIN_AUDIT.json"
DIRECT_CERTS = (
    ROOT / "data/certificates/SSZ_P5_DIRECT_GLOBAL_KRGM_CERTIFICATE.json",
    ROOT / "data/certificates/DIRECT_GLOBAL_KRGM_CERTIFICATE.json",
)
GLOBAL_EXPORT = ROOT / "data/generated/spectral/GLOBAL_CANONICAL_KRGM_SPECTRAL_EXPORT.npz"

REL_CLUSTER_GAP = 2e-6
ROBUST_MARGIN = 0.05


def sym(a):
    return (a + np.swapaxes(a, -1, -2)) / 2.0


def solve_modes(K, G):
    K = sym(K)
    G = sym(G)
    kw = np.linalg.eigvalsh(K)
    if np.min(kw) <= 0:
        raise ValueError("K_NOT_POSITIVE")
    lam, V = eigh(G, K, check_finite=True)
    if np.min(lam) <= 0:
        raise ValueError("RADIAL_PRINCIPAL_NOT_POSITIVE")
    return lam, np.sqrt(lam), V


def metric_overlap(Va, Vb, Ka, Kb):
    Kmid = (Ka + Kb) / 2.0
    M = Va.T @ Kmid @ Vb
    aa = np.sqrt(np.maximum(np.diag(Va.T @ Kmid @ Va), 1e-300))
    bb = np.sqrt(np.maximum(np.diag(Vb.T @ Kmid @ Vb), 1e-300))
    return M / (aa[:, None] * bb[None, :])


def best_perm(prevV, curV, prevK, curK):
    O = np.abs(metric_overlap(prevV, curV, prevK, curK))
    n = O.shape[0]
    best = None
    for p in itertools.permutations(range(n)):
        score = float(sum(O[i, p[i]] for i in range(n)))
        if best is None or score > best[0]:
            best = (score, p)
    return list(best[1])


def track(lams, omegas, vecs, Ks):
    tl = [lams[0].copy()]
    to = [omegas[0].copy()]
    tv = [vecs[0].copy()]
    minov = [1.0]
    for i in range(1, len(omegas)):
        p = best_perm(tv[-1], vecs[i], Ks[i-1], Ks[i])
        v = vecs[i][:, p].copy()
        l = lams[i][p].copy()
        o = omegas[i][p].copy()
        O = metric_overlap(tv[-1], v, Ks[i-1], Ks[i])
        for j in range(v.shape[1]):
            if O[j, j] < 0:
                v[:, j] *= -1.0
        d = np.abs(np.diag(metric_overlap(tv[-1], v, Ks[i-1], Ks[i])))
        minov.append(float(np.min(d)))
        tl.append(l); to.append(o); tv.append(v)
    return np.asarray(tl), np.asarray(to), np.asarray(tv), np.asarray(minov)


def clusters(lam):
    """Contiguous clusters in ascending eigenvalue order."""
    lam = np.asarray(lam, float)
    groups = []
    start = 0
    for i in range(len(lam)-1):
        gap = abs(lam[i+1]-lam[i]) / max(abs(lam[i+1]), abs(lam[i]), 1.0)
        if gap > REL_CLUSTER_GAP:
            groups.append(tuple(range(start, i+1)))
            start = i+1
    groups.append(tuple(range(start, len(lam))))
    return groups


def fixed_observables(n):
    out = {f"field_{i}": np.eye(n)[i] for i in range(n)}
    out["equal_weight"] = np.ones(n)/np.sqrt(n)
    q = np.arange(1,n+1,dtype=float)
    out["ramp_weight"] = q/np.linalg.norm(q)
    return out


def raw_residues(omega, V, q):
    return np.abs(q @ V)**2 / np.maximum(2.0*omega, 1e-300)


def cluster_weights(lam, omega, V, q):
    rr = raw_residues(omega, V, q)
    gs = clusters(lam)
    z = np.array([np.sum(rr[list(g)]) for g in gs], float)
    z /= max(float(np.sum(z)), 1e-300)
    centers = np.array([float(np.mean(lam[list(g)])) for g in gs])
    return gs, centers, z


def resolvent(K,G,q,z):
    A = G.astype(complex) - z*z*K.astype(complex)
    return complex(q @ np.linalg.solve(A,q))


def robust_pair_inversions(W, u, margin=ROBUST_MARGIN):
    """W N x M, fixed mode/cluster labels; require sign reversal by finite margin."""
    out=[]
    for a in range(W.shape[1]):
        for b in range(a+1,W.shape[1]):
            d=W[:,a]-W[:,b]
            if np.max(d)>margin and np.min(d)<-margin:
                ia=int(np.argmax(d)); ib=int(np.argmin(d))
                out.append({
                    "a":a,"b":b,
                    "max_diff":float(np.max(d)),"u_max":float(u[ia]),
                    "min_diff":float(np.min(d)),"u_min":float(u[ib]),
                })
    return out


def basis_resolvent_control(Ks,Gs,obs,u_indices):
    n=Ks[0].shape[0]
    A=np.eye(n)
    for i in range(n):
        for j in range(i+1,n):
            A[i,j]=0.17/(1+j-i)
    max_eval=0.0
    max_res_rel=0.0
    max_cluster=0.0
    checked_clusters=0

    for k in u_indices:
        lam0,w0,V0=solve_modes(Ks[k],Gs[k])
        K1=A.T@Ks[k]@A
        G1=A.T@Gs[k]@A
        lam1,w1,V1=solve_modes(K1,G1)
        max_eval=max(max_eval,float(np.max(np.abs(lam0-lam1))))

        # Scalar resolvent is exact basis invariant.
        probes=np.linspace(float(w0.min()*0.8),float(w0.max()*1.2),7)
        eta=max(1e-3*float(np.median(w0)),1e-12)
        for q in obs.values():
            q1=A.T@q
            for wr in probes:
                z=complex(wr,eta)
                r0=resolvent(Ks[k],Gs[k],q,z)
                r1=resolvent(K1,G1,q1,z)
                rel=abs(r0-r1)/max(abs(r0),abs(r1),1.0)
                max_res_rel=max(max_res_rel,float(rel))

            # Compare only cluster-summed residues, matching clusters by eigenvalues.
            g0,c0,z0=cluster_weights(lam0,w0,V0,q)
            g1,c1,z1=cluster_weights(lam1,w1,V1,q1)
            if len(g0)==len(g1) and np.max(np.abs(c0-c1))<1e-7:
                max_cluster=max(max_cluster,float(np.max(np.abs(z0-z1))))
                checked_clusters+=1

    return {
        "max_abs_generalized_eigenvalue_delta":max_eval,
        "max_relative_resolvent_delta":max_res_rel,
        "max_abs_cluster_residue_delta":max_cluster,
        "cluster_comparisons":checked_clusters,
        "pass":bool(max_eval<1e-7 and max_res_rel<1e-8 and max_cluster<2e-6),
    }


def scale_control(Ks,Gs,obs,u_indices):
    max_eval=0.0
    max_w=0.0
    for kk,k in enumerate(u_indices):
        s=0.83+0.31*(kk/max(len(u_indices)-1,1))
        l0,w0,V0=solve_modes(Ks[k],Gs[k])
        l1,w1,V1=solve_modes(s*Ks[k],s*Gs[k])
        max_eval=max(max_eval,float(np.max(np.abs(l0-l1))))
        for q in obs.values():
            # K-normalized eigenvectors scale, but normalized residue weights should not.
            z0=raw_residues(w0,V0,q); z0/=max(float(np.sum(z0)),1e-300)
            z1=raw_residues(w1,V1,q); z1/=max(float(np.sum(z1)),1e-300)
            max_w=max(max_w,float(np.max(np.abs(z0-z1))))
    return {
        "max_abs_generalized_eigenvalue_delta":max_eval,
        "max_abs_normalized_residue_delta":max_w,
        "pass":bool(max_eval<1e-7 and max_w<2e-6),
    }


def global_gate():
    blockers=[]
    if FULL_AUDIT.exists():
        try:
            d=json.loads(FULL_AUDIT.read_text())
            if not d.get("full_available_domain_K_pass",False):
                blockers.append("FULL_NATIVE_DOMAIN_K_FAIL")
            if not d.get("full_available_domain_radial_pass",False):
                blockers.append("FULL_NATIVE_DOMAIN_RADIAL_FAIL")
        except Exception:
            blockers.append("INVALID_FULL_NATIVE_DOMAIN_AUDIT")
    else:
        blockers.append("MISSING_FULL_NATIVE_DOMAIN_AUDIT")
    if not any(p.exists() for p in DIRECT_CERTS):
        blockers.append("MISSING_DIRECT_GLOBAL_KRGM_CERTIFICATE")
    if not GLOBAL_EXPORT.exists():
        blockers.append("MISSING_GLOBAL_CANONICAL_KRGM_SPECTRAL_EXPORT")
    return {"status":"READY" if not blockers else "BLOCKED","blockers":blockers}


def per_L(stream,reducer,L):
    a=reducer.canonical_audit(stream,int(L))
    Kall=sym(np.asarray(a["K"],float))
    Gall=sym(np.asarray(a["G"],float))
    uall=stream.u.to_numpy(float)
    ids=np.flatnonzero((uall>0.62)&(uall<0.70))
    ids=ids[np.unique(np.linspace(0,len(ids)-1,min(151,len(ids))).round().astype(int))]
    u=uall[ids]; Ks=Kall[ids]; Gs=Gall[ids]

    Ls=[]; Ws=[]; Vs=[]
    for K,G in zip(Ks,Gs):
        lam,w,V=solve_modes(K,G)
        Ls.append(lam); Ws.append(w); Vs.append(V)
    Ls=np.asarray(Ls); Ws=np.asarray(Ws); Vs=np.asarray(Vs)
    tl,tw,tv,ov=track(Ls,Ws,Vs,Ks)
    n=tw.shape[1]
    obs=fixed_observables(n)

    # Only compare individual mode residues where no local degeneracy occurs.
    nongenerate=np.ones(len(u),dtype=bool)
    min_gap=np.full(len(u),np.inf)
    for k in range(len(u)):
        gaps=np.diff(np.sort(tl[k]))
        rel=gaps/np.maximum(np.maximum(np.abs(tl[k][:-1]),np.abs(tl[k][1:])),1.0)
        min_gap[k]=float(np.min(rel)) if len(rel) else np.inf
        nongenerate[k]=bool(min_gap[k]>REL_CLUSTER_GAP)

    ores={}
    any_robust=False
    for name,q in obs.items():
        W=np.asarray([
            raw_residues(tw[k],tv[k],q)/max(float(np.sum(raw_residues(tw[k],tv[k],q))),1e-300)
            for k in range(len(u))
        ])
        inv=robust_pair_inversions(W[nongenerate],u[nongenerate]) if np.sum(nongenerate)>=3 else []
        any_robust|=bool(inv)

        # Basis-independent cluster concentration at each radius.
        top1=[]; top2=[]; neff=[]; nclusters=[]
        for k in range(len(u)):
            gs,cc,z=cluster_weights(tl[k],tw[k],tv[k],q)
            zz=np.sort(z)[::-1]
            top1.append(float(zz[0]))
            top2.append(float(np.sum(zz[:min(2,len(zz))])))
            neff.append(float(1.0/np.sum(z*z)))
            nclusters.append(len(gs))
        ores[name]={
            "robust_nondegenerate_pair_inversions":inv[:30],
            "robust_nondegenerate_inversion_count":len(inv),
            "cluster_top1_range":[float(min(top1)),float(max(top1))],
            "cluster_top2_range":[float(min(top2)),float(max(top2))],
            "cluster_neff_range":[float(min(neff)),float(max(neff))],
            "cluster_count_range":[int(min(nclusters)),int(max(nclusters))],
        }

    cidx=list(range(0,len(u),max(1,len(u)//14)))
    bc=basis_resolvent_control(Ks,Gs,obs,cidx)
    sc=scale_control(Ks,Gs,obs,cidx)

    # Coarse robustness of robust inversions.
    coarse_ok=True
    for name,q in obs.items():
        W=np.asarray([
            raw_residues(tw[k],tv[k],q)/max(float(np.sum(raw_residues(tw[k],tv[k],q))),1e-300)
            for k in range(len(u))
        ])
        fine=bool(robust_pair_inversions(W[nongenerate],u[nongenerate])) if np.sum(nongenerate)>=3 else False
        csel=np.arange(0,len(u),2)
        csel=csel[nongenerate[csel]]
        coarse=bool(robust_pair_inversions(W[csel],u[csel])) if len(csel)>=3 else False
        if fine!=coarse:
            coarse_ok=False

    ratio_spans={}
    for i in range(n):
        for j in range(i+1,n):
            rr=tw[:,i]/tw[:,j]
            ratio_spans[f"{i}/{j}"]=float((np.max(rr)-np.min(rr))/max(float(np.mean(rr)),1e-300))

    return {
        "L":int(L),
        "u_min":float(u.min()),"u_max":float(u.max()),"samples":len(u),
        "min_tracking_overlap":float(np.min(ov)),
        "min_relative_eigenvalue_gap":float(np.min(min_gap)),
        "nondegenerate_sample_fraction":float(np.mean(nongenerate)),
        "any_robust_fixed_observable_inversion":bool(any_robust),
        "non_common_frequency_scaling_detected_1e3":bool(any(v>1e-3 for v in ratio_spans.values())),
        "frequency_ratio_relative_spans":ratio_spans,
        "observables":ores,
        "controls":{
            "basis_resolvent_and_cluster":bc,
            "common_KG_scale":sc,
            "coarse_grid_same_robust_inversion_boolean":coarse_ok,
            "pass":bool(bc["pass"] and sc["pass"] and coarse_ok and np.min(ov)>0.9),
        },
    }


def main():
    b=build_onshell_central(ROOT)
    stream=b.direct41.sort_values("x").reset_index(drop=True)
    reducer=module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    per={str(int(L)):per_L(stream,reducer,int(L)) for L in DEFAULT_L}
    all_controls=all(v["controls"]["pass"] for v in per.values())
    any_inv=any(v["any_robust_fixed_observable_inversion"] for v in per.values())
    any_noncommon=any(v["non_common_frequency_scaling_detected_1e3"] for v in per.values())
    if not all_controls:
        status="LOCAL_SPECTROSCOPY_V3_CONTROL_FAIL"
    elif any_inv:
        status="LOCAL_FIXED_OBSERVABLE_ROBUST_RESIDUE_REORDERING_DETECTED"
    else:
        status="LOCAL_FIXED_OBSERVABLE_ROBUST_RESIDUE_REORDERING_NULL"

    report={
        "global_physical_spectroscopy":global_gate(),
        "local_principal_spectroscopy_v3":{
            "status":status,
            "scope":"registered healthy production window only",
            "method":"generalized eigenproblem + fixed observables + degeneracy clusters + invariant resolvent controls",
            "all_controls_pass":all_controls,
            "any_robust_residue_reordering":any_inv,
            "any_non_common_frequency_scaling":any_noncommon,
            "per_L":per,
        },
        "interpretation_guard":"Global QNM/retarded-Green residue claims remain blocked until healthy direct-global KRGSM and certificate exist.",
    }
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps({
        "global":report["global_physical_spectroscopy"],
        "status":status,
        "all_controls":all_controls,
        "any_robust_residue_reordering":any_inv,
        "any_non_common_frequency_scaling":any_noncommon,
        "per_L":{
            L:{
                "control":v["controls"]["pass"],
                "basis":v["controls"]["basis_resolvent_and_cluster"],
                "robust_inversion":v["any_robust_fixed_observable_inversion"],
                "nondeg_frac":v["nondegenerate_sample_fraction"],
                "min_overlap":v["min_tracking_overlap"],
            } for L,v in per.items()
        }
    },indent=2,allow_nan=False))
    return 0 if all_controls else 1

if __name__=="__main__":
    raise SystemExit(main())
