#!/usr/bin/env python3
"""Convergence and k-scan audit for SSZ production-window spectroscopy.

This is a second-stage falsification audit for the diagnostic spectroscopy suite.

It asks two questions:

1. Box convergence / boundary robustness:
   Are negative omega^2 modes and the low positive spectrum stable when the
   artificial Dirichlet box resolution and subwindow are changed?

2. Local frozen-coefficient k scan:
   At fixed healthy production radii, does the local Hermitian pencil
       H(k)=k^2 G + i k S - M
   enter an all-positive omega^2 band at finite radial k, and do tracked
   frequency ratios / probe weights change there?

All K,G,S,M arrays are reduced on the full production-resolution grid before
any sampling.  No result is a center-to-infinity QNM spectrum.
"""

from __future__ import annotations

import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ssz_p5.numerics import module  # noqa:E402
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central  # noqa:E402

OUT = ROOT / "data/generated/spectral/PRODUCTION_SPECTROSCOPY_CONVERGENCE_AUDIT.json"
LS = (6, 42, 110, 420, 1000)
NODE_COUNTS = (81, 121, 161)
WINDOWS = ((0.628, 0.692), (0.632, 0.688))
RADII = (0.632, 0.648, 0.664, 0.680, 0.692)
K_FACTORS = (0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0)


def load_spec_module():
    p = ROOT / "tools/run_production_spectroscopy_suite.py"
    s = importlib.util.spec_from_file_location("prod_spec", p)
    m = importlib.util.module_from_spec(s)
    assert s.loader is not None
    s.loader.exec_module(m)
    return m


def sample_ids(u, lo, hi, n):
    ids = np.flatnonzero((u > lo) & (u < hi))
    if len(ids) < n:
        raise ValueError(f"window {lo},{hi} has only {len(ids)} rows for n={n}")
    pick = np.linspace(0, len(ids)-1, n).round().astype(int)
    return ids[np.unique(pick)]


def positive_omega_first(rec, n=8):
    return rec["omega_first"][:n]


def relative_vector_difference(a, b):
    n = min(len(a), len(b))
    if n == 0:
        return None
    aa = np.asarray(a[:n], float)
    bb = np.asarray(b[:n], float)
    return float(np.max(np.abs(aa-bb) / np.maximum(1.0, np.abs(bb))))


def box_audit(prod, A, specmod, L):
    uall = prod.u.to_numpy(float)
    xall = prod.x.to_numpy(float)
    fall = prod.f.to_numpy(float)
    Kf, Gf, Sf, Mf = (np.asarray(A[k], float) for k in ("K","G","S","M"))

    grid = {}
    for lo, hi in WINDOWS:
        wkey = f"{lo:.3f}_{hi:.3f}"
        grid[wkey] = {}
        for n in NODE_COUNTS:
            ids = sample_ids(uall, lo, hi, n)
            K,G,S,M = Kf[ids],Gf[ids],Sf[ids],Mf[ids]
            x,u,f = xall[ids],uall[ids],fall[ids]
            phys = specmod.boxed_spectrum_from_arrays(K,G,S,M,x,u,f)
            mid = len(ids)//2
            ctrl = specmod.boxed_spectrum_from_arrays(
                np.repeat(K[mid][None,:,:], len(ids), 0),
                np.repeat(G[mid][None,:,:], len(ids), 0),
                np.repeat(S[mid][None,:,:], len(ids), 0),
                np.repeat(M[mid][None,:,:], len(ids), 0),
                x,u,np.repeat(f[mid],len(ids)),
            )
            grid[wkey][str(n)] = {
                "nodes": int(len(ids)),
                "physical": {
                    "omega2_min": float(phys["zero_or_negative_omega2_min"]),
                    "negative_count": int(phys["negative_omega2_count"]),
                    "negative_fraction": float(phys["negative_omega2_count"]/phys["dimension"]),
                    "omega_first_positive": positive_omega_first(phys),
                },
                "control": {
                    "omega2_min": float(ctrl["zero_or_negative_omega2_min"]),
                    "negative_count": int(ctrl["negative_omega2_count"]),
                    "negative_fraction": float(ctrl["negative_omega2_count"]/ctrl["dimension"]),
                    "omega_first_positive": positive_omega_first(ctrl),
                },
            }

    # convergence relative to highest node count within each window
    convergence = {}
    for wkey, rows in grid.items():
        ref = rows[str(max(NODE_COUNTS))]
        convergence[wkey] = {}
        for n in NODE_COUNTS[:-1]:
            cur = rows[str(n)]
            convergence[wkey][str(n)] = {
                "physical_first_positive_max_rel_diff_vs_161": relative_vector_difference(
                    cur["physical"]["omega_first_positive"],
                    ref["physical"]["omega_first_positive"],
                ),
                "physical_omega2_min_abs_diff_vs_161": float(
                    abs(cur["physical"]["omega2_min"]-ref["physical"]["omega2_min"])
                ),
                "control_first_positive_max_rel_diff_vs_161": relative_vector_difference(
                    cur["control"]["omega_first_positive"],
                    ref["control"]["omega_first_positive"],
                ),
            }

    # boundary robustness: compare two windows at the common 161-node resolution
    a = grid[f"{WINDOWS[0][0]:.3f}_{WINDOWS[0][1]:.3f}"]["161"]["physical"]
    b = grid[f"{WINDOWS[1][0]:.3f}_{WINDOWS[1][1]:.3f}"]["161"]["physical"]
    boundary = {
        "omega2_min_sign_same": bool(np.sign(a["omega2_min"]) == np.sign(b["omega2_min"])),
        "negative_fraction_abs_diff": float(abs(a["negative_fraction"]-b["negative_fraction"])),
        "first_positive_max_rel_diff": relative_vector_difference(
            a["omega_first_positive"], b["omega_first_positive"]
        ),
    }
    return {"grid":grid,"convergence":convergence,"boundary_robustness":boundary}


def local_k_scan(prod, A, specmod, L):
    u = prod.u.to_numpy(float)
    x = prod.x.to_numpy(float)
    K,G,S,M = (np.asarray(A[k]) for k in ("K","G","S","M"))
    rows = specmod.nearest_rows(u, RADII)
    dx = abs(float(x[rows[-1]]-x[rows[0]]))
    k0 = math.pi/max(dx,1e-12)

    records = []
    for row in rows:
        radius = {
            "u": float(u[row]),
            "x": float(x[row]),
            "scan": [],
            "first_all_positive_factor": None,
        }
        prev = None
        for fac in K_FACTORS:
            kval = fac*k0
            o2, Uw, V, hres = specmod.local_spectrum(
                K[row],G[row],S[row],M[row],kval
            )
            if prev is not None:
                perm = specmod.best_overlap_permutation(prev,Uw)
                o2=o2[perm]; Uw=Uw[:,perm]; V=V[:,perm]
            prev=Uw
            weights=specmod.normalized_probe_weights(V,o2)
            allpos=bool(np.all(o2>0))
            if allpos and radius["first_all_positive_factor"] is None:
                radius["first_all_positive_factor"]=float(fac)
            radius["scan"].append({
                "k_factor":float(fac),
                "k":float(kval),
                "omega2":[float(z) for z in o2],
                "omega":[float(np.sqrt(z)) if z>0 else None for z in o2],
                "all_positive":allpos,
                "probe":weights,
                "hermiticity_residual":float(hres),
            })
        records.append(radius)

    # Only use all-positive points to summarize observable weight switching.
    stable_switches={p:0 for p in specmod.PROBES}
    stable_points=0
    for rad in records:
        prev_dom={p:None for p in specmod.PROBES}
        for rec in rad["scan"]:
            if not rec["all_positive"]:
                continue
            stable_points += 1
            for p in specmod.PROBES:
                dom=rec["probe"][p]["dominant_mode"]
                if prev_dom[p] is not None and dom!=prev_dom[p]:
                    stable_switches[p]+=1
                prev_dom[p]=dom

    return {
        "k_unit":float(k0),
        "records":records,
        "stable_all_positive_scan_points":int(stable_points),
        "probe_dominance_switches_within_all_positive_points":stable_switches,
    }


def main():
    specmod=load_spec_module()
    build=build_onshell_central(ROOT)
    d=build.direct41.sort_values("x").reset_index(drop=True)
    prod=d[(d.u>0.62)&(d.u<0.70)].reset_index(drop=True)
    red=module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")

    box={}
    kscan={}
    for L in LS:
        A=red.canonical_audit(prod,int(L))
        box[str(L)]=box_audit(prod,A,specmod,L)
        kscan[str(L)]=local_k_scan(prod,A,specmod,L)

    # classification is descriptive, not a claim of physical QNM instability.
    robust_negative = {}
    for L,rec in box.items():
        signs=[]
        for w in rec["grid"].values():
            for n in w.values():
                signs.append(n["physical"]["omega2_min"]<0)
        robust_negative[L]=bool(all(signs))

    any_stable_weight_switch=any(
        any(v>0 for v in rec["probe_dominance_switches_within_all_positive_points"].values())
        for rec in kscan.values()
    )

    report={
        "scope":"diagnostic spectroscopy convergence; full-resolution reduction first; not QNM",
        "box":box,
        "local_k_scan":kscan,
        "summary":{
            "negative_omega2_persists_all_box_windows_and_resolutions":robust_negative,
            "any_probe_dominance_switch_in_all_positive_local_k_points":any_stable_weight_switch,
        },
        "interpretation_guard":(
            "Persistent boxed negative omega^2 is a lower-order/box-spectrum diagnostic. "
            "It is not a physical instability claim until boundary-condition and global "
            "same-action KRGM/QNM convergence gates are satisfied."
        ),
    }
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps(report["summary"],indent=2))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
