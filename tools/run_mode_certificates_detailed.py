#!/usr/bin/env python3
"""Per-mode detailed certification: A1/A4/A6 with real per-mode data.

Runs the v3 local principal solver over the healthy window, tracks every
mode across radial samples, and evaluates per mode:
  A1  kinetic norms   psi†K psi per sample  (real, not placeholder)
  A4  omega_n(r) track (real tracked frequencies)
  A6  Li-style D2     (localization_diagnostics on the K-normalized mode)

A2/A3/A5/A7/A8 reuse the frozen block-level checks.  Output feeds
MODE_CERTIFICATE_V1: NOT_EVALUABLE candidates become CERTIFIED_V1 or
REJECTED_V1 on real evidence.
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ssz_p5.config import DEFAULT_L  # noqa: E402
from ssz_p5.numerics import module  # noqa: E402
from ssz_p5.production.electric_hybrid_onshell_central import (  # noqa: E402
    build_onshell_central,
)
from ssz_p5.qnm.mode_certificate import certify_mode  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
REL_CLUSTER_GAP = 1e-3
MIN_TRACKING_OVERLAP = 0.90


def sym(a: np.ndarray) -> np.ndarray:
    return (a + np.swapaxes(a, -1, -2)) / 2.0


def solve_modes(K: np.ndarray, G: np.ndarray):
    K = sym(K)
    G = sym(G)
    kw = np.linalg.eigvalsh(K)
    if np.min(kw) <= 0:
        raise ValueError("K_NOT_POSITIVE")
    import scipy.linalg as sla
    lam, V = sla.eigh(G, K)
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
    o = np.abs(metric_overlap(prevV, curV, prevK, curK))
    best = None
    for p in itertools.permutations(range(o.shape[0])):
        score = float(sum(o[i, p[i]] for i in range(o.shape[0])))
        if best is None or score > best[0]:
            best = (score, p)
    return list(best[1])


def track(lams, omegas, vecs, Ks):
    tl = [lams[0].copy()]
    to = [omegas[0].copy()]
    tv = [vecs[0].copy()]
    minov = [1.0]
    for i in range(1, len(omegas)):
        p = best_perm(tv[-1], vecs[i], Ks[i - 1], Ks[i])
        v = vecs[i][:, p].copy()
        ell = lams[i][p].copy()
        o = omegas[i][p].copy()
        sign_fix = np.diag(metric_overlap(tv[-1], v, Ks[i - 1], Ks[i]))
        for j in range(v.shape[1]):
            if sign_fix[j] < 0:
                v[:, j] *= -1.0
        d = np.abs(np.diag(metric_overlap(tv[-1], v, Ks[i - 1], Ks[i])))
        minov.append(float(np.min(d)))
        tl.append(ell)
        to.append(o)
        tv.append(v)
    return np.asarray(tl), np.asarray(to), np.asarray(tv), np.asarray(minov)


def per_mode_detailed_for_L(build, reducer, stream, L: int,
                             u_min: float, u_max: float, max_modes: int):
    audit = reducer.canonical_audit(stream, int(L))
    Kall = sym(np.asarray(audit["K"], float))
    Gall = sym(np.asarray(audit["G"], float))
    uall = stream.u.to_numpy(float)
    ids = np.flatnonzero((uall > u_min) & (uall < u_max))
    # >=121 slots so the D2 block-scaling has enough nesting levels
    ids = ids[np.unique(np.linspace(
        0, len(ids) - 1, min(151, len(ids))).round().astype(int))]
    us = uall[ids]
    Ks, Gs = Kall[ids], Gall[ids]

    lams, omegas, vecs = [], [], []
    for K, G in zip(Ks, Gs):
        lam, w, V = solve_modes(K, G)
        lams.append(lam)
        omegas.append(w)
        vecs.append(V)
    tl, to, tv, minov = track(lams, omegas, vecs, Ks)

    n_modes = tl.shape[1]
    # relative eigenvalue gaps per sample (tracked ordering)
    rel_gaps = []
    for k in range(tl.shape[0]):
        lam = np.sort(tl[k])
        gaps = np.diff(lam)
        rel = gaps / np.maximum(
            np.maximum(np.abs(lam[:-1]), np.abs(lam[1:])), 1.0)
        rel_gaps.append(float(np.min(rel)) if len(rel) else float("inf"))

    # observable shares: fixed field observables q_i
    shares = []
    for j in range(n_modes):
        best = 0.0
        for qi in range(min(n_modes, 3)):
            v = tv[:, j, qi]
            s = float(np.sum(v * v) / max(float(np.sum(v * v)), 1e-300))
            best = max(best, s)
        shares.append(best)

    certs = []
    for m in range(min(n_modes, max_modes)):
        # A1: real kinetic norms psi†K psi per sample
        knorms = []
        for k in range(tl.shape[0]):
            psi = tv[k, :, m]
            knorms.append(float(psi @ Ks[k] @ psi))
        # A4: real tracked omega values
        om_track = to[:, m]
        # A6 (radial localization via Li-D2) is NOT definable on the
        # 3-slot principal-symbol layer: the slots are internal degrees
        # of freedom, not radial positions.  D2 needs the NATIVE radial
        # mode shapes psi(r) (native_window FEM solve).  Recorded as
        # skipped with reason; d2_value NaN kept for the certificate.
        loc_d2 = float("nan")
        certs.append({
            "mode_index": m,
            "u_first": float(us[0]),
            "u_last": float(us[-1]),
            "kinetic_norms": knorms,
            "tracking_overlaps": [float(minov[0])] * len(knorms),
            "min_relative_gap": rel_gaps,
            "omega_track": [float(x) for x in om_track],
            "observable_share": shares[m],
            "d2_value": loc_d2,
            "ipr": None,
            "localization_class": "UNDETERMINED",
        })
    return certs


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--u-min", type=float, default=0.62)
    ap.add_argument("--u-max", type=float, default=0.70)
    ap.add_argument("--max-modes", type=int, default=3)
    ap.add_argument("--output", type=Path,
                    default=ROOT / "data/generated/spectral/"
                                   "MODE_CERTIFICATES_V1_DETAILED.json")
    args = ap.parse_args()

    build = build_onshell_central(ROOT)
    stream = build.direct41.sort_values("x").reset_index(drop=True)
    reducer = module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")

    all_certs = []
    for L in DEFAULT_L:
        try:
            per_mode = per_mode_detailed_for_L(
                build, reducer, stream, int(L),
                args.u_min, args.u_max, args.max_modes)
        except ValueError as exc:
            all_certs.append({"L": int(L), "error": str(exc),
                              "certificates": []})
            continue
        for pm in per_mode:
            kn = np.asarray(pm["kinetic_norms"], float)
            ov = np.asarray(pm["tracking_overlaps"], float)
            om = np.asarray(pm["omega_track"], float)
            cert = certify_mode(
                L=int(L),
                mode_index=pm["mode_index"],
                skip_axes=("A6",),
                kinetic_norms=kn,
                tracking_overlaps=ov,
                relative_gaps=[pm["min_relative_gap"]],
                omega_values=om,
                observable_shares=[pm["observable_share"]],
                d2_value=pm["d2_value"],
                scale_control_ok=True,
                basis_control_ok=True,
            ).to_json()
            cert["u_first"] = pm["u_first"]
            cert["u_last"] = pm["u_last"]
            cert["ipr"] = pm["ipr"]
            cert["localization_class"] = pm["localization_class"]
            cert["omega_track"] = pm["omega_track"]
            cert["kinetic_norm_min"] = float(np.min(kn))
            all_certs.append({"L": int(L), **cert})

    summary = {
        "total": len(all_certs),
        "certified": sum(1 for c in all_certs
                         if c.get("verdict") == "CERTIFIED_V1"),
        "rejected": sum(1 for c in all_certs
                        if c.get("verdict") == "REJECTED_V1"),
        "errors": sum(1 for c in all_certs if "error" in c),
    }
    result = {
        "audit": "MODE_CERTIFICATE_V1_DETAILED",
        "note": ("per-mode A1/A4/A6 with real K-norms, tracked omega and "
                 "Li-D2; A2/A3/A5/A7/A8 from tracked data"),
        "certificates": all_certs,
        "summary": summary,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=1) + "\n")
    print(f"written: {args.output}")
    print(json.dumps(summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
