#!/usr/bin/env python3
"""SSZ spectroscopy suite: guarded global gate + honest local principal spectroscopy.

This suite finally executes the spectroscopy question without conflating three
different objects:

1. PHYSICAL GLOBAL SPECTROSCOPY
   Requires a hash-bound DIRECT_GLOBAL_KRGM certificate/export and a healthy
   full native member.  If either is absent, this gate is BLOCKED, not guessed.

2. LOCAL PRINCIPAL SPECTROSCOPY
   Executed only on the registered healthy production window 0.62<u<0.70.
   It diagonalizes the canonically normalized principal pencil
       C(u)=K(u)^(-1/2) G(u) K(u)^(-1/2)
   and computes observable-channel weights
       W_{a n}(u)=| e_a^T v_n(u) |^2
   with continuous mode tracking by overlap.  This is a local characteristic
   response diagnostic, not a QNM spectrum.

3. CONTROLS
   - common positive rescaling K,G -> s(u)K,s(u)G must leave C and weights invariant;
   - coarse-grid/subsample tracking must preserve the qualitative reordering result;
   - K positivity is enforced before any local spectral statement.

No physical QNM/residue claim is emitted unless the direct global gate is open.
"""

from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ssz_p5.config import DEFAULT_L  # noqa:E402
from ssz_p5.numerics import module  # noqa:E402
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central  # noqa:E402
from ssz_p5.provenance.manifest import sha256  # noqa:E402
from ssz_p5.qnm.gate import require_direct_krgm_certificate  # noqa:E402

OUT = ROOT / "data/generated/spectral/SSZ_SPECTROSCOPY_SUITE_REPORT.json"
DIRECT = ROOT / "data/certificates/DIRECT_GLOBAL_KRGM_CERTIFICATE.json"
EXPORT = ROOT / "data/generated/spectral/GLOBAL_CANONICAL_KRGM_SPECTRAL_EXPORT.npz"
FULL_AUDIT = ROOT / "data/generated/spectral/CURRENT_MEMBER_FULL_AVAILABLE_DOMAIN_AUDIT.json"


def sym(a):
    return (a + np.swapaxes(a, -1, -2)) / 2.0


def canon(K, G):
    w, U = np.linalg.eigh(sym(K))
    if np.min(w) <= 0:
        raise ValueError("K not positive")
    inv = U @ np.diag(1.0 / np.sqrt(w)) @ U.T
    C = sym(inv @ sym(G) @ inv)
    ev, V = np.linalg.eigh(C)
    return ev, V, C


def best_perm(prev, cur):
    """Small-F exact overlap matching."""
    F = prev.shape[1]
    ov = np.abs(prev.T @ cur)
    best = None
    best_score = -1.0
    for p in itertools.permutations(range(F)):
        score = float(sum(ov[i, p[i]] for i in range(F)))
        if score > best_score:
            best_score = score
            best = p
    return np.array(best, dtype=int), best_score / F


def track_modes(Vs):
    tracked = [Vs[0].copy()]
    overlaps = [1.0]
    perms = [list(range(Vs[0].shape[1]))]
    for cur in Vs[1:]:
        prev = tracked[-1]
        p, score = best_perm(prev, cur)
        cc = cur[:, p].copy()
        # continuous signs
        for j in range(cc.shape[1]):
            if float(prev[:, j] @ cc[:, j]) < 0:
                cc[:, j] *= -1
        tracked.append(cc)
        overlaps.append(float(score))
        perms.append(p.tolist())
    return np.asarray(tracked), np.asarray(overlaps), perms


def channel_weights(tracked):
    # tracked shape N,F,F; columns=modes, rows=canonical observable channels
    W = np.abs(tracked) ** 2
    # output channel x mode x radius
    return np.transpose(W, (1, 2, 0))


def any_rank_inversion(weights, tol=1e-8):
    # weights shape A,M,N.  For each observable channel, compare mode pairs.
    witnesses = []
    A, M, N = weights.shape
    for a in range(A):
        for i in range(M):
            for j in range(i + 1, M):
                d = weights[a, i] - weights[a, j]
                if np.any(d > tol) and np.any(d < -tol):
                    witnesses.append({"channel": a, "mode_a": i, "mode_b": j})
    return bool(witnesses), witnesses[:20]


def local_for_L(stream, reducer, L):
    a = reducer.canonical_audit(stream, int(L))
    K = np.asarray(a["K"], float)
    G = np.asarray(a["G"], float)
    u = stream.u.to_numpy(float)
    mask = (u > 0.62) & (u < 0.70)
    idx = np.flatnonzero(mask)
    # deterministic thinning for tractability, preserving endpoints.
    pick = np.unique(np.linspace(0, len(idx) - 1, min(81, len(idx))).round().astype(int))
    idx = idx[pick]
    uu = u[idx]

    eigvals, Vs, Cs = [], [], []
    kmins = []
    for i in idx:
        kmins.append(float(np.min(np.linalg.eigvalsh(sym(K[i])))))
        ev, V, C = canon(K[i], G[i])
        eigvals.append(ev)
        Vs.append(V)
        Cs.append(C)
    eigvals = np.asarray(eigvals)
    Vs = np.asarray(Vs)
    Cs = np.asarray(Cs)
    tracked, overlaps, perms = track_modes(Vs)
    weights = channel_weights(tracked)
    reordered, witnesses = any_rank_inversion(weights)

    # Common-rescaling null control.
    s = 0.7 + 0.6 * (uu - uu.min()) / max(uu.max() - uu.min(), 1e-15)
    max_c_delta = 0.0
    max_w_delta = 0.0
    V2 = []
    for q, i in enumerate(idx):
        _, vv, cc = canon(s[q] * K[i], s[q] * G[i])
        max_c_delta = max(max_c_delta, float(np.max(np.abs(cc - Cs[q]))))
        V2.append(vv)
    V2t, _, _ = track_modes(np.asarray(V2))
    W2 = channel_weights(V2t)
    max_w_delta = float(np.max(np.abs(W2 - weights)))

    # Coarse-grid qualitative control.
    coarse = np.arange(0, len(idx), 2)
    coarse_reordered, _ = any_rank_inversion(weights[:, :, coarse])

    # Principal characteristic ratios after removing common scale.
    positive = np.maximum(eigvals, 1e-300)
    ratios = positive / positive[:, :1]
    ratio_span = np.max(ratios, axis=0) - np.min(ratios, axis=0)

    return {
        "L": int(L),
        "nodes": int(len(idx)),
        "u_min": float(uu.min()),
        "u_max": float(uu.max()),
        "min_K": float(min(kmins)),
        "mode_count": int(eigvals.shape[1]),
        "min_tracking_overlap": float(np.min(overlaps)),
        "rank_reordering_detected": reordered,
        "rank_reordering_witnesses": witnesses,
        "coarse_grid_same_reordering_boolean": bool(coarse_reordered == reordered),
        "common_rescaling_control": {
            "max_abs_C_delta": max_c_delta,
            "max_abs_weight_delta": max_w_delta,
            "pass": bool(max_c_delta < 1e-9 and max_w_delta < 1e-8),
        },
        "normalized_characteristic_ratio_span": [float(x) for x in ratio_span],
        "weights_minmax_by_channel_mode": [
            [
                [float(np.min(weights[a, m])), float(np.max(weights[a, m]))]
                for m in range(weights.shape[1])
            ]
            for a in range(weights.shape[0])
        ],
    }


def global_gate():
    blockers = []
    direct_hash = None
    export_hash = None

    if FULL_AUDIT.is_file():
        try:
            d = json.loads(FULL_AUDIT.read_text())
            if not d.get("full_available_domain_K_pass", False):
                blockers.append("FULL_NATIVE_DOMAIN_K_FAIL")
            if not d.get("full_available_domain_radial_pass", False):
                blockers.append("FULL_NATIVE_DOMAIN_RADIAL_FAIL")
        except Exception as exc:
            blockers.append(f"INVALID_FULL_DOMAIN_AUDIT:{type(exc).__name__}")
    else:
        blockers.append("MISSING_FULL_DOMAIN_HEALTH_AUDIT")

    if not DIRECT.is_file():
        blockers.append("MISSING_DIRECT_GLOBAL_KRGM_CERTIFICATE")
    else:
        try:
            require_direct_krgm_certificate(DIRECT, ROOT)
            direct_hash = sha256(DIRECT)
        except Exception as exc:
            blockers.append(f"INVALID_DIRECT_GLOBAL_KRGM_CERTIFICATE:{type(exc).__name__}")

    if not EXPORT.is_file():
        blockers.append("MISSING_GLOBAL_CANONICAL_KRGM_SPECTRAL_EXPORT")
    else:
        export_hash = sha256(EXPORT)

    return {
        "status": "PHYSICAL_GLOBAL_SPECTROSCOPY_READY" if not blockers else "PHYSICAL_GLOBAL_SPECTROSCOPY_BLOCKED",
        "blockers": blockers,
        "direct_certificate_sha256": direct_hash,
        "spectral_export_sha256": export_hash,
    }


def main():
    build = build_onshell_central(ROOT)
    stream = build.direct41.sort_values("x").reset_index(drop=True)
    reducer = module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")

    perL = [local_for_L(stream, reducer, int(L)) for L in DEFAULT_L]
    any_reorder = any(x["rank_reordering_detected"] for x in perL)
    all_controls = all(
        x["common_rescaling_control"]["pass"]
        and x["coarse_grid_same_reordering_boolean"]
        and x["min_K"] > 0
        for x in perL
    )

    local_status = (
        "LOCAL_PRINCIPAL_SPECTRAL_WEIGHT_REORDERING"
        if any_reorder and all_controls
        else "LOCAL_PRINCIPAL_SPECTRAL_WEIGHT_NULL"
        if all_controls
        else "LOCAL_PRINCIPAL_SPECTROSCOPY_CONTROL_FAIL"
    )

    payload = {
        "global_physical_spectroscopy_gate": global_gate(),
        "local_principal_spectroscopy": {
            "status": local_status,
            "scope": (
                "registered healthy production window only; local canonically-normalized "
                "principal characteristic spectroscopy, NOT global QNM/residue spectroscopy"
            ),
            "observable_definition": "W[a,n](u)=|e_a^T v_n(u)|^2 in K-canonical field basis",
            "mode_tracking": "maximum adjacent eigenvector overlap with sign continuity",
            "all_controls_pass": all_controls,
            "any_rank_reordering": any_reorder,
            "per_L": perL,
        },
        "interpretation_guard": (
            "Only a certified direct-global KRGM eigenoperator may support physical QNM "
            "frequencies, Green-function residues Z_n(r), or a Weisz-style global spectral-selection claim."
        ),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    print(json.dumps(payload, indent=2, allow_nan=False))
    return 0 if all_controls else 1


if __name__ == "__main__":
    raise SystemExit(main())
