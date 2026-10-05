#!/usr/bin/env python3
"""SSZ spectroscopy v2: fixed-observable generalized-eigenmode residue audit.

This is the strongest spectroscopy test currently admissible without a
center-to-infinity direct-global KRGM certificate.

Inside the registered healthy production window only, it solves

    G(u) v_n(u) = omega_n(u)^2 K(u) v_n(u)

with K-normalized eigenvectors v_n^T K v_m = delta_nm.  For a fixed original
field-space observable covector q, the local pole-residue proxy is

    Z_n^q(u) = |q^T v_n(u)|^2 / (2 omega_n(u)).

The normalized weights Z/sum(Z) are tracked continuously in radius.  This is
basis-covariant under constant invertible field redefinitions and much closer
to the Green-function residue question than canonical-basis component weights.

It remains a LOCAL PRINCIPAL proxy, not a global QNM spectrum.  Global physical
spectroscopy is hard-gated on full health + direct-global KRGM export.
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

OUT = ROOT / "data/generated/spectral/SSZ_SPECTROSCOPY_FIXED_OBSERVABLE_V2.json"
DIRECT_CERTS = (
    ROOT / "data/certificates/SSZ_P5_DIRECT_GLOBAL_KRGM_CERTIFICATE.json",
    ROOT / "data/certificates/DIRECT_GLOBAL_KRGM_CERTIFICATE.json",
)
GLOBAL_EXPORT = ROOT / "data/generated/spectral/GLOBAL_CANONICAL_KRGM_SPECTRAL_EXPORT.npz"
FULL_AUDIT = ROOT / "data/generated/spectral/CURRENT_MEMBER_FULL_AVAILABLE_DOMAIN_AUDIT.json"


def sym(a):
    return (a + np.swapaxes(a, -1, -2)) / 2.0


def solve_modes(K, G):
    K = (K + K.T) / 2.0
    G = (G + G.T) / 2.0
    kw = np.linalg.eigvalsh(K)
    if np.min(kw) <= 0:
        raise ValueError("K_NOT_POSITIVE")
    lam, V = eigh(G, K, check_finite=True)
    if np.min(lam) <= 0:
        raise ValueError("RADIAL_PRINCIPAL_NOT_POSITIVE")
    return np.sqrt(lam), V


def metric_overlap(Va, Vb, Ka, Kb):
    Kmid = (Ka + Kb) / 2.0
    M = Va.T @ Kmid @ Vb
    aa = np.sqrt(np.maximum(np.diag(Va.T @ Kmid @ Va), 1e-300))
    bb = np.sqrt(np.maximum(np.diag(Vb.T @ Kmid @ Vb), 1e-300))
    return M / (aa[:, None] * bb[None, :])


def best_perm(prevV, curV, prevK, curK):
    obs_o = np.abs(metric_overlap(prevV, curV, prevK, curK))
    n = obs_o.shape[0]
    best = None
    for p in itertools.permutations(range(n)):
        score = float(sum(obs_o[i, p[i]] for i in range(n)))
        if best is None or score > best[0]:
            best = (score, p)
    return list(best[1]), float(best[0] / n)


def track(omegas, vecs, Ks):
    to = [omegas[0].copy()]
    tv = [vecs[0].copy()]
    minovs = [1.0]
    perms = [list(range(vecs[0].shape[1]))]
    for i in range(1, len(omegas)):
        p, sc = best_perm(tv[-1], vecs[i], Ks[i - 1], Ks[i])
        v = vecs[i][:, p].copy()
        o = omegas[i][p].copy()
        obs_o = metric_overlap(tv[-1], v, Ks[i - 1], Ks[i])
        for j in range(v.shape[1]):
            if obs_o[j, j] < 0:
                v[:, j] *= -1.0
        to.append(o)
        tv.append(v)
        minovs.append(float(np.min(np.abs(np.diag(metric_overlap(tv[-2], v, Ks[i - 1], Ks[i]))))))
        perms.append(p)
    return np.asarray(to), np.asarray(tv), np.asarray(minovs), perms


def residue_weights(omega, V, q):
    raw = np.abs(q @ V) ** 2 / np.maximum(2.0 * omega, 1e-300)
    return raw / max(float(np.sum(raw)), 1e-300), raw


def rank_changes(W, u):
    ranks = [tuple(np.argsort(-w).tolist()) for w in W]
    changes = []
    for i in range(1, len(ranks)):
        if ranks[i] != ranks[i - 1]:
            changes.append({
                "u": float(u[i]),
                "from": list(ranks[i - 1]),
                "to": list(ranks[i]),
            })
    return ranks, changes


def ratio_spans(omega):
    out = {}
    n = omega.shape[1]
    for i in range(n):
        for j in range(i + 1, n):
            r = omega[:, i] / omega[:, j]
            out[f"{i}/{j}"] = {
                "min": float(np.min(r)),
                "max": float(np.max(r)),
                "relative_span": float((np.max(r) - np.min(r)) / max(np.mean(r), 1e-300)),
            }
    return out


def generalized_resolvent_peak_proxy(K, G, q, omega, eta_frac=0.01):
    """Evaluate -Im q^T [G-(w+i eta)^2 K]^-1 q / pi at modal w."""
    scale = float(np.median(omega))
    eta = max(eta_frac * scale, 1e-12)
    vals = []
    for w in omega:
        z = complex(float(w), eta)
        A = G.astype(complex) - z * z * K.astype(complex)
        val = q @ np.linalg.solve(A, q)
        vals.append(float(-np.imag(val) / np.pi))
    return vals


def build_observables(n):
    obs = {f"field_{i}": np.eye(n)[i] for i in range(n)}
    obs["equal_weight"] = np.ones(n) / np.sqrt(n)
    # A deterministic mixed channel that is not symmetry-aligned.
    q = np.arange(1, n + 1, dtype=float)
    obs["ramp_weight"] = q / np.linalg.norm(q)
    return obs


def basis_control(Ks, Gs, qs, uidx):
    n = Ks[0].shape[0]
    A = np.eye(n)
    for i in range(n):
        for j in range(i + 1, n):
            A[i, j] = 0.17 / (1 + j - i)
    # x = A x'; K'=A^T K A, G'=A^T G A, q'=A^T q
    max_freq = 0.0
    max_weight = 0.0
    for k in uidx:
        w0, V0 = solve_modes(Ks[k], Gs[k])
        Kp = A.T @ Ks[k] @ A
        Gp = A.T @ Gs[k] @ A
        w1, V1 = solve_modes(Kp, Gp)
        max_freq = max(max_freq, float(np.max(np.abs(w0 - w1))))
        for q in qs.values():
            q1 = A.T @ q
            W0, _ = residue_weights(w0, V0, q)
            W1, _ = residue_weights(w1, V1, q1)
            max_weight = max(max_weight, float(np.max(np.abs(W0 - W1))))
    return {
        "max_abs_frequency_delta": max_freq,
        "max_abs_normalized_residue_delta": max_weight,
        "pass": bool(max_freq < 1e-8 and max_weight < 2e-7),
    }


def common_scale_control(K, G, obs, u):
    s = 0.8 + 0.4 * (u - u.min()) / max(float(u.max() - u.min()), 1e-15)
    maxf = 0.0
    maxw = 0.0
    for k in range(len(u)):
        w0, V0 = solve_modes(K[k], G[k])
        w1, V1 = solve_modes(s[k] * K[k], s[k] * G[k])
        maxf = max(maxf, float(np.max(np.abs(w0 - w1))))
        for q in obs.values():
            W0, _ = residue_weights(w0, V0, q)
            W1, _ = residue_weights(w1, V1, q)
            maxw = max(maxw, float(np.max(np.abs(W0 - W1))))
    return {
        "max_abs_frequency_delta": maxf,
        "max_abs_normalized_residue_delta": maxw,
        "pass": bool(maxf < 1e-8 and maxw < 2e-7),
    }


def global_gate():
    blockers = []
    if FULL_AUDIT.exists():
        try:
            d = json.loads(FULL_AUDIT.read_text())
            if not d.get("full_available_domain_K_pass", False):
                blockers.append("FULL_NATIVE_DOMAIN_K_FAIL")
            if not d.get("full_available_domain_radial_pass", False):
                blockers.append("FULL_NATIVE_DOMAIN_RADIAL_FAIL")
        except Exception:
            blockers.append("INVALID_FULL_NATIVE_DOMAIN_AUDIT")
    else:
        blockers.append("MISSING_FULL_NATIVE_DOMAIN_AUDIT")
    if not any(p.exists() for p in DIRECT_CERTS):
        blockers.append("MISSING_DIRECT_GLOBAL_KRGM_CERTIFICATE")
    if not GLOBAL_EXPORT.exists():
        blockers.append("MISSING_GLOBAL_CANONICAL_KRGM_SPECTRAL_EXPORT")
    return {
        "status": "READY" if not blockers else "BLOCKED",
        "blockers": blockers,
    }


def per_L(stream, reducer, L):
    a = reducer.canonical_audit(stream, int(L))
    Kall = sym(np.asarray(a["K"], float))
    Gall = sym(np.asarray(a["G"], float))
    uall = stream.u.to_numpy(float)
    idx = np.flatnonzero((uall > 0.62) & (uall < 0.70))
    idx = idx[np.unique(np.linspace(0, len(idx) - 1, min(121, len(idx))).round().astype(int))]
    u = uall[idx]
    Ks = Kall[idx]
    Gs = Gall[idx]

    omegas = []
    vecs = []
    for K, G in zip(Ks, Gs):
        w, V = solve_modes(K, G)
        omegas.append(w)
        vecs.append(V)
    omegas = np.asarray(omegas)
    vecs = np.asarray(vecs)
    tom, tv, overlaps, perms = track(omegas, vecs, Ks)

    n = tom.shape[1]
    obs = build_observables(n)
    ores = {}
    any_reorder = False

    for name, q in obs.items():
        W = []
        raw = []
        peaks = []
        for k in range(len(u)):
            ww, rr = residue_weights(tom[k], tv[k], q)
            W.append(ww)
            raw.append(rr)
            peaks.append(generalized_resolvent_peak_proxy(Ks[k], Gs[k], q, tom[k]))
        W = np.asarray(W)
        raw = np.asarray(raw)
        peaks = np.asarray(peaks)
        ranks, changes = rank_changes(W, u)
        any_reorder |= bool(changes)
        entropy = -np.sum(np.clip(W, 1e-300, 1) * np.log(np.clip(W, 1e-300, 1)), axis=1)
        ores[name] = {
            "ranking_change_count": len(changes),
            "ranking_changes": changes[:30],
            "initial_ranking": list(ranks[0]),
            "final_ranking": list(ranks[-1]),
            "normalized_residue_min_by_mode": [float(x) for x in np.min(W, axis=0)],
            "normalized_residue_max_by_mode": [float(x) for x in np.max(W, axis=0)],
            "raw_residue_dynamic_range_by_mode": [
                float(np.max(raw[:, m]) / max(np.min(raw[:, m]), 1e-300))
                for m in range(n)
            ],
            "spectral_entropy_min": float(np.min(entropy)),
            "spectral_entropy_max": float(np.max(entropy)),
            "resolvent_peak_proxy_min_by_mode": [float(x) for x in np.min(peaks, axis=0)],
            "resolvent_peak_proxy_max_by_mode": [float(x) for x in np.max(peaks, axis=0)],
        }

    # Controls on every 10th radius.
    cidx = list(range(0, len(u), max(1, len(u) // 12)))
    bc = basis_control(Ks, Gs, obs, cidx)
    sc = common_scale_control(Ks[cidx], Gs[cidx], obs, u[cidx])

    # Coarse-grid re-run of rankings from already tracked physical residues.
    coarse_same = True
    for name, q in obs.items():
        W = np.asarray([residue_weights(tom[k], tv[k], q)[0] for k in range(len(u))])
        fine = bool(rank_changes(W, u)[1])
        coarse = bool(rank_changes(W[::2], u[::2])[1])
        if fine != coarse:
            coarse_same = False

    ratios = ratio_spans(tom)
    noncommon = any(v["relative_span"] > 1e-3 for v in ratios.values())

    return {
        "L": int(L),
        "u_min": float(u.min()),
        "u_max": float(u.max()),
        "samples": len(u),
        "mode_count": n,
        "min_tracking_overlap": float(np.min(overlaps)),
        "any_fixed_observable_residue_reordering": any_reorder,
        "non_common_frequency_scaling_detected_1e3": noncommon,
        "frequency_ratio_stats": ratios,
        "observables": ores,
        "controls": {
            "constant_field_basis_change": bc,
            "common_KG_scale": sc,
            "coarse_grid_same_reordering_boolean": coarse_same,
            "pass": bool(bc["pass"] and sc["pass"] and coarse_same and np.min(overlaps) > 0.9),
        },
    }


def main():
    build = build_onshell_central(ROOT)
    stream = build.direct41.sort_values("x").reset_index(drop=True)
    reducer = module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")

    per = {}
    for L in DEFAULT_L:
        per[str(int(L))] = per_L(stream, reducer, int(L))

    controls = all(v["controls"]["pass"] for v in per.values())
    any_reorder = any(v["any_fixed_observable_residue_reordering"] for v in per.values())
    any_noncommon = any(v["non_common_frequency_scaling_detected_1e3"] for v in per.values())

    if not controls:
        status = "LOCAL_FIXED_OBSERVABLE_SPECTROSCOPY_CONTROL_FAIL"
    elif any_reorder:
        status = "LOCAL_FIXED_OBSERVABLE_RESIDUE_REORDERING_DETECTED"
    else:
        status = "LOCAL_FIXED_OBSERVABLE_RESIDUE_REORDERING_NULL"

    report = {
        "global_physical_spectroscopy": global_gate(),
        "local_fixed_observable_principal_spectroscopy": {
            "status": status,
            "scope": "registered healthy production window only; generalized local principal operator",
            "residue_definition": "Z_n^q=|q^T v_n|^2/(2 omega_n), with G v=omega^2 K v and v^T K v=1",
            "fixed_observables": "original field-space basis covectors plus equal/ramp mixtures",
            "all_controls_pass": controls,
            "any_residue_ranking_reordering": any_reorder,
            "any_non_common_frequency_scaling": any_noncommon,
            "per_L": per,
        },
        "interpretation": {
            "allowed": (
                "radius-dependent local principal mode character and fixed-observable residue-proxy hierarchy "
                "inside the certified production window"
            ),
            "not_allowed_yet": (
                "global QNM spectrum, physical retarded Green-function poles/residues, or a global "
                "Weisz-style spectral-selection claim"
            ),
        },
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps({
        "global": report["global_physical_spectroscopy"],
        "local_status": status,
        "controls": controls,
        "any_reorder": any_reorder,
        "any_noncommon": any_noncommon,
        "per_L_summary": {
            L: {
                "min_overlap": v["min_tracking_overlap"],
                "reorder": v["any_fixed_observable_residue_reordering"],
                "noncommon": v["non_common_frequency_scaling_detected_1e3"],
                "controls": v["controls"]["pass"],
            }
            for L, v in per.items()
        },
    }, indent=2, allow_nan=False))
    return 0 if controls else 1


if __name__ == "__main__":
    raise SystemExit(main())
