#!/usr/bin/env python3
"""Production-window spectroscopy diagnostics for the current SSZ P5 member.

This module finally performs the spectroscopy tests that are admissible *before*
a direct center-to-infinity KRGM/QNM certificate exists.

It deliberately separates three levels:

1. local frozen-coefficient spectroscopy on the healthy production window,
2. finite boxed normal-mode spectroscopy on an interior production subwindow,
3. the physical center-to-infinity QNM problem, which remains guarded and is
   NOT bypassed here.

The reduced project convention is

    L = Ydot^T K Ydot - Y'^T G Y' + Y'^T S Y + Y^T M Y,

with S^T=-S.  Hence the frozen-k Hermitian pencil is

    H(k) = k^2 G + i k S - M,    H v = omega^2 K v.

For the boxed problem the corresponding self-adjoint spatial operator is

    H = -d_r(G d_r) + S d_r + (1/2) S' - M,

implemented in weak/variational form.  Dirichlet box boundaries are artificial
and therefore the boxed spectrum is a diagnostic normal-mode spectrum, not QNM.

Outputs are designed to test:
* non-common-redshift frequency structure,
* radial changes in observable/probe spectral weights,
* Weisz-like spectral selectivity (few modes carrying most local weight),
* mode identity tracked by eigenvector overlap,
* contrast against a frozen-coefficient box control.
"""

from __future__ import annotations

import itertools
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy import linalg, sparse

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ssz_p5.config import DEFAULT_L  # noqa:E402
from ssz_p5.numerics import module  # noqa:E402
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central  # noqa:E402
from ssz_p5.qnm.descriptor_pencil import local_poly_differentiation_matrix  # noqa:E402
from ssz_p5.qnm.gate import require_direct_krgm_certificate  # noqa:E402

OUT = ROOT / "data/generated/spectral/PRODUCTION_WINDOW_SPECTROSCOPY_SUITE.json"

LOCAL_LS = (6, 42, 110, 420, 1000)
BOX_LS = tuple(DEFAULT_L)
LOCAL_U = np.linspace(0.628, 0.692, 17)
BOX_U_LO = 0.628
BOX_U_HI = 0.692
BOX_NODES = 121
BOX_KEEP_MODES = 36
PROBES = ("psi", "dphi", "V")


def sym(a):
    return (a + np.swapaxes(a, -1, -2)) / 2.0


def herm(a):
    return (a + np.conjugate(np.swapaxes(a, -1, -2))) / 2.0


def nearest_rows(u, targets):
    return np.array([int(np.argmin(np.abs(u - t))) for t in targets], dtype=int)


def k_whiten(K):
    w, U = np.linalg.eigh(sym(K))
    if np.min(w) <= 0:
        raise ValueError(f"K not positive: min={np.min(w)}")
    inv = U @ np.diag(1.0 / np.sqrt(w)) @ U.T
    return inv, w


def best_overlap_permutation(prev, cur):
    # columns are mode vectors in whitened coordinates
    n = prev.shape[1]
    best = None
    for p in itertools.permutations(range(n)):
        score = sum(abs(np.vdot(prev[:, j], cur[:, p[j]])) for j in range(n))
        if best is None or score > best[0]:
            best = (score, p)
    return list(best[1])


def local_spectrum(K, G, S, M, k):
    Ki, _ = k_whiten(K)
    H = k * k * sym(G) + 1j * k * ((S - S.T) / 2.0) - sym(M)
    Hw = herm(Ki @ H @ Ki)
    omega2, U = np.linalg.eigh(Hw)
    # original-coordinate eigenvectors, K-normalized
    V = Ki @ U
    return omega2.real, U, V, float(np.max(np.abs(Hw - Hw.conj().T)))


def normalized_probe_weights(V, omega2):
    out = {}
    good = omega2 > 0
    om = np.sqrt(np.clip(omega2, 0, None))
    for a, name in enumerate(PROBES):
        raw = np.zeros(len(omega2), float)
        raw[good] = np.abs(V[a, good]) ** 2 / np.maximum(2.0 * om[good], 1e-300)
        s = float(raw.sum())
        p = raw / s if s > 0 else raw
        entropy = float(-np.sum(p[p > 0] * np.log(p[p > 0]))) if s > 0 else None
        out[name] = {
            "weights": [float(x) for x in p],
            "dominant_mode": int(np.argmax(p)) if s > 0 else None,
            "top_fraction": float(np.max(p)) if s > 0 else None,
            "effective_mode_count": float(np.exp(entropy)) if entropy is not None else None,
        }
    return out


def local_suite(d, red, L):
    a = red.canonical_audit(d, int(L))
    K, G, S, M = (np.asarray(a[x]) for x in ("K", "G", "S", "M"))
    u = d.u.to_numpy(float)
    x = d.x.to_numpy(float)
    rows = nearest_rows(u, LOCAL_U)
    dx_box = abs(float(x[rows[-1]] - x[rows[0]]))
    ks = {"k0": 0.0, "k_boxfund": math.pi / max(dx_box, 1e-12)}

    result = {}
    for kname, kval in ks.items():
        records = []
        prev = None
        tracked_frequencies = []
        tracked_probe = {name: [] for name in PROBES}
        tracked_vectors = []
        herm_res = 0.0
        for ii, row in enumerate(rows):
            o2, Uw, V, hres = local_spectrum(K[row], G[row], S[row], M[row], kval)
            herm_res = max(herm_res, hres)
            # initial sort ascending omega2; afterwards track by overlap
            if prev is not None:
                perm = best_overlap_permutation(prev, Uw)
                o2 = o2[perm]
                Uw = Uw[:, perm]
                V = V[:, perm]
            prev = Uw
            tracked_vectors.append(Uw)
            weights = normalized_probe_weights(V, o2)
            om = [float(np.sqrt(v)) if v > 0 else None for v in o2]
            tracked_frequencies.append(om)
            for name in PROBES:
                tracked_probe[name].append(weights[name])
            records.append({
                "u": float(u[row]),
                "x": float(x[row]),
                "omega2": [float(v) for v in o2],
                "omega": om,
                "positive_mode_count": int(np.sum(o2 > 0)),
                "probe": weights,
            })

        # Redshift-only null: all positive tracked frequency ratios would be constant.
        ratio_stats = {}
        nmode = 3
        for p in range(nmode):
            for q in range(p + 1, nmode):
                arr = []
                for om in tracked_frequencies:
                    if om[p] is not None and om[q] is not None and om[q] != 0:
                        arr.append(om[p] / om[q])
                key = f"{p}/{q}"
                if len(arr) >= 2:
                    arr = np.asarray(arr)
                    ratio_stats[key] = {
                        "n": int(len(arr)),
                        "mean": float(np.mean(arr)),
                        "fractional_range": float((np.max(arr) - np.min(arr)) / max(abs(np.mean(arr)), 1e-300)),
                    }
                else:
                    ratio_stats[key] = {"n": int(len(arr)), "mean": None, "fractional_range": None}

        probe_summary = {}
        for name in PROBES:
            dom = [r["dominant_mode"] for r in tracked_probe[name] if r["dominant_mode"] is not None]
            neff = [r["effective_mode_count"] for r in tracked_probe[name] if r["effective_mode_count"] is not None]
            switches = sum(dom[i] != dom[i - 1] for i in range(1, len(dom)))
            probe_summary[name] = {
                "dominant_mode_switches": int(switches),
                "unique_dominant_modes": sorted(set(dom)),
                "mean_effective_mode_count": float(np.mean(neff)) if neff else None,
                "min_effective_mode_count": float(np.min(neff)) if neff else None,
                "max_top_fraction": float(max(r["top_fraction"] for r in tracked_probe[name] if r["top_fraction"] is not None)) if neff else None,
            }

        result[kname] = {
            "k": float(kval),
            "hermiticity_residual_max": herm_res,
            "records": records,
            "frequency_ratio_stats": ratio_stats,
            "probe_summary": probe_summary,
        }
    return result


def trap_weights(x):
    x = np.asarray(x, float)
    w = np.zeros(len(x))
    w[1:-1] = (x[2:] - x[:-2]) / 2.0
    w[0] = (x[1] - x[0]) / 2.0
    w[-1] = (x[-1] - x[-2]) / 2.0
    return w


def blockdiag_weighted(arr, w):
    return sparse.block_diag([w[i] * arr[i] for i in range(len(w))], format="csr")


def boxed_matrices(K, G, S, M, x):
    n = len(x)
    D = local_poly_differentiation_matrix(x, 1, window=9, degree=8)
    D3 = sparse.kron(D, sparse.identity(3, format="csr"), format="csr")
    w = trap_weights(x)
    BK = blockdiag_weighted(sym(K), w)
    WG = blockdiag_weighted(sym(G), w)
    WS = blockdiag_weighted((S - np.swapaxes(S,1,2))/2.0, w)
    WM = blockdiag_weighted(sym(M), w)

    H_G = D3.T @ WG @ D3
    H_S = 0.5 * (WS @ D3 - D3.T @ WS)
    H = (H_G + H_S - WM).tocsr()

    # Dirichlet: remove first/last node, all three fields.
    keep = np.arange(3, 3 * (n - 1), dtype=int)
    H = H[keep][:, keep]
    BK = BK[keep][:, keep]
    # Explicit symmetric cleanup at roundoff level.
    H = ((H + H.T) * 0.5).tocsr()
    BK = ((BK + BK.T) * 0.5).tocsr()
    return H, BK, keep, w


def boxed_spectrum_from_arrays(K, G, S, M, x, u, f):
    H, BK, keep, w = boxed_matrices(K, G, S, M, x)
    Hd = H.toarray()
    Bd = BK.toarray()
    be = np.linalg.eigvalsh(Bd)
    if np.min(be) <= 0:
        raise ValueError(f"boxed K mass not positive: {np.min(be)}")
    vals, vecs = linalg.eigh(Hd, Bd, check_finite=True)
    order = np.argsort(vals)
    vals = vals[order]
    vecs = vecs[:, order]
    pos_idx = np.flatnonzero(vals > 0)[:BOX_KEEP_MODES]
    vals_p = vals[pos_idx]
    vecs_p = vecs[:, pos_idx]

    n = len(x)
    nm = len(vals_p)
    full = np.zeros((3*n, nm), float)
    full[keep, :] = vecs_p

    # node-major vectors -> (node,field,mode)
    Y = full.reshape(n, 3, nm)
    omega = np.sqrt(vals_p)

    # local kinetic density sums to ~1 mode-by-mode.
    kin = np.zeros((n, nm))
    for i in range(n):
        for m in range(nm):
            y = Y[i,:,m]
            kin[i,m] = max(0.0, w[i] * float(y @ sym(K[i]) @ y))

    # selected radii and Weisz-like local spectral selectivity metrics.
    targets = (0.632, 0.648, 0.664, 0.680, 0.690)
    rows = nearest_rows(u, targets)
    local = []
    for i in rows:
        rec = {"u": float(u[i]), "x": float(x[i]), "probes": {}}
        for a, name in enumerate(PROBES):
            raw = np.abs(Y[i,a,:]) ** 2 / np.maximum(2.0*omega, 1e-300)
            p = raw / raw.sum() if raw.sum() > 0 else raw
            ent = -np.sum(p[p>0]*np.log(p[p>0])) if np.any(p>0) else np.nan
            rec["probes"][name] = {
                "dominant_mode": int(np.argmax(p)) if len(p) else None,
                "top_fraction": float(np.max(p)) if len(p) else None,
                "effective_mode_count": float(np.exp(ent)) if np.isfinite(ent) else None,
                "top5": [
                    {"mode": int(j), "omega": float(omega[j]), "weight": float(p[j])}
                    for j in np.argsort(-p)[:5]
                ],
            }
        pk = kin[i] / kin[i].sum() if kin[i].sum() > 0 else kin[i]
        entk = -np.sum(pk[pk>0]*np.log(pk[pk>0])) if np.any(pk>0) else np.nan
        rec["kinetic_density_spectral"] = {
            "dominant_mode": int(np.argmax(pk)) if len(pk) else None,
            "top_fraction": float(np.max(pk)) if len(pk) else None,
            "effective_mode_count": float(np.exp(entk)) if np.isfinite(entk) else None,
        }
        local.append(rec)

    ipr = []
    for m in range(nm):
        p = kin[:,m] / max(kin[:,m].sum(), 1e-300)
        ipr.append(float(np.sum(p*p)))

    return {
        "dimension": int(H.shape[0]),
        "negative_omega2_count": int(np.sum(vals < 0)),
        "zero_or_negative_omega2_min": float(vals[0]),
        "positive_mode_count_total": int(np.sum(vals > 0)),
        "kept_positive_modes": int(nm),
        "omega_first": [float(x) for x in omega[:min(12,nm)]],
        "omega_redshift_div_sqrt_f_mid": [
            float(x / math.sqrt(float(f[len(f)//2]))) for x in omega[:min(12,nm)]
        ],
        "mode_kinetic_ipr": ipr,
        "local_spectral_weights": local,
        "operator_symmetry_residual": float(np.max(np.abs(Hd-Hd.T))),
        "mass_symmetry_residual": float(np.max(np.abs(Bd-Bd.T))),
    }


def boxed_suite(d, red, L):
    # IMPORTANT: reduce K,G,S,M on the full-resolution production grid first.
    # S and M contain profile derivatives; re-reducing an already downsampled
    # coefficient table changes the lower-order operator and is not admissible.
    a = red.canonical_audit(d, int(L))
    Kf, Gf, Sf, Mf = (np.asarray(a[x], float) for x in ("K","G","S","M"))

    uall = d.u.to_numpy(float)
    ids = np.flatnonzero((uall > BOX_U_LO) & (uall < BOX_U_HI))
    pick = np.linspace(0, len(ids)-1, BOX_NODES).round().astype(int)
    ids = ids[np.unique(pick)]

    K, G, S, M = Kf[ids], Gf[ids], Sf[ids], Mf[ids]
    x = d.x.to_numpy(float)[ids]
    u = uall[ids]
    f = d.f.to_numpy(float)[ids]

    physical = boxed_spectrum_from_arrays(K,G,S,M,x,u,f)

    # Frozen-coefficient control: same box and coordinates, matrices frozen at midpoint.
    mid = len(ids)//2
    K0 = np.repeat(K[mid][None,:,:], len(ids), axis=0)
    G0 = np.repeat(G[mid][None,:,:], len(ids), axis=0)
    S0 = np.repeat(S[mid][None,:,:], len(ids), axis=0)
    M0 = np.repeat(M[mid][None,:,:], len(ids), axis=0)
    control = boxed_spectrum_from_arrays(K0,G0,S0,M0,x,u,np.repeat(f[mid],len(ids)))

    # Summary comparison.
    phys_neff = []
    ctrl_neff = []
    for rec in physical["local_spectral_weights"]:
        for name in PROBES:
            phys_neff.append(rec["probes"][name]["effective_mode_count"])
    for rec in control["local_spectral_weights"]:
        for name in PROBES:
            ctrl_neff.append(rec["probes"][name]["effective_mode_count"])

    return {
        "window_u": [float(u.min()), float(u.max())],
        "nodes": int(len(ids)),
        "reduction_semantics": "full-resolution reduction before matrix sampling",
        "physical": physical,
        "frozen_matrix_control": control,
        "comparison": {
            "mean_probe_effective_mode_count_physical": float(np.mean(phys_neff)),
            "mean_probe_effective_mode_count_control": float(np.mean(ctrl_neff)),
            "mean_mode_ipr_physical": float(np.mean(physical["mode_kinetic_ipr"])) if physical["mode_kinetic_ipr"] else None,
            "mean_mode_ipr_control": float(np.mean(control["mode_kinetic_ipr"])) if control["mode_kinetic_ipr"] else None,
            "negative_omega2_count_physical": int(physical["negative_omega2_count"]),
            "negative_omega2_count_control": int(control["negative_omega2_count"]),
        },
    }


def main():
    build = build_onshell_central(ROOT)
    d = build.direct41.sort_values("x").reset_index(drop=True)
    prod = d[(d.u > 0.62) & (d.u < 0.70)].reset_index(drop=True)
    red = module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")

    local = {str(L): local_suite(prod, red, L) for L in LOCAL_LS}
    boxed = {str(L): boxed_suite(prod, red, L) for L in BOX_LS}

    cert = ROOT / "data/certificates/DIRECT_GLOBAL_KRGM_CERTIFICATE.json"
    try:
        require_direct_krgm_certificate(cert, ROOT)
        qnm_admissible = True
        qnm_guard = "VALID_DIRECT_GLOBAL_CERTIFICATE_PRESENT"
    except RuntimeError as exc:
        qnm_admissible = False
        qnm_guard = str(exc)

    # Cross-L compact result flags. These are diagnostic, not theory claims.
    any_local_probe_switch = any(
        rec["probe_summary"][p]["dominant_mode_switches"] > 0
        for L in local.values()
        for rec in L.values()
        for p in PROBES
    )
    max_ratio_variation = 0.0
    for L in local.values():
        for rec in L.values():
            for rr in rec["frequency_ratio_stats"].values():
                if rr["fractional_range"] is not None:
                    max_ratio_variation = max(max_ratio_variation, rr["fractional_range"])

    boxed_selection = {}
    for L, rec in boxed.items():
        p = rec["comparison"]["mean_probe_effective_mode_count_physical"]
        c = rec["comparison"]["mean_probe_effective_mode_count_control"]
        boxed_selection[L] = {
            "physical_mean_Neff": p,
            "control_mean_Neff": c,
            "physical_over_control_Neff": float(p/c) if c else None,
            "physical_mean_IPR": rec["comparison"]["mean_mode_ipr_physical"],
            "control_mean_IPR": rec["comparison"]["mean_mode_ipr_control"],
            "negative_omega2_count_physical": rec["comparison"]["negative_omega2_count_physical"],
            "negative_omega2_count_control": rec["comparison"]["negative_omega2_count_control"],
        }

    report = {
        "scope": {
            "local": "frozen-coefficient healthy production-window spectroscopy",
            "boxed": "full-resolution K,G,S,M reduction followed by finite Dirichlet production-window normal-mode spectroscopy",
            "global_qnm": "not attempted unless direct-global KRGM certificate validates",
        },
        "semantics": {
            "local_pencil": "H(k)=k^2 G + i k S - M; H v = omega^2 K v",
            "boxed_operator": "-d(G d)+S d+0.5 S_prime-M in variational discretization",
            "probe_residue_proxy": "|e_a^T Y_n(r)|^2/(2 omega_n), normalized over retained positive modes",
            "warning": "Local/boxed spectra are diagnostics. They do not replace center-regular/outgoing-infinity QNM residues.",
        },
        "qnm_admissible": qnm_admissible,
        "qnm_guard": qnm_guard,
        "local": local,
        "boxed": boxed,
        "summary": {
            "any_local_probe_dominance_switch": any_local_probe_switch,
            "max_tracked_frequency_ratio_fractional_range": max_ratio_variation,
            "boxed_selection_vs_frozen_control": boxed_selection,
        },
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps({
        "qnm_admissible": qnm_admissible,
        "any_local_probe_dominance_switch": any_local_probe_switch,
        "max_tracked_frequency_ratio_fractional_range": max_ratio_variation,
        "boxed_selection_vs_frozen_control": boxed_selection,
    }, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
