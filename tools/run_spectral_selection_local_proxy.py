#!/usr/bin/env python3
"""Finite-domain canonical spectral-selection probe for the CURRENT frozen P5 member.

This is deliberately NOT a global QNM calculation.  It uses the current
on-shell electric production member on its certified strong-field production
window (0.62 < u < 0.70), derives the existing reduced K,G,S,M matrices, and
asks whether radial operator variation changes canonical mode localization /
local modal composition relative to a constant-coefficient null control.

Outputs are explicitly tagged LOCAL_PROXY and MUST NOT be promoted to a
center-to-infinity physical spectral-selection verdict.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import linalg, sparse

from ssz_p5.numerics import module
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central
from ssz_p5.qnm.descriptor_pencil import local_poly_differentiation_matrix

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/generated/spectral_selection_probe"
LS = (6, 20, 42, 110)
RESOLUTIONS = (81, 121, 161)
NMODES = 18


def trap_weights(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, float)
    w = np.empty_like(x)
    w[1:-1] = 0.5 * (x[2:] - x[:-2])
    w[0] = 0.5 * (x[1] - x[0])
    w[-1] = 0.5 * (x[-1] - x[-2])
    return w


def block_diag_profile(a: np.ndarray) -> sparse.csr_matrix:
    return sparse.block_diag([a[i] for i in range(len(a))], format="csr")


def solve_box(x, K, G, S, M, nmodes=NMODES):
    """Variational finite-box normal modes with Dirichlet endpoint values.

    The action convention is
      L = ydot^T K ydot - y'^T G y' + y'^T S y - y^T M y,
      S^T=-S.
    We discretize the quadratic form directly, symmetrize at roundoff level,
    and solve A v = omega^2 B v.
    """
    n = len(x)
    D = local_poly_differentiation_matrix(x, 1, window=9, degree=8)
    D3 = sparse.kron(D, sparse.identity(3, format="csr"), format="csr")
    w = trap_weights(x)
    W3 = sparse.diags(np.repeat(w, 3), format="csr")
    Kb, Gb, Sb, Mb = map(block_diag_profile, (K, G, S, M))

    A = D3.T @ W3 @ Gb @ D3 - D3.T @ W3 @ Sb + W3 @ Mb
    B = W3 @ Kb
    A = 0.5 * (A + A.T)
    B = 0.5 * (B + B.T)

    keep_nodes = np.arange(1, n - 1)
    keep = np.concatenate([np.arange(3*i, 3*i+3) for i in keep_nodes])
    Ad = A[keep][:, keep].toarray()
    Bd = B[keep][:, keep].toarray()

    bmin = float(np.linalg.eigvalsh(Bd)[0])
    if bmin <= 0:
        raise RuntimeError(f"finite-box B not positive: {bmin}")

    hi = min(len(keep) - 1, max(nmodes * 4, 48))
    evals, evecs = linalg.eigh(Ad, Bd, subset_by_index=(0, hi))
    good = np.flatnonzero(np.isfinite(evals) & (evals > 1e-10))
    if len(good) < nmodes:
        raise RuntimeError(f"only {len(good)} positive finite modes")
    good = good[:nmodes]
    lam = evals[good]
    vec = evecs[:, good]

    full = np.zeros((3*n, nmodes), float)
    full[keep, :] = vec
    fields = full.reshape(n, 3, nmodes)

    density = np.empty((n, nmodes), float)
    for i in range(n):
        # canonical kinetic density; eigenvectors are globally B-normalized
        density[i] = np.einsum("an,ab,bn->n", fields[i], K[i], fields[i])
    mass = np.maximum(density * w[:, None], 0.0)
    norm = mass.sum(axis=0)
    mass = mass / np.maximum(norm[None, :], 1e-300)
    ipr = np.sum(mass**2, axis=0)

    # Local modal composition, using canonical kinetic density as a
    # basis-invariant positive local weight on this finite box.
    local = np.maximum(density, 0.0)
    probs = local / np.maximum(local.sum(axis=1, keepdims=True), 1e-300)
    dom = np.argmax(probs, axis=1)
    omega = np.sqrt(lam)
    centroid = probs @ omega

    return {
        "omega": omega,
        "ipr": ipr,
        "density": density,
        "probs": probs,
        "dominant": dom,
        "centroid": centroid,
        "bmin": bmin,
    }


def constant_control(K, G, S, M):
    """Freeze all radial matrices at the median interior profile value."""
    n = len(K)
    mid = slice(max(1, n//4), min(n-1, 3*n//4))
    return tuple(np.repeat(np.median(a[mid], axis=0)[None, :, :], n, axis=0)
                 for a in (K, G, S, M))


def js_div(p, q):
    p = np.asarray(p, float)
    q = np.asarray(q, float)
    p = p / p.sum()
    q = q / q.sum()
    m = 0.5 * (p + q)
    def kl(a, b):
        z = a > 0
        return float(np.sum(a[z] * np.log(a[z] / b[z])))
    return 0.5 * kl(p, m) + 0.5 * kl(q, m)


def summarize(x, sol, ctl):
    n = len(x)
    lo, hi = int(0.15*n), int(0.85*n)
    dom = sol["dominant"][lo:hi]
    domc = ctl["dominant"][lo:hi]
    trans = int(np.sum(dom[1:] != dom[:-1]))
    transc = int(np.sum(domc[1:] != domc[:-1]))
    q1, q3 = int(0.30*n), int(0.70*n)
    js = js_div(sol["probs"][q1], sol["probs"][q3])
    jsc = js_div(ctl["probs"][q1], ctl["probs"][q3])
    ipr_med = float(np.median(sol["ipr"][:12]))
    ipr_ctl = float(np.median(ctl["ipr"][:12]))
    return {
        "omega_first6": [float(v) for v in sol["omega"][:6]],
        "median_ipr_first12": ipr_med,
        "control_median_ipr_first12": ipr_ctl,
        "ipr_ratio": float(ipr_med / max(ipr_ctl, 1e-300)),
        "interior_dominance_transitions": trans,
        "control_interior_dominance_transitions": transc,
        "js_q30_q70": float(js),
        "control_js_q30_q70": float(jsc),
        "js_excess": float(js - jsc),
        "centroid_relative_span": float(
            (np.max(sol["centroid"][lo:hi]) - np.min(sol["centroid"][lo:hi]))
            / max(np.mean(sol["centroid"][lo:hi]), 1e-300)
        ),
        "control_centroid_relative_span": float(
            (np.max(ctl["centroid"][lo:hi]) - np.min(ctl["centroid"][lo:hi]))
            / max(np.mean(ctl["centroid"][lo:hi]), 1e-300)
        ),
        "B_min_eigenvalue": float(sol["bmin"]),
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    build = build_onshell_central(ROOT)
    d = build.direct41.reset_index(drop=True)
    prod = np.flatnonzero((d.u.to_numpy(float) > 0.62) & (d.u.to_numpy(float) < 0.70))
    reducer = module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")

    report = {
        "schema_version": "1.0",
        "status": "LOCAL_PROXY_ONLY",
        "global_physical_verdict": "NOT_EVALUATED",
        "reason_global_not_evaluated": (
            "current frozen electric member is certified only on the strong-field "
            "production window; no same-member center-to-infinity canonical operator "
            "is supplied"
        ),
        "member_hash": "8bd460ef022a9cdbcc3644abd8aecbfbb910f8e364ac1410378d2641291559cf",
        "observable": "canonical kinetic local mode density on finite strong-field box",
        "boundary_condition": "Dirichlet at finite production-window endpoints",
        "controls": "same grid/boundaries with K,G,S,M frozen to median interior matrices",
        "runs": [],
    }
    mode_rows = []

    for L in LS:
        audit = reducer.canonical_audit(d, float(L))
        K0, G0, S0, M0 = (np.asarray(audit[k], float) for k in ("K","G","S","M"))
        for N in RESOLUTIONS:
            choose = prod[np.unique(np.linspace(0, len(prod)-1, N, dtype=int))]
            # x in source is decreasing with u; spectral differentiator requires increasing x.
            choose = choose[np.argsort(d.x.to_numpy(float)[choose])]
            x = d.x.to_numpy(float)[choose]
            K, G, S, M = (a[choose] for a in (K0,G0,S0,M0))
            sol = solve_box(x, K, G, S, M)
            Kc,Gc,Sc,Mc = constant_control(K,G,S,M)
            ctl = solve_box(x, Kc,Gc,Sc,Mc)
            s = summarize(x, sol, ctl)
            s.update({"L": int(L), "N": int(len(x)), "x_min": float(x.min()), "x_max": float(x.max())})
            report["runs"].append(s)
            for n in range(min(12, len(sol["omega"]))):
                mode_rows.append({
                    "L": L, "N": len(x), "mode": n,
                    "omega": float(sol["omega"][n]),
                    "ipr": float(sol["ipr"][n]),
                    "control_ipr": float(ctl["ipr"][n]),
                })

    # Convergence summaries: compare highest two radial resolutions.
    conv = []
    for L in LS:
        rr = [r for r in report["runs"] if r["L"] == L]
        rr = sorted(rr, key=lambda z: z["N"])
        a, b = rr[-2], rr[-1]
        oa, ob = np.array(a["omega_first6"]), np.array(b["omega_first6"])
        conv.append({
            "L": L,
            "max_rel_omega_change_last_two": float(np.max(np.abs(ob-oa)/np.maximum(np.abs(ob),1e-12))),
            "ipr_ratio_change": float(abs(b["ipr_ratio"]-a["ipr_ratio"])),
            "js_excess_change": float(abs(b["js_excess"]-a["js_excess"])),
        })
    report["convergence"] = conv

    # Deliberately conservative local-proxy verdict.
    hi = [r for r in report["runs"] if r["N"] == max(RESOLUTIONS)]
    positive = all(
        r["ipr_ratio"] > 1.10
        and r["js_excess"] > 0.0
        and r["centroid_relative_span"] > r["control_centroid_relative_span"]
        for r in hi
    )
    converged = all(c["max_rel_omega_change_last_two"] < 0.08 for c in conv)
    report["local_proxy_verdict"] = (
        "LOCAL_PROXY_SIGNAL" if positive and converged else "LOCAL_PROXY_NULL_OR_INCONCLUSIVE"
    )
    report["interpretation_guard"] = (
        "This finite-window result cannot establish physical SSZ spectral-weight "
        "permutation. It is a falsification-oriented precursor only; global regular-center "
        "and outgoing/asymptotic boundary conditions remain required."
    )

    (OUT / "SPECTRAL_SELECTION_LOCAL_PROXY.json").write_text(json.dumps(report, indent=2) + "\n")
    pd.DataFrame(mode_rows).to_csv(OUT / "SPECTRAL_SELECTION_LOCAL_MODES.csv", index=False)
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
