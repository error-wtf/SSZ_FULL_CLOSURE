#!/usr/bin/env python3
"""Local strong-field spectral-selection audit for the frozen current SSZ P5 member.

This is deliberately NOT a global QNM certificate.  It consumes the current
same-action Central member and the accepted profile-aware reducer to obtain the
full local K,R,G,S,M operator on 0.62<u<0.70, then performs two independent
checks:

1. local principal-symbol mode tracking from the generalized (G,K) problem;
2. a finite-window self-adjoint generalized eigenproblem using the complete
   reduced K,G,S,M matrices with Dirichlet endpoints.

The finite-window problem is a falsification/diagnostic proxy for radial
spectral reordering.  Its boundary conditions are artificial and therefore its
frequencies MUST NOT be published as physical QNMs.

Outputs:
  build/SPECTRAL_SELECTION_LOCAL_PROXY.json
  build/SPECTRAL_SELECTION_LOCAL_PROXY.md
"""
from __future__ import annotations

import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.linalg import eigh
from scipy.sparse import coo_matrix, csc_matrix
from scipy.sparse.linalg import eigsh

from ssz_p5.numerics import module
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central

ROOT = Path(__file__).resolve().parents[1]
OUT_JSON = ROOT / "build/SPECTRAL_SELECTION_LOCAL_PROXY.json"
OUT_MD = ROOT / "build/SPECTRAL_SELECTION_LOCAL_PROXY.md"
L_VALUES = (6, 42, 1000)
RESOLUTIONS = (100, 160, 220)
NMODES = 12
POSITIVE_MODES = 6
U_LO = 0.62
U_HI = 0.70


def sym(a):
    return 0.5 * (a + np.swapaxes(a, -1, -2))


def matrix_sqrt(a):
    w, q = np.linalg.eigh(sym(a))
    if np.min(w) <= 0:
        raise ValueError(f"non-positive K in matrix_sqrt: {np.min(w)}")
    return (q * np.sqrt(w)) @ q.T


def matrix_invsqrt(a):
    w, q = np.linalg.eigh(sym(a))
    if np.min(w) <= 0:
        raise ValueError(f"non-positive K in matrix_invsqrt: {np.min(w)}")
    return (q * (1.0 / np.sqrt(w))) @ q.T


def best_perm(prev, cur):
    best = None
    best_score = -1.0
    for p in itertools.permutations(range(3)):
        score = sum(abs(float(prev[:, i] @ cur[:, p[i]])) for i in range(3))
        if score > best_score:
            best_score = score
            best = p
    return best


def principal_track(r, K, G):
    n = len(r)
    vals = np.empty((n, 3))
    vecs = np.empty((n, 3, 3))
    for i in range(n):
        ki = matrix_invsqrt(K[i])
        c = sym(ki @ sym(G[i]) @ ki)
        w, q = np.linalg.eigh(c)
        vals[i] = w
        vecs[i] = q

    # Track physical branches by overlap, never by eigenvalue order alone.
    for i in range(1, n):
        p = best_perm(vecs[i - 1], vecs[i])
        vals[i] = vals[i, list(p)]
        vecs[i] = vecs[i][:, list(p)]
        for j in range(3):
            if float(vecs[i - 1, :, j] @ vecs[i, :, j]) < 0:
                vecs[i, :, j] *= -1.0

    branches = []
    for j in range(3):
        q0 = vecs[0, :, j]
        ref_angles = []
        adj_angles = []
        comps = []
        for i in range(n):
            dot = min(1.0, max(-1.0, abs(float(q0 @ vecs[i, :, j]))))
            ref_angles.append(math.degrees(math.acos(dot)))
            comps.append(int(np.argmax(np.abs(vecs[i, :, j]))))
            if i:
                d = min(
                    1.0,
                    max(-1.0, abs(float(vecs[i - 1, :, j] @ vecs[i, :, j]))),
                )
                adj_angles.append(math.degrees(math.acos(d)))
        switches = sum(a != b for a, b in zip(comps[:-1], comps[1:]))
        branches.append(
            {
                "branch": j,
                "c2_min": float(np.min(vals[:, j])),
                "c2_max": float(np.max(vals[:, j])),
                "max_rotation_from_outer_deg": float(max(ref_angles)),
                "total_adjacent_rotation_deg": float(sum(adj_angles)),
                "max_adjacent_rotation_deg": float(max(adj_angles) if adj_angles else 0.0),
                "dominant_component_switches": int(switches),
                "dominant_components_seen": sorted(set(comps)),
            }
        )
    return vals, vecs, branches


def build_interpolators(r, arrays):
    return [
        [
            [PchipInterpolator(r, arrays[:, i, j], extrapolate=False) for j in range(3)]
            for i in range(3)
        ]
        for arrays in arrays
    ]


def eval_matrix(interp, x):
    return np.array([[interp[i][j](x) for j in range(3)] for i in range(3)], float)


def assemble_fem(r_src, K_src, G_src, S_src, M_src, n_nodes, grid_kind="r", frozen=False):
    rmin, rmax = float(r_src[0]), float(r_src[-1])
    if grid_kind == "r":
        r = np.linspace(rmin, rmax, n_nodes)
    elif grid_kind == "u":
        u = np.linspace(1.0 / rmin, 1.0 / rmax, n_nodes)
        r = 1.0 / u
    else:
        raise ValueError(grid_kind)

    interps = build_interpolators(r_src, [K_src, G_src, S_src, M_src])
    if frozen:
        rmid = 0.5 * (rmin + rmax)
        frozen_mats = [eval_matrix(ip, rmid) for ip in interps]

    rows = []
    cols = []
    avals = []
    bvals = []

    def add_block(ni, nj, Ablock, Bblock):
        for a in range(3):
            for b in range(3):
                rows.append(3 * ni + a)
                cols.append(3 * nj + b)
                avals.append(float(Ablock[a, b]))
                bvals.append(float(Bblock[a, b]))

    for e in range(n_nodes - 1):
        h = float(r[e + 1] - r[e])
        rm = 0.5 * (r[e + 1] + r[e])
        if frozen:
            Ke, Ge, Se, Me = frozen_mats
        else:
            Ke, Ge, Se, Me = [eval_matrix(ip, rm) for ip in interps]
        Ke, Ge, Me = sym(Ke), sym(Ge), sym(Me)
        Se = 0.5 * (Se - Se.T)
        shape_mass = (h / 6.0) * np.array([[2.0, 1.0], [1.0, 2.0]])
        shape_stiff = (1.0 / h) * np.array([[1.0, -1.0], [-1.0, 1.0]])
        shape_s = np.array([[0.0, 0.5], [-0.5, 0.0]])
        nodes = (e, e + 1)
        for p in range(2):
            for q in range(2):
                ab = (
                    shape_stiff[p, q] * Ge
                    + shape_mass[p, q] * Me
                    + shape_s[p, q] * Se
                )
                bb = shape_mass[p, q] * Ke
                add_block(nodes[p], nodes[q], ab, bb)

    ndof = 3 * n_nodes
    A = coo_matrix((avals, (rows, cols)), shape=(ndof, ndof)).tocsc()
    B = coo_matrix((bvals, (rows, cols)), shape=(ndof, ndof)).tocsc()
    keep = np.arange(3, ndof - 3)
    A = csc_matrix(0.5 * (A[keep][:, keep] + A[keep][:, keep].T))
    B = csc_matrix(0.5 * (B[keep][:, keep] + B[keep][:, keep].T))

    k = min(NMODES, A.shape[0] - 2)
    try:
        w, v = eigsh(A, k=k, M=B, sigma=0.0, which="LM", tol=1e-9, maxiter=50000)
        order = np.argsort(w)
        w, v = w[order], v[:, order]
    except Exception:
        # Deterministic dense fallback for CI and small matrices.
        wd, vd = eigh(A.toarray(), B.toarray(), check_finite=True)
        idx = np.argsort(np.abs(wd))[:k]
        idx = idx[np.argsort(wd[idx])]
        w, v = wd[idx], vd[:, idx]

    full = np.zeros((ndof, len(w)))
    full[keep] = v

    K_nodes = np.empty((n_nodes, 3, 3))
    for i, rv in enumerate(r):
        if frozen:
            K_nodes[i] = frozen_mats[0]
        else:
            K_nodes[i] = sym(eval_matrix(interps[0], rv))

    return r, w, full.reshape(n_nodes, 3, len(w)), K_nodes


def mode_diagnostics(r, omega2, modes, K_nodes):
    pos = np.flatnonzero(omega2 > 0)
    use = pos[:POSITIVE_MODES]
    if len(use) < 2:
        return {
            "positive_modes_found": int(len(pos)),
            "usable_modes": int(len(use)),
            "rank_reordering": {},
            "modes": [],
        }

    z = np.empty((len(r), len(use)))
    canonical = np.empty((len(r), 3, len(use)))
    for i in range(len(r)):
        ks = matrix_sqrt(K_nodes[i])
        for j, n in enumerate(use):
            y = modes[i, :, n]
            q = ks @ y
            canonical[i, :, j] = q
            z[i, j] = float(q @ q)

    for j in range(len(use)):
        norm = float(np.trapezoid(z[:, j], r))
        if norm > 0:
            z[:, j] /= norm
            canonical[:, :, j] /= math.sqrt(norm)

    mode_rows = []
    for j, n in enumerate(use):
        zz = z[:, j]
        ipr = float(np.trapezoid(zz * zz, r))
        imax = int(np.argmax(zz))
        ref = canonical[imax, :, j]
        ref /= max(np.linalg.norm(ref), 1e-300)
        mask = zz > 0.10 * float(np.max(zz))
        angles = []
        comps = []
        for i in np.flatnonzero(mask):
            q = canonical[i, :, j]
            nq = np.linalg.norm(q)
            if nq <= 1e-14:
                continue
            q = q / nq
            angles.append(
                math.degrees(
                    math.acos(min(1.0, max(-1.0, abs(float(ref @ q)))))
                )
            )
            comps.append(int(np.argmax(np.abs(q))))
        switches = sum(a != b for a, b in zip(comps[:-1], comps[1:]))
        mode_rows.append(
            {
                "mode_index": int(n),
                "omega2": float(omega2[n]),
                "omega": float(math.sqrt(omega2[n])),
                "ipr_r": ipr,
                "max_canonical_polarization_rotation_deg": float(max(angles) if angles else 0.0),
                "dominant_component_switches": int(switches),
                "dominant_components_seen": sorted(set(comps)),
            }
        )

    # Weisz-like source residues for each basis source, using the kinetic inner product.
    rank_reordering = {}
    third = max(2, len(r) // 3)
    inner_slice = slice(0, third)
    outer_slice = slice(len(r) - third, len(r))
    for source in range(3):
        residues = []
        for j in range(len(use)):
            # e_source^T K y = e_source^T K^(1/2) q
            vals = []
            for i in range(len(r)):
                y = modes[i, :, use[j]]
                vals.append(float(K_nodes[i][source] @ y))
            amp = float(np.trapezoid(vals, r))
            residues.append(amp * amp)
        residues = np.asarray(residues)
        W = z * residues[None, :]
        inner = np.trapezoid(W[inner_slice], r[inner_slice], axis=0)
        outer = np.trapezoid(W[outer_slice], r[outer_slice], axis=0)
        inner_top = int(use[int(np.argmax(inner))])
        outer_top = int(use[int(np.argmax(outer))])
        rank_reordering[str(source)] = {
            "inner_top_mode": inner_top,
            "outer_top_mode": outer_top,
            "top_mode_changes": bool(inner_top != outer_top),
            "residues": [float(x) for x in residues],
        }

    return {
        "positive_modes_found": int(len(pos)),
        "usable_modes": int(len(use)),
        "rank_reordering": rank_reordering,
        "modes": mode_rows,
    }


def summarize_run(run):
    return [
        float(x)
        for x in run["omega"][: min(POSITIVE_MODES, len(run["omega"]))]
    ]


def main():
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    build = build_onshell_central(ROOT)
    d = build.direct41.reset_index(drop=True)
    red = module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")

    u_all = d.u.to_numpy(float)
    mask = (u_all > U_LO) & (u_all < U_HI)
    idx = np.flatnonzero(mask)
    r_raw = d.x.to_numpy(float)[idx]
    order = np.argsort(r_raw)
    r = r_raw[order]

    report = {
        "schema_version": "1.0",
        "scope": "current-member local strong-field full-KRGSM spectral proxy",
        "member_hash": "8bd460ef022a9cdbcc3644abd8aecbfbb910f8e364ac1410378d2641291559cf",
        "domain_u": [U_LO, U_HI],
        "boundary_conditions": "Dirichlet at finite-window endpoints; diagnostic only, NOT physical QNM",
        "L_values": list(L_VALUES),
        "resolutions": list(RESOLUTIONS),
        "principal": {},
        "finite_window": {},
        "global_physical_verdict": "NOT_EVALUATED_BY_THIS_LOCAL_PROXY",
    }

    principal_signal = False
    cavity_signal = False

    for L in L_VALUES:
        audit = red.canonical_audit(d, float(L))
        K = np.asarray(audit["K"], float)[idx][order]
        G = np.asarray(audit["G"], float)[idx][order]
        S = np.asarray(audit["S"], float)[idx][order]
        M = np.asarray(audit["M"], float)[idx][order]
        R = np.asarray(audit["R"], float)[idx][order]

        pvals, pvecs, pbranches = principal_track(r, K, G)
        _ = pvecs
        p_signal = any(
            b["max_rotation_from_outer_deg"] > 10.0
            and b["dominant_component_switches"] > 0
            for b in pbranches
        )
        principal_signal = principal_signal or p_signal
        report["principal"][str(L)] = {
            "operator_checks": {
                "min_K_eigenvalue": float(
                    min(np.min(np.linalg.eigvalsh(sym(x))) for x in K)
                ),
                "min_c2": float(np.min(pvals)),
                "max_R_abs": float(np.max(np.abs(R))),
                "max_S_sym": float(np.max(np.abs(S + np.swapaxes(S, 1, 2)))),
            },
            "branches": pbranches,
            "mode_character_signal": bool(p_signal),
        }

        runs = {}
        for n in RESOLUTIONS:
            rr, w2, y, Kn = assemble_fem(r, K, G, S, M, n, grid_kind="r", frozen=False)
            diag = mode_diagnostics(rr, w2, y, Kn)
            runs[f"r_{n}"] = {
                "omega2": [float(x) for x in w2],
                "omega": [float(math.sqrt(x)) for x in w2 if x > 0],
                "diagnostics": diag,
            }
        # Coordinate/grid robustness probe at the middle resolution.
        rr, w2, y, Kn = assemble_fem(
            r, K, G, S, M, RESOLUTIONS[1], grid_kind="u", frozen=False
        )
        runs[f"u_{RESOLUTIONS[1]}"] = {
            "omega2": [float(x) for x in w2],
            "omega": [float(math.sqrt(x)) for x in w2 if x > 0],
            "diagnostics": mode_diagnostics(rr, w2, y, Kn),
        }
        # Frozen-coefficient null control.
        rr, w2, y, Kn = assemble_fem(
            r, K, G, S, M, RESOLUTIONS[1], grid_kind="r", frozen=True
        )
        frozen_diag = mode_diagnostics(rr, w2, y, Kn)
        runs[f"frozen_{RESOLUTIONS[1]}"] = {
            "omega2": [float(x) for x in w2],
            "omega": [float(math.sqrt(x)) for x in w2 if x > 0],
            "diagnostics": frozen_diag,
        }

        # Frequency convergence on the first positive modes.
        om160 = summarize_run(runs[f"r_{RESOLUTIONS[1]}"])
        om220 = summarize_run(runs[f"r_{RESOLUTIONS[2]}"])
        kcmp = min(len(om160), len(om220))
        freq_rel = [
            abs(om160[i] - om220[i]) / max(1.0, abs(om220[i])) for i in range(kcmp)
        ]
        omu = summarize_run(runs[f"u_{RESOLUTIONS[1]}"])
        kg = min(len(om160), len(omu))
        grid_rel = [
            abs(om160[i] - omu[i]) / max(1.0, abs(omu[i])) for i in range(kg)
        ]

        ss = runs[f"r_{RESOLUTIONS[2]}"]["diagnostics"]["rank_reordering"]
        fc = frozen_diag["rank_reordering"]
        source_diffs = {}
        for s in ("0", "1", "2"):
            source_diffs[s] = {
                "ssz_changes": bool(ss.get(s, {}).get("top_mode_changes", False)),
                "frozen_changes": bool(fc.get(s, {}).get("top_mode_changes", False)),
                "ssz_inner_top": ss.get(s, {}).get("inner_top_mode"),
                "ssz_outer_top": ss.get(s, {}).get("outer_top_mode"),
                "frozen_inner_top": fc.get(s, {}).get("inner_top_mode"),
                "frozen_outer_top": fc.get(s, {}).get("outer_top_mode"),
            }
        l_cavity_signal = any(
            v["ssz_changes"] and not v["frozen_changes"] for v in source_diffs.values()
        )
        cavity_signal = cavity_signal or l_cavity_signal

        report["finite_window"][str(L)] = {
            "runs": runs,
            "convergence": {
                "max_rel_omega_r160_vs_r220": float(max(freq_rel) if freq_rel else math.inf),
                "max_rel_omega_uniform_r_vs_uniform_u": float(max(grid_rel) if grid_rel else math.inf),
            },
            "source_rank_comparison": source_diffs,
            "local_cavity_selection_signal": bool(l_cavity_signal),
        }

    report["summary"] = {
        "principal_mode_character_signal_any_L": bool(principal_signal),
        "finite_window_weight_reordering_beyond_frozen_control_any_L": bool(cavity_signal),
        "local_proxy_verdict": (
            "LOCAL_KRGSM_SELECTION_SIGNAL"
            if principal_signal and cavity_signal
            else "LOCAL_KRGSM_NULL_OR_AMBIGUOUS"
        ),
        "publication_guard": (
            "This does not upgrade the preregistered global verdict. "
            "A center-to-infinity same-operator K,R,G,S,M export with physical "
            "boundary conditions remains required for SPECTRAL_SELECTION_PASS/NULL."
        ),
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")

    lines = [
        "# SSZ local KRGSM spectral-selection proxy",
        "",
        f"- member hash: `{report['member_hash']}`",
        f"- domain: {U_LO} < u < {U_HI}",
        "- boundary conditions: finite-window Dirichlet (diagnostic only)",
        f"- local proxy verdict: **{report['summary']['local_proxy_verdict']}**",
        "",
        "## Principal-symbol tracking",
        "",
    ]
    for L in L_VALUES:
        p = report["principal"][str(L)]
        lines.append(f"### L={L}")
        lines.append(
            f"min eig(K)={p['operator_checks']['min_K_eigenvalue']:.6e}, "
            f"min c^2={p['operator_checks']['min_c2']:.6e}, "
            f"mode-character signal={p['mode_character_signal']}"
        )
        for b in p["branches"]:
            lines.append(
                f"- branch {b['branch']}: c²=[{b['c2_min']:.6g},{b['c2_max']:.6g}], "
                f"max rotation={b['max_rotation_from_outer_deg']:.3f}°, "
                f"component switches={b['dominant_component_switches']}"
            )
        lines.append("")

    lines += ["## Finite-window full K,G,S,M eigenproblem", ""]
    for L in L_VALUES:
        f = report["finite_window"][str(L)]
        lines.append(
            f"### L={L}: local selection={f['local_cavity_selection_signal']}, "
            f"resolution Δω={f['convergence']['max_rel_omega_r160_vs_r220']:.3e}, "
            f"r/u-grid Δω={f['convergence']['max_rel_omega_uniform_r_vs_uniform_u']:.3e}"
        )
        for s, v in f["source_rank_comparison"].items():
            lines.append(
                f"- source {s}: SSZ {v['ssz_inner_top']}→{v['ssz_outer_top']} "
                f"(change={v['ssz_changes']}), frozen "
                f"{v['frozen_inner_top']}→{v['frozen_outer_top']} "
                f"(change={v['frozen_changes']})"
            )
        lines.append("")

    lines += [
        "## Guard",
        "",
        report["summary"]["publication_guard"],
        "",
    ]
    OUT_MD.write_text("\n".join(lines))
    print(json.dumps(report["summary"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())