#!/usr/bin/env python3
"""Run finite-window SSZ spectroscopy now, while global QNM stays fail-closed.

This is an actual eigenmode/residue calculation on the healthy registered
production window of the current locked Electric member.  It is deliberately
reported as FINITE_WINDOW_DIAGNOSTIC, not as a global QNM spectrum.

For every default L it:
  * regenerates canonical K,G,S,M from the locked current member;
  * restricts to 0.62<u<0.70 where K and radial principal gates pass;
  * solves the finite-box generalized eigenproblem at two resolutions;
  * computes local K-normalized residues Z_n(r);
  * detects robust pairwise residue-rank inversions;
  * checks eigenvalue convergence;
  * exports representative local spectral densities rho(r,lambda).

The global direct-KRGM/QNM gate is independently checked and remains BLOCKED
unless its real certificate exists and validates.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ssz_p5.config import DEFAULT_L  # noqa:E402
from ssz_p5.numerics import module  # noqa:E402
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central  # noqa:E402
from ssz_p5.qnm.native_window import (  # noqa:E402
    resample_operator,
    robust_pairwise_inversions,
    solve_box_spectrum,
    spectral_density_lambda,
)
from ssz_p5.qnm.gate import require_direct_krgm_certificate  # noqa:E402

OUTDIR = ROOT / "data/generated/spectral/native_window"
REPORT = ROOT / "data/generated/spectral/NATIVE_WINDOW_SPECTROSCOPY_REPORT.json"
CERT = ROOT / "data/certificates/DIRECT_GLOBAL_KRGM_CERTIFICATE.json"

RESOLUTIONS = (96, 128)
NMODES = 14
COMPARE_MODES = 8


def relerr(a, b):
    return np.abs(np.asarray(a) - np.asarray(b)) / np.maximum(1.0, np.abs(b))


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    build = build_onshell_central(ROOT)
    d = build.direct41.sort_values("x").reset_index(drop=True)
    mask = (d.u.to_numpy(float) > 0.62) & (d.u.to_numpy(float) < 0.70)
    d = d.loc[mask].sort_values("x").reset_index(drop=True)
    red = module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")

    perL = {}
    any_reordering = False
    all_converged = True
    any_tachyon = False

    for L in DEFAULT_L:
        a = red.canonical_audit(d, int(L))
        base = dict(
            r=d.x.to_numpy(float),
            u=d.u.to_numpy(float),
            K=np.asarray(a["K"], float),
            G=np.asarray(a["G"], float),
            S=np.asarray(a["S"], float),
            M=np.asarray(a["M"], float),
        )
        runs = {}
        spectra = {}
        for n in RESOLUTIONS:
            rg, ug, K, G, S, M = resample_operator(**base, n=n)
            sp = solve_box_spectrum(rg, ug, K, G, S, M, modes=NMODES)
            spectra[n] = sp
            positive = sp.omega2 > 0
            inv = robust_pairwise_inversions(sp.weights[positive][:COMPARE_MODES])
            runs[str(n)] = {
                "negative_omega2_count": sp.negative_omega2_count,
                "omega2_first": [float(x) for x in sp.omega2[:COMPARE_MODES]],
                "omega_first_positive": [
                    float(x) for x in sp.omega[positive][:COMPARE_MODES]
                ],
                "ipr_first_positive": [
                    float(x) for x in sp.ipr[positive][:COMPARE_MODES]
                ],
                "robust_inversion_pairs": [list(x) for x in inv],
                "A_symmetry_scaled": sp.symmetry_error_A,
                "B_symmetry_scaled": sp.symmetry_error_B,
            }

        lo, hi = spectra[RESOLUTIONS[0]], spectra[RESOLUTIONS[1]]
        pos_lo = lo.omega2[lo.omega2 > 0][:COMPARE_MODES]
        pos_hi = hi.omega2[hi.omega2 > 0][:COMPARE_MODES]
        m = min(len(pos_lo), len(pos_hi))
        conv = float(np.max(relerr(pos_lo[:m], pos_hi[:m]))) if m else float("inf")
        pairs_lo = {tuple(x) for x in robust_pairwise_inversions(lo.weights[lo.omega2 > 0][:COMPARE_MODES])}
        pairs_hi = {tuple(x) for x in robust_pairwise_inversions(hi.weights[hi.omega2 > 0][:COMPARE_MODES])}
        robust_pairs = sorted(pairs_lo & pairs_hi)
        converged = bool(m >= 4 and conv < 0.05)
        tachyon = bool(lo.negative_omega2_count or hi.negative_omega2_count)

        any_reordering |= bool(robust_pairs)
        all_converged &= converged
        any_tachyon |= tachyon

        # Export high-resolution modal data and a few local spectral-density slices.
        good = hi.omega2 > 0
        lam_grid, rho = spectral_density_lambda(
            hi.omega2[good][:COMPARE_MODES],
            hi.weights[good][:COMPARE_MODES],
        )
        probe_idx = np.array([
            int(np.argmin(np.abs(hi.u - q))) for q in (0.625, 0.65, 0.675, 0.695)
        ])
        np.savez_compressed(
            OUTDIR / f"L{int(L):04d}_finite_window_spectroscopy.npz",
            r=hi.r,
            u=hi.u,
            omega2=hi.omega2,
            omega=hi.omega,
            weights=hi.weights,
            ipr=hi.ipr,
            lambda_grid=lam_grid,
            probe_u=hi.u[probe_idx],
            rho_probe=rho[probe_idx],
        )

        perL[str(L)] = {
            "resolutions": runs,
            "positive_mode_match_count": int(m),
            "omega2_resolution_max_relative_error": conv,
            "resolution_converged_5pct": converged,
            "robust_inversion_pairs_common_to_both_resolutions": [list(x) for x in robust_pairs],
            "finite_box_tachyon_present": tachyon,
        }

    global_gate = {"status": "BLOCKED", "reason": "MISSING_DIRECT_GLOBAL_KRGM_CERTIFICATE"}
    if CERT.is_file():
        try:
            cert = require_direct_krgm_certificate(CERT, ROOT)
            global_gate = {"status": "READY", "certificate_release": cert.get("release")}
        except RuntimeError as exc:
            global_gate = {"status": "BLOCKED", "reason": str(exc)}

    if any_tachyon:
        finite_status = "FINITE_WINDOW_SPECTROSCOPY_TACHYONIC"
    elif all_converged and any_reordering:
        finite_status = "FINITE_WINDOW_SPECTRAL_WEIGHT_REORDERING"
    elif all_converged:
        finite_status = "FINITE_WINDOW_SPECTRAL_SELECTION_NULL"
    else:
        finite_status = "FINITE_WINDOW_SPECTROSCOPY_NOT_CONVERGED"

    payload = {
        "scope": "finite healthy production-window spectroscopy diagnostic; not global QNM",
        "member_hash": getattr(build, "member_hash", None),
        "u_window": [0.62, 0.70],
        "boundary_condition": "homogeneous Dirichlet at finite-box endpoints",
        "canonical_action_operator": "K,G,S,M including antisymmetric S via quadratic-form FEM",
        "resolutions": list(RESOLUTIONS),
        "modes_requested": NMODES,
        "per_L": perL,
        "finite_window_status": finite_status,
        "global_qnm_gate": global_gate,
        "interpretation_guard": (
            "A finite-window residue reordering is a diagnostic of local modal-weight structure. "
            "It is not a physical center-to-infinity QNM claim and cannot replace the direct-global "
            "KRGM certificate, outgoing boundary condition, or global retarded Green function."
        ),
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    print(json.dumps(payload, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
