#!/usr/bin/env python3
"""Finite-window SSZ spectroscopy with strict global-QNM fail-closed semantics.

This executes a real generalized radial eigenproblem using the canonical
K,G,S,M reduced operator on the registered healthy window.  It is *not* a
center-to-infinity QNM calculation.  Physical global spectroscopy remains
blocked unless DIRECT_GLOBAL_KRGM_CERTIFICATE exists and validates.

Diagnostics:
- finite-box eigenfrequencies/eigenfunctions for every required L;
- native vs stride-2 convergence by K-overlap mode matching;
- residue/radial-weight inversions;
- binned residue concentration and effective mode count;
- Lorentz-broadened local spectral density export;
- IR/UV boundary-window sensitivity;
- tachyon/negative-omega2 detection;
- explicit global QNM gate.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ssz_p5.config import DEFAULT_L
from ssz_p5.numerics import module
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central
from ssz_p5.qnm.gate import require_direct_krgm_certificate
from ssz_p5.qnm.native_window import (
    binned_residue_summary,
    downsample_native,
    match_modes_by_kinetic_overlap,
    robust_pairwise_inversions,
    solve_near_zero_box_spectrum,
    spectral_density_lambda,
)

OUTDIR = ROOT / "data/generated/spectral/native_window"
REPORT = ROOT / "data/generated/spectral/NATIVE_WINDOW_SPECTROSCOPY_REPORT.json"
CERT = ROOT / "data/certificates/DIRECT_GLOBAL_KRGM_CERTIFICATE.json"

NMODES = 14
COMPARE = 8
WINDOWS = {
    "registered": (0.62, 0.70),
    "inner_trim": (0.625, 0.70),
    "outer_trim": (0.62, 0.695),
    "both_trim": (0.625, 0.695),
}


def _solve_arrays(r, u, K, G, S, M, mask=None, stride=1, modes=NMODES):
    if mask is not None:
        r, u, K, G, S, M = (np.asarray(A)[mask] for A in (r, u, K, G, S, M))
    args = (r, u, K, G, S, M)
    if stride != 1:
        args = downsample_native(*args, stride)
    return solve_near_zero_box_spectrum(*args, modes=modes), args[2]


def _global_gate():
    if not CERT.is_file():
        return {
            "status": "BLOCKED",
            "blocker": "MISSING_DIRECT_GLOBAL_KRGM_CERTIFICATE",
            "detail": "Finite-window spectroscopy is allowed; physical global QNM/residue claims are not.",
        }
    try:
        require_direct_krgm_certificate(CERT, ROOT)
    except Exception as exc:
        return {
            "status": "BLOCKED",
            "blocker": "INVALID_DIRECT_GLOBAL_KRGM_CERTIFICATE",
            "detail": str(exc),
        }
    return {"status": "READY", "blocker": None}


def _window_mask(u, lo, hi):
    return (u > lo) & (u < hi)


def _pairset(pairs):
    return {tuple(map(int, p)) for p in pairs}


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    REPORT.parent.mkdir(parents=True, exist_ok=True)

    build = build_onshell_central(ROOT)
    stream = build.direct41.sort_values("x").reset_index(drop=True)
    r = stream.x.to_numpy(float)
    u = stream.u.to_numpy(float)
    reducer = module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")

    reg = _window_mask(u, *WINDOWS["registered"])
    if np.sum(reg) < 100:
        raise RuntimeError("registered production window is unexpectedly small")

    perL = {}
    all_bulk_converged = True
    any_robust_inversion = False
    any_boundary_sensitive = False

    for L in DEFAULT_L:
        op = reducer.canonical_audit(stream, int(L))
        K = np.asarray(op["K"], float)
        G = np.asarray(op["G"], float)
        S = np.asarray(op["S"], float)
        M = np.asarray(op["M"], float)

        # Health guard before spectroscopy.
        Ks = (K[reg] + np.swapaxes(K[reg], 1, 2)) / 2
        minK = float(np.min(np.linalg.eigvalsh(Ks)))
        if minK <= 0:
            raise RuntimeError(f"L={L}: registered window K is not positive: {minK}")

        native, Kreg = _solve_arrays(r, u, K, G, S, M, mask=reg, stride=1)
        coarse, _ = _solve_arrays(r, u, K, G, S, M, mask=reg, stride=2)

        matches = match_modes_by_kinetic_overlap(coarse, native, Kreg)
        compared = [m for m in matches if m["native_index"] < COMPARE and m["coarse_index"] < COMPARE]
        conv_err = max((m["relative_omega2_error"] for m in compared), default=float("inf"))
        min_overlap = min((m["overlap"] for m in compared), default=0.0)
        converged = bool(conv_err < 5e-2 and min_overlap > 0.70)
        all_bulk_converged &= converged

        native_pairs = _pairset(robust_pairwise_inversions(native.weights[:COMPARE], rel_margin=0.05))
        coarse_pairs = _pairset(robust_pairwise_inversions(coarse.weights[:COMPARE], rel_margin=0.05))
        common_pairs = sorted(native_pairs & coarse_pairs)
        any_robust_inversion |= bool(common_pairs)

        # Boundary sensitivity: compare first positive modes across trimmed boxes.
        boundary = {}
        base_pos = native.omega2[native.omega2 > 0][:COMPARE]
        for name, (lo, hi) in WINDOWS.items():
            mask = _window_mask(u, lo, hi)
            sp, _ = _solve_arrays(r, u, K, G, S, M, mask=mask, stride=2, modes=max(COMPARE + 4, 12))
            pos = sp.omega2[sp.omega2 > 0][:COMPARE]
            m = min(len(base_pos), len(pos))
            rel = (
                float(np.max(np.abs(pos[:m] - base_pos[:m]) / np.maximum(1.0, np.abs(base_pos[:m]))))
                if m else float("inf")
            )
            boundary[name] = {"positive_modes_compared": int(m), "max_relative_shift": rel}
        ir_sensitive = max(v["max_relative_shift"] for k, v in boundary.items() if k != "registered") > 0.10
        any_boundary_sensitive |= ir_sensitive

        centers, Z, P, top2, neff = binned_residue_summary(native.r, native.weights[:COMPARE], bins=10)
        lam_grid, rho = spectral_density_lambda(native.omega2, native.weights)

        np.savez_compressed(
            OUTDIR / f"L{int(L)}_native_window_spectroscopy.npz",
            r=native.r,
            u=native.u,
            omega2=native.omega2,
            omega=native.omega,
            psi=native.psi,
            weights=native.weights,
            ipr=native.ipr,
            spectral_lambda=lam_grid,
            spectral_density=rho,
            binned_centers=centers,
            binned_Z=Z,
            binned_P=P,
            binned_top2=top2,
            binned_neff=neff,
        )

        perL[str(int(L))] = {
            "min_K_registered": minK,
            "omega2": [float(x) for x in native.omega2],
            "positive_mode_count": int(np.sum(native.omega2 > 0)),
            "finite_box_tachyon_present": bool(np.any(native.omega2 < 0)),
            "omega2_resolution_max_relative_error": float(conv_err),
            "min_native_coarse_K_overlap": float(min_overlap),
            "bulk_resolution_converged": converged,
            "robust_inversion_pairs_native": [list(x) for x in sorted(native_pairs)],
            "robust_inversion_pairs_coarse": [list(x) for x in sorted(coarse_pairs)],
            "robust_inversion_pairs_common_to_both_resolutions": [list(x) for x in common_pairs],
            "boundary_sensitivity": boundary,
            "ir_boundary_sensitive": bool(ir_sensitive),
            "ipr_first_modes": [float(x) for x in native.ipr[:COMPARE]],
            "binned_top2_fraction_minmax": [float(np.min(top2)), float(np.max(top2))],
            "binned_effective_mode_count_minmax": [float(np.min(neff)), float(np.max(neff))],
            "pencil_symmetry_error_A": float(native.symmetry_error_A),
            "pencil_symmetry_error_B": float(native.symmetry_error_B),
        }

    if not all_bulk_converged:
        status = "FINITE_WINDOW_SPECTROSCOPY_BULK_NOT_CONVERGED"
    elif any_robust_inversion and any_boundary_sensitive:
        status = "FINITE_WINDOW_BULK_RESIDUE_REORDERING_WITH_IR_BOUNDARY_SENSITIVITY"
    elif any_boundary_sensitive:
        status = "FINITE_WINDOW_SPECTROSCOPY_IR_BOUNDARY_SENSITIVE"
    elif any_robust_inversion:
        status = "FINITE_WINDOW_BULK_RESIDUE_REORDERING"
    else:
        status = "FINITE_WINDOW_SPECTRAL_SELECTION_NULL"

    payload = {
        "finite_window_status": status,
        "finite_window_semantics": (
            "Real finite-box K,G,S,M generalized eigenproblem on the healthy registered window. "
            "Useful for falsification, convergence and residue-structure diagnostics; not a physical global QNM spectrum."
        ),
        "global_qnm_gate": _global_gate(),
        "windows": WINDOWS,
        "modes_requested": NMODES,
        "modes_compared": COMPARE,
        "all_bulk_resolution_converged": bool(all_bulk_converged),
        "any_robust_residue_reordering": bool(any_robust_inversion),
        "any_boundary_sensitivity": bool(any_boundary_sensitive),
        "per_L": perL,
        "interpretation_guard": (
            "A Weisz-style physical spectral-selection claim requires a certified center-to-infinity direct-global "
            "KRGM eigenproblem and boundary-condition-stable residues. Finite-box reorderings are candidate structure only."
        ),
    }
    REPORT.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    print(json.dumps(payload, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
