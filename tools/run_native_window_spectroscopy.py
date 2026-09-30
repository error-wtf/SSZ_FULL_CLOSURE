#!/usr/bin/env python3
"""Actual finite-window SSZ spectroscopy, with global-QNM semantics fail-closed."""

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


def _relerr(a, b):
    a, b = np.asarray(a), np.asarray(b)
    m = min(len(a), len(b))
    if m == 0:
        return float("inf")
    return float(np.max(np.abs(a[:m] - b[:m]) / np.maximum(1.0, np.abs(b[:m]))))


def _solve_frame(d, red, L, stride=1):
    a = red.canonical_audit(d, int(L))
    args = (
        d.x.to_numpy(float),
        d.u.to_numpy(float),
        *(np.asarray(a[k], float) for k in ("K", "G", "S", "M")),
    )
    if stride != 1:
        args = downsample_native(*args, stride)
    return solve_near_zero_box_spectrum(*args, modes=NMODES)


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    build = build_onshell_central(ROOT)
    all41 = build.direct41.sort_values("x").reset_index(drop=True)
    red = module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")

    perL = {}
    any_binned_reordering = False
    boundary_sensitive = False
    all_resolution_ok = True

    for L in DEFAULT_L:
        lo, hi = WINDOWS["registered"]
        d = all41[(all41.u > lo) & (all41.u < hi)].sort_values("x").reset_index(drop=True)
        coarse = _solve_frame(d, red, L, stride=2)
        native = _solve_frame(d, red, L, stride=1)

        pos_c = coarse.omega2[coarse.omega2 > 0][:COMPARE]
        pos_n = native.omega2[native.omega2 > 0][:COMPARE]
        conv = _relerr(pos_c, pos_n)
        resolution_ok = bool(conv < 0.15)
        all_resolution_ok &= resolution_ok

        good = native.omega2 > 0
        W = native.weights[good][:COMPARE]
        point_pairs = robust_pairwise_inversions(W)
        centers, Zbin, Pbin, top2, neff = binned_residue_summary(native.r, W, bins=10)
        bin_pairs = robust_pairwise_inversions(Zbin, rel_margin=0.03)
        any_binned_reordering |= bool(bin_pairs)

        crops = {}
        sign_pattern = []
        for name, (wlo, whi) in WINDOWS.items():
            dd = all41[(all41.u > wlo) & (all41.u < whi)].sort_values("x").reset_index(drop=True)
            sp = _solve_frame(dd, red, L, stride=1)
            near = sp.omega2[: min(5, len(sp.omega2))]
            crops[name] = {
                "rows": int(len(dd)),
                "omega2_near_zero": [float(x) for x in near],
                "negative_near_zero_count": int(np.sum(sp.omega2 <= 0)),
            }
            sign_pattern.append(bool(np.any(sp.omega2 <= 0)))
        bc_sensitive = len(set(sign_pattern)) > 1
        boundary_sensitive |= bc_sensitive

        lam_grid, rho = spectral_density_lambda(native.omega2[good][:COMPARE], W)
        probe_idx = np.array([
            int(np.argmin(np.abs(native.u - q))) for q in (0.625, 0.65, 0.675, 0.695)
        ])
        np.savez_compressed(
            OUTDIR / f"L{int(L):04d}_finite_window_spectroscopy.npz",
            r=native.r,
            u=native.u,
            omega2=native.omega2,
            omega=native.omega,
            weights=native.weights,
            ipr=native.ipr,
            binned_r=centers,
            binned_residues=Zbin,
            binned_probabilities=Pbin,
            top2_share=top2,
            effective_mode_number=neff,
            lambda_grid=lam_grid,
            probe_u=native.u[probe_idx],
            rho_probe=rho[probe_idx],
        )

        perL[str(L)] = {
            "native_omega2_near_zero": [float(x) for x in native.omega2[:COMPARE]],
            "native_negative_near_zero_count": int(np.sum(native.omega2 <= 0)),
            "stride2_to_native_positive_omega2_max_relative_error": conv,
            "resolution_converged_15pct": resolution_ok,
            "pointwise_residue_inversion_pairs": [list(x) for x in point_pairs],
            "binned_residue_inversion_pairs": [list(x) for x in bin_pairs],
            "top2_spectral_share_range": [float(np.min(top2)), float(np.max(top2))],
            "effective_mode_number_range": [float(np.min(neff)), float(np.max(neff))],
            "ipr_first_positive_modes": [float(x) for x in native.ipr[good][:COMPARE]],
            "boundary_crop_checks": crops,
            "near_zero_sign_is_boundary_sensitive": bc_sensitive,
            "A_symmetry_scaled": native.symmetry_error_A,
            "B_symmetry_scaled": native.symmetry_error_B,
        }

    global_gate = {
        "status": "BLOCKED",
        "reason": "MISSING_DIRECT_GLOBAL_KRGM_CERTIFICATE",
    }
    if CERT.is_file():
        try:
            cert = require_direct_krgm_certificate(CERT, ROOT)
            global_gate = {
                "status": "READY",
                "certificate_release": cert.get("release"),
            }
        except RuntimeError as exc:
            global_gate = {"status": "BLOCKED", "reason": str(exc)}

    if not all_resolution_ok:
        status = "FINITE_WINDOW_SPECTROSCOPY_PARTIALLY_CONVERGED"
    elif boundary_sensitive:
        status = "FINITE_WINDOW_SPECTROSCOPY_BOUNDARY_SENSITIVE"
    elif any_binned_reordering:
        status = "FINITE_WINDOW_BINNED_WEIGHT_REORDERING"
    else:
        status = "FINITE_WINDOW_SPECTRAL_SELECTION_NULL"

    payload = {
        "scope": "registered finite-window spectroscopy diagnostic; not global QNM",
        "member_hash": getattr(build, "member_hash", None),
        "registered_window": list(WINDOWS["registered"]),
        "finite_boundary_condition": "homogeneous Dirichlet",
        "operator": "canonical K,G,S,M with antisymmetric S retained in FEM action",
        "mode_count": NMODES,
        "per_L": perL,
        "finite_window_status": status,
        "global_qnm_gate": global_gate,
        "weiss_li_semantics": (
            "Local/binned K-residue structure is a spectroscopy diagnostic. "
            "It is not yet a Weisz-style physical observable residue nor a Li-style "
            "mobility edge until the center-to-infinity retarded problem is certified."
        ),
        "interpretation_guard": (
            "Do not publish these finite-box frequencies as SSZ QNMs. Physical QNM/"
            "retarded-Green spectroscopy remains gated by DIRECT_GLOBAL_KRGM_EXPORT."
        ),
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    print(json.dumps(payload, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
