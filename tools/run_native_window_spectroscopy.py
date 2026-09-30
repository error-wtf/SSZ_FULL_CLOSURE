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


def _relerr(a, b):
    a, b = np.asarray(a), np.asarray(b)
    m = min(len(a), len(b))
    if m == 0:
        return float("inf")
    return float(np.max(np.abs(a[:m] - b[:m]) / np.maximum(1.0, np.abs(b[:m]))))


def _solve_arrays(r, u, K, G, S, M, mask=None, stride=1, modes=NMODES):
    if mask is not None:
        r, u, K, G, S, M = (
            np.asarray(A)[mask] for A in (r, u, K, G, S, M)
        )
    args = (r, u, K, G, S, M)
    if stride != 1:
        args = downsample_native(*args, stride)
    return solve_near_zero_box_spectrum(*args, modes=modes), args[2]

