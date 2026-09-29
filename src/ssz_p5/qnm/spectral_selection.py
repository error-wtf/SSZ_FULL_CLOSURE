"""Fail-closed radial spectral-weight selection diagnostics.

This module intentionally does not manufacture a global eigenproblem.  It only
post-processes a hash-bound spectral export produced from the certified direct
global KRGM operator.

Expected mode export convention:
    r            shape (N,)
    omega2       shape (M,)
    psi          shape (M,N,F), real or complex
    K            shape (N,F,F)

The local spectral weight is the positive kinetic density
    Z_n(r) = psi_n(r)^\dagger K(r) psi_n(r)
after global K-normalization.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class SpectralSelectionResult:
    status: str
    reordered: bool
    rankings: np.ndarray
    weights: np.ndarray
    normalized_modes: np.ndarray


def _validate_inputs(r, K, psi):
    r = np.asarray(r, float)
    K = np.asarray(K)
    psi = np.asarray(psi)
    if r.ndim != 1 or len(r) < 2 or np.any(~np.isfinite(r)) or np.any(np.diff(r) <= 0):
        raise ValueError("r must be finite, unique, and strictly increasing")
    if K.ndim != 3 or K.shape[0] != len(r) or K.shape[1] != K.shape[2]:
        raise ValueError("K must have shape (N,F,F)")
    if psi.ndim != 3 or psi.shape[1] != len(r) or psi.shape[2] != K.shape[1]:
        raise ValueError("psi must have shape (M,N,F)")
    if not np.all(np.isfinite(K)):
        raise ValueError("K contains non-finite values")
    if not np.all(np.isfinite(psi)):
        raise ValueError("psi contains non-finite values")
    herm = K - np.swapaxes(K.conj(), 1, 2)
    scale = max(1.0, float(np.max(np.abs(K))))
    if float(np.max(np.abs(herm))) / scale > 1e-10:
        raise ValueError("K is not Hermitian")
    evals = np.linalg.eigvalsh((K + np.swapaxes(K.conj(), 1, 2)) / 2)
    if float(np.min(evals)) <= 0:
        raise ValueError("K is not positive definite")
    return r, K, psi


def local_kinetic_weights(r, K, psi):
    """Return globally K-normalized modes and local kinetic densities Z_n(r)."""
    r, K, psi = _validate_inputs(r, K, psi)
    density = np.einsum("mni,nij,mnj->mn", psi.conj(), K, psi).real
    if np.min(density) < -1e-10:
        raise ValueError("negative local kinetic density")
    density = np.maximum(density, 0.0)
    norms = np.trapezoid(density, r, axis=1)
    if np.any(~np.isfinite(norms)) or np.any(norms <= 0):
        raise ValueError("non-positive K norm")
    normalized = psi / np.sqrt(norms)[:, None, None]
    weights = density / norms[:, None]
    return normalized, weights


def ranking_by_radius(weights):
    """Return mode indices sorted by descending Z_n(r) at every radial node."""
    weights = np.asarray(weights, float)
    if weights.ndim != 2 or not np.all(np.isfinite(weights)):
        raise ValueError("weights must be a finite (M,N) array")
    return np.argsort(-weights, axis=0, kind="stable").T


def has_pairwise_reordering(weights, *, atol=1e-12, rtol=1e-9):
    """Detect a robust pairwise rank inversion anywhere in radius."""
    weights = np.asarray(weights, float)
    if weights.ndim != 2:
        raise ValueError("weights must be two-dimensional")
    modes = weights.shape[0]
    for a in range(modes):
        for b in range(a + 1, modes):
            delta = weights[a] - weights[b]
            tol = atol + rtol * np.maximum(np.abs(weights[a]), np.abs(weights[b]))
            if np.any(delta > tol) and np.any(delta < -tol):
                return True
    return False


def evaluate_spectral_selection(r, K, psi):
    """Evaluate the pre-registered PASS/NULL spectral-selection observable."""
    normalized, weights = local_kinetic_weights(r, K, psi)
    rankings = ranking_by_radius(weights)
    reordered = has_pairwise_reordering(weights)
    return SpectralSelectionResult(
        status="SPECTRAL_SELECTION_PASS" if reordered else "SPECTRAL_SELECTION_NULL",
        reordered=reordered,
        rankings=rankings,
        weights=weights,
        normalized_modes=normalized,
    )
