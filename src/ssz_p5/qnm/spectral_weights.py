"""Fail-closed modal-weight diagnostics.

For self-adjoint bound/normal-mode problems with positive K, the local kinetic
weight is

    Z_n(r) = psi_n(r)^dagger K(r) psi_n(r)

after global K-normalization.

Open QNM/resonance problems require a Green-function/left-right residue
definition and must not be silently mapped onto this bound-mode formula.
"""

from __future__ import annotations

import numpy as np


def normalized_kinetic_weights(r, K, psi):
    r = np.asarray(r, float)
    K = np.asarray(K)
    psi = np.asarray(psi)
    if r.ndim != 1 or len(r) < 2 or np.any(np.diff(r) <= 0):
        raise ValueError("r must be strictly increasing")
    if K.shape != (len(r), psi.shape[-1], psi.shape[-1]):
        raise ValueError("K shape mismatch")
    if psi.ndim != 3 or psi.shape[1] != len(r):
        raise ValueError("psi must have shape (modes, radii, fields)")

    ks = (K + np.swapaxes(K.conj(), 1, 2)) / 2
    if np.max(np.abs(K - np.swapaxes(K.conj(), 1, 2))) > 1e-9 * max(
        1.0, float(np.max(np.abs(K)))
    ):
        raise ValueError("K is not Hermitian")
    if np.min(np.linalg.eigvalsh(ks)) <= 0:
        raise ValueError("K is not positive definite")

    density = np.einsum("mni,nij,mnj->mn", psi.conj(), ks, psi).real
    if np.min(density) < -1e-10:
        raise ValueError("negative kinetic density")
    density = np.maximum(density, 0.0)
    norm = np.trapezoid(density, r, axis=1)
    if np.any(norm <= 0):
        raise ValueError("nonpositive mode norm")
    weights = density / norm[:, None]
    return weights


def ranking(weights):
    w = np.asarray(weights, float)
    if w.ndim != 2 or not np.all(np.isfinite(w)):
        raise ValueError("weights must be finite 2D")
    return np.argsort(-w, axis=0, kind="stable").T


def has_pairwise_reordering(weights, atol=1e-12, rtol=1e-9):
    w = np.asarray(weights, float)
    for a in range(w.shape[0]):
        for b in range(a + 1, w.shape[0]):
            d = w[a] - w[b]
            tol = atol + rtol * np.maximum(np.abs(w[a]), np.abs(w[b]))
            if np.any(d > tol) and np.any(d < -tol):
                return True
    return False
