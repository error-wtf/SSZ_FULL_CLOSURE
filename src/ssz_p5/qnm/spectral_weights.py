"""Spectral-weight observables for bound modes and non-Hermitian resonances.

The two cases are deliberately separated:

* A positive, self-adjoint bound-mode problem may use the local kinetic density
  psi^dagger K psi after global K normalization.
* A QNM/resonance problem is non-Hermitian.  Its observable pole weight must be
  formed from left/right modes and the frequency derivative of the operator.

This module contains no SSZ-specific fitting parameters.
"""

from __future__ import annotations

import numpy as np


def kinetic_density_weights(r, k_metric, modes):
    """Globally K-normalize bound modes and return local positive densities.

    Parameters
    ----------
    r:
        Strictly increasing radial grid, shape (N,).
    k_metric:
        Hermitian positive matrix K(r), shape (N,F,F).
    modes:
        Mode profiles, shape (M,N,F).
    """
    r = np.asarray(r, float)
    k_metric = np.asarray(k_metric)
    modes = np.asarray(modes)

    if r.ndim != 1 or len(r) < 2 or np.any(~np.isfinite(r)) or np.any(np.diff(r) <= 0):
        raise ValueError("r must be finite and strictly increasing")
    if k_metric.ndim != 3 or k_metric.shape[0] != len(r):
        raise ValueError("K must have shape (N,F,F)")
    if k_metric.shape[1] != k_metric.shape[2]:
        raise ValueError("K must be square in field space")
    if modes.ndim != 3 or modes.shape[1:] != (len(r), k_metric.shape[1]):
        raise ValueError("modes must have shape (M,N,F)")

    kh = (k_metric + np.swapaxes(k_metric.conj(), 1, 2)) / 2
    asym = np.max(np.abs(k_metric - np.swapaxes(k_metric.conj(), 1, 2)))
    scale = max(1.0, float(np.max(np.abs(k_metric))))
    if asym / scale > 1e-10:
        raise ValueError("K is not Hermitian")
    if np.min(np.linalg.eigvalsh(kh)) <= 0:
        raise ValueError("K is not positive definite")

    density = np.einsum("mni,nij,mnj->mn", modes.conj(), kh, modes).real
    if np.min(density) < -1e-10:
        raise ValueError("negative kinetic density")
    density = np.maximum(density, 0.0)
    norms = np.trapezoid(density, r, axis=1)
    if np.any(~np.isfinite(norms)) or np.any(norms <= 0):
        raise ValueError("invalid K norm")
    return modes / np.sqrt(norms)[:, None, None], density / norms[:, None]


def qnm_pole_weight(left, right, d_operator_domega, probe):
    """Return the absolute observable residue factor for one simple QNM pole.

    For a nonlinear frequency-domain operator L(omega), a simple pole has
    residue proportional to

        (p^dagger r) (l^dagger p) /
        (l^dagger [dL/domega] r).

    This definition is invariant under reciprocal rescaling of left/right
    eigenvectors and does not pretend that QNMs have a positive Hilbert norm.
    """
    left = np.asarray(left, complex)
    right = np.asarray(right, complex)
    deriv = np.asarray(d_operator_domega, complex)
    probe = np.asarray(probe, complex)

    n = right.size
    if left.shape != (n,) or probe.shape != (n,) or deriv.shape != (n, n):
        raise ValueError("incompatible residue shapes")
    if not (
        np.all(np.isfinite(left))
        and np.all(np.isfinite(right))
        and np.all(np.isfinite(deriv))
        and np.all(np.isfinite(probe))
    ):
        raise ValueError("non-finite residue input")

    denominator = np.vdot(left, deriv @ right)
    if abs(denominator) <= 1e-14:
        raise ValueError("singular QNM residue normalization")
    numerator = np.vdot(probe, right) * np.vdot(left, probe)
    return float(abs(numerator / denominator))


def ranking_by_radius(weights):
    """Mode indices ordered by descending real non-negative local weight."""
    weights = np.asarray(weights, float)
    if weights.ndim != 2 or np.any(~np.isfinite(weights)) or np.any(weights < 0):
        raise ValueError("weights must be finite and non-negative with shape (M,N)")
    return np.argsort(-weights, axis=0, kind="stable").T


def has_pairwise_reordering(weights, *, atol=1e-12, rtol=1e-9):
    """Detect a robust pairwise ranking inversion somewhere on the radial grid."""
    weights = np.asarray(weights, float)
    if weights.ndim != 2 or np.any(~np.isfinite(weights)):
        raise ValueError("weights must be a finite 2D array")
    for a in range(weights.shape[0]):
        for b in range(a + 1, weights.shape[0]):
            delta = weights[a] - weights[b]
            tol = atol + rtol * np.maximum(abs(weights[a]), abs(weights[b]))
            if np.any(delta > tol) and np.any(delta < -tol):
                return True
    return False
