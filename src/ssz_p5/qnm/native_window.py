"""Finite-window spectroscopy diagnostics for healthy SSZ reduced operators.

This module does NOT implement physical center-to-infinity QNMs. It solves a
finite radial generalized eigenproblem derived directly from the canonical
quadratic form

    L = dot(Y)^T K dot(Y) - Y'^T G Y' + Y'^T S Y - Y^T M Y,

with K,G,M symmetric and S antisymmetric. The spatial form is discretized by
piecewise-linear finite elements on the native radial grid, preserving the
symmetry of the generalized pencil A psi = omega^2 B psi.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh


@dataclass(frozen=True)
class BoxSpectrum:
    r: np.ndarray
    u: np.ndarray
    omega2: np.ndarray
    omega: np.ndarray
    psi: np.ndarray
    weights: np.ndarray
    ipr: np.ndarray
    symmetry_error_A: float
    symmetry_error_B: float


def _trapezoid_weights(r):
    r = np.asarray(r, float)
    h = np.diff(r)
    if len(r) < 5 or np.any(h <= 0):
        raise ValueError("r must be strictly increasing with at least five nodes")
    w = np.empty(len(r), float)
    w[0] = h[0] / 2
    w[-1] = h[-1] / 2
    w[1:-1] = (h[:-1] + h[1:]) / 2
    return w


def assemble_dirichlet_fem(r, K, G, S, M):
    """Return sparse symmetric finite-element matrices and node quadrature."""
    r = np.asarray(r, float)
    K, G, S, M = (np.asarray(A, float) for A in (K, G, S, M))
    n, f = K.shape[:2]
    if any(A.shape != (n, f, f) for A in (G, S, M)):
        raise ValueError("operator matrices must share shape (N,F,F)")
    w = _trapezoid_weights(r)

    rows, cols, avals, bvals = [], [], [], []

    def add_block(i, j, A, B=None):
        bi, bj = i * f, j * f
        for a in range(f):
            for b in range(f):
                rows.append(bi + a)
                cols.append(bj + b)
                avals.append(float(A[a, b]))
                bvals.append(0.0 if B is None else float(B[a, b]))

    for i in range(n):
        Ks = (K[i] + K[i].T) / 2
        Ms = (M[i] + M[i].T) / 2
        add_block(i, i, w[i] * Ms, w[i] * Ks)

    h = np.diff(r)
    for i in range(n - 1):
        Gi = (G[i] + G[i + 1]) / 2
        Gi = (Gi + Gi.T) / 2
        Si = (S[i] + S[i + 1]) / 2
        Si = (Si - Si.T) / 2
        hh = h[i]
        add_block(i, i, Gi / hh)
        add_block(i + 1, i + 1, Gi / hh)
        add_block(i, i + 1, -Gi / hh + 0.5 * Si)
        add_block(i + 1, i, -Gi / hh + 0.5 * Si.T)

    nd = n * f
    A = coo_matrix((avals, (rows, cols)), shape=(nd, nd)).tocsr()
    B = coo_matrix((bvals, (rows, cols)), shape=(nd, nd)).tocsr()

    keep = np.arange(f, (n - 1) * f)
    A = A[keep][:, keep]
    B = B[keep][:, keep]
    return A, B, w, keep


def solve_near_zero_box_spectrum(r, u, K, G, S, M, *, modes=14):
    """Solve the modes closest to omega^2=0 using sparse shift-invert."""
    r = np.asarray(r, float)
    u = np.asarray(u, float)
    A, B, w, keep = assemble_dirichlet_fem(r, K, G, S, M)

    da = A - A.T
    db = B - B.T
    asym_a = float(np.max(np.abs(da.data))) if da.nnz else 0.0
    asym_b = float(np.max(np.abs(db.data))) if db.nnz else 0.0
    scale_a = max(1.0, float(np.max(np.abs(A.data)))) if A.nnz else 1.0
    scale_b = max(1.0, float(np.max(np.abs(B.data)))) if B.nnz else 1.0
    asym_a /= scale_a
    asym_b /= scale_b
    if asym_a > 1e-10 or asym_b > 1e-10:
        raise ValueError("assembled generalized pencil is not symmetric")

    k = min(int(modes), A.shape[0] - 2)
    vals, vecs = eigsh(
        A,
        M=B,
        k=k,
        sigma=0.0,
        which="LM",
        tol=1e-10,
        maxiter=10000,
    )
    order = np.argsort(vals)
    vals = np.asarray(vals[order], float)
    vecs = np.asarray(vecs[:, order], float)

    n = len(r)
    f = K.shape[1]
    psi = np.zeros((k, n, f), float)
    for q in range(k):
        full = np.zeros(n * f, float)
        full[keep] = vecs[:, q]
        psi[q] = full.reshape(n, f)

    density = np.einsum("mni,nij,mnj->mn", psi, K, psi)
    norms = np.sum(density * w[None, :], axis=1)
    if np.any(~np.isfinite(norms)) or np.any(norms <= 0):
        raise ValueError("non-positive modal K norm")
    psi /= np.sqrt(norms)[:, None, None]
    weights = np.einsum("mni,nij,mnj->mn", psi, K, psi)

    span = r[-1] - r[0]
    ipr = np.trapezoid((weights * span) ** 2, r, axis=1) / span
    omega = np.where(vals > 0, np.sqrt(vals), np.nan)
    return BoxSpectrum(
        r=r,
        u=u,
        omega2=vals,
        omega=omega,
        psi=psi,
        weights=weights,
        ipr=ipr,
        symmetry_error_A=asym_a,
        symmetry_error_B=asym_b,
    )


def downsample_native(r, u, K, G, S, M, stride):
    """Stride the native grid without interpolation; always retain the last row."""
    idx = np.arange(0, len(r), int(stride))
    idx = np.unique(np.r_[idx, len(r) - 1])
    return (
        np.asarray(r)[idx],
        np.asarray(u)[idx],
        np.asarray(K)[idx],
        np.asarray(G)[idx],
        np.asarray(S)[idx],
        np.asarray(M)[idx],
    )


def robust_pairwise_inversions(weights, *, rel_margin=0.05):
    """Mode pairs with a robust reversal in normalized local residue."""
    W = np.asarray(weights, float)
    total = np.sum(W, axis=0)
    P = W / np.maximum(total[None, :], 1e-300)
    pairs = []
    for a in range(len(P)):
        for b in range(a + 1, len(P)):
            d = P[a] - P[b]
            if float(np.max(d)) > rel_margin and float(np.min(d)) < -rel_margin:
                pairs.append((a, b))
    return pairs


def binned_residue_summary(r, weights, *, bins=10):
    """Integrate local residues in radial bins to suppress nodal point artifacts."""
    r = np.asarray(r, float)
    W = np.asarray(weights, float)
    edges = np.linspace(r[0], r[-1], int(bins) + 1)
    Z = np.zeros((len(W), bins), float)
    centers = np.empty(bins, float)
    for j in range(bins):
        m = (r >= edges[j]) & (r <= edges[j + 1])
        centers[j] = 0.5 * (edges[j] + edges[j + 1])
        if np.sum(m) >= 2:
            Z[:, j] = np.trapezoid(W[:, m], r[m], axis=1)
    P = Z / np.maximum(np.sum(Z, axis=0, keepdims=True), 1e-300)
    top2 = np.sort(P, axis=0)[-2:].sum(axis=0)
    neff = 1.0 / np.maximum(np.sum(P**2, axis=0), 1e-300)
    return centers, Z, P, top2, neff


def spectral_density_lambda(omega2, weights, *, n_lambda=600, eta_fraction=2e-3):
    """Lorentz-broadened modal spectral density in lambda=omega^2."""
    lam = np.asarray(omega2, float)
    good = lam > 0
    lam = lam[good]
    W = np.asarray(weights, float)[good]
    if len(lam) < 2:
        raise ValueError("too few positive modes")
    lo = max(0.0, float(lam.min()) * 0.95)
    hi = float(lam.max()) * 1.05
    grid = np.linspace(lo, hi, int(n_lambda))
    eta = max((hi - lo) * float(eta_fraction), 1e-12)
    rho = np.sum(
        W[:, :, None]
        * (eta / np.pi)
        / ((grid[None, None, :] - lam[:, None, None]) ** 2 + eta**2),
        axis=0,
    )
    return grid, rho
