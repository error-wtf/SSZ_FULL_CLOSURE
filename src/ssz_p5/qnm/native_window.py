"""Finite-window spectroscopy diagnostics for a healthy reduced SSZ operator.

This module is intentionally separate from physical global QNM spectroscopy.
It discretizes the canonical quadratic action on a finite radial interval with
Dirichlet endpoints.  It is useful for falsification, convergence and local
spectral-weight diagnostics while the direct-global KRGM gate remains closed.

For the accepted canonical action,
    L = dot(Y)^T K dot(Y) - Y'^T G Y' + Y'^T S Y - Y^T M Y,
with K,G,M symmetric and S antisymmetric, the spatial quadratic form is
    Q = integral [Y'^T G Y' - Y'^T S Y + Y^T M Y] dr.
The generalized finite-box problem is A psi = omega^2 B psi.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import eigh
from scipy.interpolate import PchipInterpolator


@dataclass(frozen=True)
class BoxSpectrum:
    r: np.ndarray
    u: np.ndarray
    omega2: np.ndarray
    omega: np.ndarray
    psi: np.ndarray
    weights: np.ndarray
    ipr: np.ndarray
    negative_omega2_count: int
    symmetry_error_A: float
    symmetry_error_B: float


def _interp_matrix(x, A, xnew):
    x = np.asarray(x, float)
    A = np.asarray(A, float)
    out = np.empty((len(xnew), A.shape[1], A.shape[2]), float)
    for i in range(A.shape[1]):
        for j in range(A.shape[2]):
            out[:, i, j] = PchipInterpolator(x, A[:, i, j])(xnew)
    return out


def resample_operator(r, u, K, G, S, M, n):
    r = np.asarray(r, float)
    order = np.argsort(r)
    r = r[order]
    u = np.asarray(u, float)[order]
    if np.any(np.diff(r) <= 0):
        raise ValueError("r must be unique")
    rg = np.linspace(r[0], r[-1], int(n))
    ug = PchipInterpolator(r, u)(rg)
    mats = [_interp_matrix(r, np.asarray(A)[order], rg) for A in (K, G, S, M)]
    return rg, ug, *mats


def assemble_dirichlet_fem(r, K, G, S, M):
    """Assemble symmetric finite-element quadratic forms with zero endpoints."""
    r = np.asarray(r, float)
    K = np.asarray(K, float)
    G = np.asarray(G, float)
    S = np.asarray(S, float)
    M = np.asarray(M, float)
    n, f = K.shape[:2]
    if n < 5 or any(A.shape != (n, f, f) for A in (G, S, M)):
        raise ValueError("invalid operator shapes")
    if np.any(np.diff(r) <= 0):
        raise ValueError("r not increasing")

    nd = n * f
    A = np.zeros((nd, nd), float)
    B = np.zeros((nd, nd), float)

    # Trapezoid node weights for K and M.
    w = np.empty(n, float)
    h = np.diff(r)
    w[0] = h[0] / 2
    w[-1] = h[-1] / 2
    w[1:-1] = (h[:-1] + h[1:]) / 2

    for i in range(n):
        sl = slice(i * f, (i + 1) * f)
        Ks = (K[i] + K[i].T) / 2
        Ms = (M[i] + M[i].T) / 2
        B[sl, sl] += w[i] * Ks
        A[sl, sl] += w[i] * Ms

    for i in range(n - 1):
        hh = r[i + 1] - r[i]
        Gi = (G[i] + G[i + 1]) / 2
        Gi = (Gi + Gi.T) / 2
        Si = (S[i] + S[i + 1]) / 2
        Si = (Si - Si.T) / 2
        a = slice(i * f, (i + 1) * f)
        b = slice((i + 1) * f, (i + 2) * f)

        # integral Y'^T G Y' dr
        A[a, a] += Gi / hh
        A[b, b] += Gi / hh
        A[a, b] += -Gi / hh
        A[b, a] += -Gi / hh

        # - integral Y'^T S Y dr = Y_i^T S_mid Y_{i+1}
        A[a, b] += 0.5 * Si
        A[b, a] += 0.5 * Si.T

    # Homogeneous Dirichlet endpoints.
    keep_nodes = np.arange(1, n - 1)
    keep = np.concatenate([np.arange(i * f, (i + 1) * f) for i in keep_nodes])
    return A[np.ix_(keep, keep)], B[np.ix_(keep, keep)], w, keep


def solve_box_spectrum(r, u, K, G, S, M, *, modes=16):
    A, B, w, keep = assemble_dirichlet_fem(r, K, G, S, M)
    asymA = float(np.max(np.abs(A - A.T)) / max(1.0, np.max(np.abs(A))))
    asymB = float(np.max(np.abs(B - B.T)) / max(1.0, np.max(np.abs(B))))
    if asymA > 1e-10 or asymB > 1e-10:
        raise ValueError("assembled pencil is not symmetric")
    be = np.linalg.eigvalsh(B)
    if float(be.min()) <= 0:
        raise ValueError("finite-box kinetic mass matrix is not positive")

    m = min(int(modes), len(A) - 1)
    vals, vecs = eigh(A, B, subset_by_index=(0, m - 1), check_finite=True)
    neg = int(np.sum(vals <= 0))

    n = len(r)
    f = K.shape[1]
    psi = np.zeros((m, n, f), float)
    for q in range(m):
        full = np.zeros(n * f, float)
        full[keep] = vecs[:, q]
        psi[q] = full.reshape(n, f)

    density = np.einsum("mni,nij,mnj->mn", psi, K, psi)
    norms = np.sum(density * w[None, :], axis=1)
    if np.any(norms <= 0):
        raise ValueError("non-positive modal K norm")
    psi /= np.sqrt(norms)[:, None, None]
    density = np.einsum("mni,nij,mnj->mn", psi, K, psi)
    weights = density

    # Dimensionless IPR on s=(r-r0)/(r1-r0): integral p(s)^2 ds.
    span = r[-1] - r[0]
    ipr = np.trapezoid((weights * span) ** 2, r, axis=1) / span
    omega = np.where(vals > 0, np.sqrt(vals), np.nan)
    return BoxSpectrum(
        r=np.asarray(r),
        u=np.asarray(u),
        omega2=vals,
        omega=omega,
        psi=psi,
        weights=weights,
        ipr=ipr,
        negative_omega2_count=neg,
        symmetry_error_A=asymA,
        symmetry_error_B=asymB,
    )


def robust_pairwise_inversions(weights, *, rel_margin=2e-2):
    """Pairs whose normalized local weights robustly reverse ordering."""
    W = np.asarray(weights, float)
    if W.ndim != 2:
        raise ValueError("weights must be (modes,nodes)")
    total = np.sum(W, axis=0)
    P = W / np.maximum(total[None, :], 1e-300)
    pairs = []
    for a in range(len(P)):
        for b in range(a + 1, len(P)):
            d = P[a] - P[b]
            if np.max(d) > rel_margin and np.min(d) < -rel_margin:
                pairs.append((a, b))
    return pairs


def spectral_density_lambda(omega2, weights, *, n_lambda=600, eta_fraction=2e-3):
    """Lorentz-broadened local spectral density in lambda=omega^2."""
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
