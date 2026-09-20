"""Radial semi-discrete descriptor pencil for the generalized-psi DAE.

This module deliberately preserves singular kinetic rows.  It converts a
profile differential Euler operator

    sum_{pt,pr} P[pt,pr](r) d_t^pt d_r^pr y = 0

into the quadratic semi-discrete pencil

    (A2 lambda^2 + A1 lambda + A0) y = 0

using a local-polynomial differentiation matrix on the chosen radial nodes.
No constraint Schur complement is taken.  Physical QNM boundary conditions are
an independent layer and are not supplied by this module.
"""
from __future__ import annotations

from math import comb, factorial
from typing import Dict, Iterable, Tuple

import numpy as np
from scipy import sparse

Operator = Dict[Tuple[int, int], np.ndarray]


def local_poly_differentiation_matrix(
    x: np.ndarray,
    order: int,
    window: int = 9,
    degree: int = 8,
) -> sparse.csr_matrix:
    """Return the linear counterpart of the local degree-8 jet derivative.

    The stencil selection mirrors JET9D8: a centered window when possible and
    one-sided windows at the radial boundaries.  Coordinates are locally scaled
    before solving the Vandermonde system to avoid gratuitous conditioning loss.
    """
    xx = np.asarray(x, float)
    n = len(xx)
    if n < window:
        raise ValueError(f"need at least {window} radial nodes")
    if degree >= window:
        raise ValueError("degree must be smaller than window")
    if order < 0 or order > degree:
        raise ValueError("invalid derivative order")
    if np.any(np.diff(xx) <= 0):
        raise ValueError("radial nodes must be strictly increasing")
    if order == 0:
        return sparse.identity(n, format="csr")

    rows: list[int] = []
    cols: list[int] = []
    vals: list[float] = []
    half = window // 2
    for i in range(n):
        lo = max(0, min(i - half, n - window))
        js = np.arange(lo, lo + window)
        dx = xx[js] - xx[i]
        scale = float(max(np.max(np.abs(dx)), np.finfo(float).tiny))
        z = dx / scale
        V = np.vstack([z**p for p in range(degree + 1)])
        rhs = np.zeros(degree + 1, float)
        rhs[order] = factorial(order) / scale**order
        w = np.linalg.solve(V, rhs)
        rows.extend([i] * window)
        cols.extend(js.tolist())
        vals.extend(w.tolist())
    return sparse.coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()


def semidiscrete_polynomial(
    P: Operator,
    x: np.ndarray,
    *,
    indices: Iterable[int] | None = None,
    window: int = 9,
    degree: int = 8,
) -> tuple[sparse.csr_matrix, sparse.csr_matrix, sparse.csr_matrix, np.ndarray]:
    """Build ``A0,A1,A2`` without eliminating descriptor constraints.

    Field-major ordering is used: all radial nodes of field 0, then field 1,
    and so on.  If ``indices`` is supplied the operator coefficients are sampled
    from the full profile while perturbation derivatives are represented on the
    selected radial nodes.
    """
    if not P:
        raise ValueError("empty differential operator")
    sample = next(iter(P.values()))
    nf = int(sample.shape[1])
    full_x = np.asarray(x, float)
    if indices is None:
        idx = np.arange(len(full_x), dtype=int)
    else:
        idx = np.asarray(list(indices), dtype=int)
    xs = full_x[idx]
    n = len(xs)
    if len(np.unique(idx)) != n or np.any(np.diff(idx) <= 0):
        raise ValueError("indices must be unique and strictly increasing")

    max_pr = max(pr for _, pr in P)
    D = {k: local_poly_differentiation_matrix(xs, k, window, degree) for k in range(max_pr + 1)}
    A = [sparse.csr_matrix((nf * n, nf * n), dtype=float) for _ in range(3)]

    for (pt, pr), block in P.items():
        if pt > 2:
            raise ValueError("only quadratic-in-time descriptors are supported")
        B = np.asarray(block[idx], float)
        Dr = D[pr]
        pieces = [[None for _ in range(nf)] for __ in range(nf)]
        any_nonzero = False
        for i in range(nf):
            for j in range(nf):
                a = B[:, i, j]
                if np.any(a):
                    pieces[i][j] = sparse.diags(a, format="csr") @ Dr
                    any_nonzero = True
                else:
                    pieces[i][j] = sparse.csr_matrix((n, n), dtype=float)
        if any_nonzero:
            A[pt] = A[pt] + sparse.bmat(pieces, format="csr")
    return A[0], A[1], A[2], xs


def quadratic_descriptor_linearization(
    A0: sparse.spmatrix,
    A1: sparse.spmatrix,
    A2: sparse.spmatrix,
) -> tuple[sparse.csr_matrix, sparse.csr_matrix]:
    """Return ``L z = lambda R z`` for ``A2 l^2+A1 l+A0``.

    ``A2`` may be singular; this is intentional for a DAE descriptor and leads
    to generalized infinite eigenvalues rather than an artificial inverse.
    """
    A0 = sparse.csr_matrix(A0)
    A1 = sparse.csr_matrix(A1)
    A2 = sparse.csr_matrix(A2)
    if A0.shape != A1.shape or A0.shape != A2.shape or A0.shape[0] != A0.shape[1]:
        raise ValueError("quadratic blocks must be same-size square matrices")
    n = A0.shape[0]
    I = sparse.identity(n, format="csr")
    Z = sparse.csr_matrix((n, n), dtype=float)
    left = sparse.bmat([[-A1, -A0], [I, Z]], format="csr")
    right = sparse.bmat([[A2, Z], [Z, I]], format="csr")
    return left, right


def pencil_structure(A0: sparse.spmatrix, A1: sparse.spmatrix, A2: sparse.spmatrix) -> dict:
    """Cheap structural diagnostics for a radial descriptor pencil."""
    A0 = sparse.csr_matrix(A0)
    A1 = sparse.csr_matrix(A1)
    A2 = sparse.csr_matrix(A2)
    n = A0.shape[0]
    # Dense rank is only used for deliberately small audit pencils.
    rank_a2 = int(np.linalg.matrix_rank(A2.toarray())) if n <= 512 else None
    return {
        "dimension": int(n),
        "nnz_A0": int(A0.nnz),
        "nnz_A1": int(A1.nnz),
        "nnz_A2": int(A2.nnz),
        "rank_A2": rank_a2,
        "singular_kinetic_expected": bool(rank_a2 is None or rank_a2 < n),
    }
