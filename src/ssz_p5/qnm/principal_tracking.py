"""Overlap tracking for the local canonical principal K/G problem."""
from __future__ import annotations

import itertools

import numpy as np


def whitened_principal_modes(K: np.ndarray, G: np.ndarray):
    K = np.asarray(K, float)
    G = np.asarray(G, float)
    if K.shape != G.shape or K.ndim != 3 or K.shape[1] != K.shape[2]:
        raise ValueError("K and G must have shape (N,F,F)")
    values = []
    vectors = []
    for k, g in zip(K, G, strict=True):
        ks = (k + k.T) / 2.0
        gs = (g + g.T) / 2.0
        ew, U = np.linalg.eigh(ks)
        if np.min(ew) <= 0:
            raise ValueError("K is not positive definite")
        invsqrt = U @ np.diag(1.0 / np.sqrt(ew)) @ U.T
        c = invsqrt @ gs @ invsqrt
        ev, V = np.linalg.eigh((c + c.T) / 2.0)
        values.append(ev)
        vectors.append(V)
    return np.asarray(values), np.asarray(vectors)


def track_by_overlap(values: np.ndarray, vectors: np.ndarray):
    """Track branches by maximizing adjacent absolute eigenvector overlap.

    Returns tracked values/vectors and the per-step permutations relative to
    naive ascending eigenvalue order.
    """
    values = np.asarray(values, float)
    vectors = np.asarray(vectors, float)
    n, f = values.shape
    if vectors.shape != (n, f, f):
        raise ValueError("vector shape mismatch")
    tv = values.copy()
    te = vectors.copy()
    perms = np.tile(np.arange(f), (n, 1))
    for i in range(1, n):
        overlap = np.abs(te[i - 1].T @ vectors[i])
        best_score = -1.0
        best = None
        for p in itertools.permutations(range(f)):
            score = sum(overlap[j, p[j]] for j in range(f))
            if score > best_score:
                best_score = score
                best = np.asarray(p, int)
        tv[i] = values[i, best]
        te[i] = vectors[i][:, best]
        perms[i] = best
        # Fix arbitrary real signs so continuity metrics are meaningful.
        signs = np.sign(np.sum(te[i - 1] * te[i], axis=0))
        signs[signs == 0] = 1.0
        te[i] *= signs
    return tv, te, perms


def continuity_metrics(vectors: np.ndarray) -> dict:
    vectors = np.asarray(vectors, float)
    overlaps = np.abs(np.einsum("nij,nij->nj", vectors[:-1], vectors[1:]))
    return {
        "min_adjacent_overlap": float(np.min(overlaps)),
        "mean_adjacent_overlap": float(np.mean(overlaps)),
        "below_0_99": int(np.sum(overlaps < 0.99)),
        "below_0_999": int(np.sum(overlaps < 0.999)),
    }
