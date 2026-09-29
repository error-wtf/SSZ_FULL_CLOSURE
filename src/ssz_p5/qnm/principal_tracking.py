"""Principal-mode tracking for the frozen SSZ radial characteristic problem.

This module analyzes the local generalized eigenproblem

    G(r) x = c_r^2(r) K(r) x

with symmetric K-whitening.  It deliberately does *not* call these local
characteristic vectors global normal modes or QNMs.
"""

from __future__ import annotations

from itertools import permutations

import numpy as np


def whitened_characteristics(K: np.ndarray, G: np.ndarray):
    K = np.asarray(K, float)
    G = np.asarray(G, float)
    if K.shape != G.shape or K.ndim != 3 or K.shape[1:] != (3, 3):
        raise ValueError("K and G must have shape (N,3,3)")

    n = K.shape[0]
    values = np.empty((n, 3), float)
    whitened = np.empty((n, 3, 3), float)
    physical = np.empty((n, 3, 3), float)
    min_k = np.empty(n, float)

    for i in range(n):
        ks = (K[i] + K[i].T) / 2.0
        gs = (G[i] + G[i].T) / 2.0
        kw, kv = np.linalg.eigh(ks)
        if np.min(kw) <= 0:
            raise ValueError(f"nonpositive kinetic eigenvalue at radial row {i}")
        invsqrt = (kv * (1.0 / np.sqrt(kw))[None, :]) @ kv.T
        a = invsqrt @ gs @ invsqrt
        a = (a + a.T) / 2.0
        ew, ev = np.linalg.eigh(a)
        values[i] = ew
        whitened[i] = ev
        physical[i] = invsqrt @ ev
        min_k[i] = np.min(kw)

    return values, whitened, physical, min_k


def track_by_overlap(vectors: np.ndarray):
    """Track columns by maximum adjacent absolute overlap.

    Returns tracked vectors, the permutation chosen at each radial step, and
    per-branch adjacent absolute overlaps after sign alignment.
    """
    v = np.asarray(vectors, float)
    if v.ndim != 3 or v.shape[1:] != (3, 3):
        raise ValueError("vectors must have shape (N,3,3)")

    out = np.empty_like(v)
    out[0] = v[0]
    chosen = np.empty((len(v), 3), int)
    chosen[0] = np.arange(3)
    overlaps = np.ones((len(v), 3), float)

    perms = tuple(permutations(range(3)))
    for i in range(1, len(v)):
        raw = np.abs(out[i - 1].T @ v[i])
        best = max(perms, key=lambda p: sum(raw[j, p[j]] for j in range(3)))
        cur = v[i][:, best].copy()
        dots = np.sum(out[i - 1] * cur, axis=0)
        signs = np.where(dots < 0, -1.0, 1.0)
        cur *= signs[None, :]
        out[i] = cur
        chosen[i] = np.asarray(best)
        overlaps[i] = np.abs(np.sum(out[i - 1] * out[i], axis=0))

    return out, chosen, overlaps


def physical_component_fractions(physical_vectors: np.ndarray):
    x = np.asarray(physical_vectors, float)
    denom = np.sum(x * x, axis=1, keepdims=True)
    if np.any(denom <= 0):
        raise ValueError("zero physical vector")
    return (x * x) / denom


def analyze_principal_tracking(K: np.ndarray, G: np.ndarray):
    values, whitened, physical, min_k = whitened_characteristics(K, G)
    tracked_w, permutations_used, overlaps = track_by_overlap(whitened)

    # Apply the same tracked ordering/sign to physical vectors and eigenvalues.
    tracked_values = np.empty_like(values)
    tracked_physical = np.empty_like(physical)
    for i in range(len(values)):
        p = permutations_used[i]
        tracked_values[i] = values[i, p]
        cur = physical[i][:, p].copy()
        if i:
            # sign is irrelevant for component fractions; retain continuity by
            # aligning the whitened sign to the tracked vector.
            raw = whitened[i][:, p]
            signs = np.sign(np.sum(tracked_w[i] * raw, axis=0))
            signs[signs == 0] = 1.0
            cur *= signs[None, :]
        tracked_physical[i] = cur

    fractions = physical_component_fractions(tracked_physical)
    dominant = np.argmax(fractions, axis=1)  # (N, branch)
    dominant_changes = [
        int(np.count_nonzero(np.diff(dominant[:, j]) != 0)) for j in range(3)
    ]
    sort_relabels = int(
        np.count_nonzero(np.any(permutations_used[1:] != np.arange(3), axis=1))
    )

    endpoint_overlap = np.abs(tracked_w[0].T @ tracked_w[-1])
    return {
        "values": tracked_values,
        "min_k": min_k,
        "overlaps": overlaps,
        "permutations": permutations_used,
        "fractions": fractions,
        "dominant_components": dominant,
        "dominant_component_changes": dominant_changes,
        "naive_sort_relabel_steps": sort_relabels,
        "min_adjacent_overlap": float(np.min(overlaps[1:])),
        "min_adjacent_overlap_per_branch": [
            float(np.min(overlaps[1:, j])) for j in range(3)
        ],
        "endpoint_abs_overlap_matrix": endpoint_overlap,
    }
