"""Sagnac reference transport - N-segment chain (Route 3).

Discretize the ring into N arcs (uniform OR randomized partition).  In arc i
the signal travels at c relative to the medium, the receiver moves at v.
Per arc: dt_i = (L_i/N) / (c - s v).  Total: sum over arcs.

In units c = L = 1: t_s(N) = sum_i w_i / (1 - s*beta) = 1/(1 - s*beta)
independent of the partition weights w_i (sum w_i = 1) — the exactness of
this partition-independence IS the property under test (G106): the discrete
counted-segment methodology must reproduce the analytic transport observable
for arbitrary partitions and converge trivially exactly in N.
"""
from __future__ import annotations

import numpy as np


def segment_chain(beta: float, n_segments: int, s: int = +1,
                  seed: int | None = None) -> float:
    """Chain transport time t_s for a partition of the unit ring.

    s = +1 co-rotating, -1 counter-rotating.  seed=None -> uniform partition;
    seed=int -> randomized Dirichlet(1,1,..) partition (sum of weights = 1).
    """
    if n_segments < 1:
        raise ValueError("n_segments must be >= 1")
    if not -1.0 < beta < 1.0:
        raise ValueError(f"beta must satisfy |beta| < 1, got {beta}")
    if s not in (+1, -1):
        raise ValueError("s must be +1 or -1")
    if seed is None:
        weights = np.full(n_segments, 1.0 / n_segments)
    else:
        rng = np.random.default_rng(seed)
        w = rng.dirichlet(np.ones(n_segments))
        weights = w
    # per-arc time in units where arc length = w_i, speed = c - s*v = 1 - s*beta
    speed = 1.0 - s * beta
    return float(np.sum(weights / speed))


def chain_convergence(beta: float, s: int = +1, seed: int | None = None,
                      n_max: int = 4096) -> list[tuple[int, float]]:
    """t_s(N) for N in a geometric ladder up to n_max."""
    out = []
    n = 8
    while n <= n_max:
        out.append((n, segment_chain(beta, n, s, seed)))
        n *= 2
    return out
