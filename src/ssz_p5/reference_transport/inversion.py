"""Sagnac reference transport - branch-selected inversion (dt -> v).

Inverts dt = 2 beta / (1 - beta^2) for beta given dt (units c = L = 1).

Quadratic: dt * beta^2 + 2 beta - dt = 0
    beta = (-1 +/- sqrt(1 + dt^2)) / dt

Physical branch: |beta| < 1 and sign(beta) = sign(dt).
As |dt| -> inf, beta -> sign(dt) (|v| -> c).  No real solution only at dt=0
ambiguity (beta=0 exactly), which is returned as 0.0.
"""
from __future__ import annotations

import math


def invert_delta_t(dt: float) -> float:
    """Recover beta = v/c from the measured dt (units c = L = 1)."""
    if dt == 0.0:
        return 0.0
    disc = math.sqrt(1.0 + dt * dt)
    r1 = (-1.0 + disc) / dt
    r2 = (-1.0 - disc) / dt
    # exactly one root lies in (-1, 1) for dt != 0 (the other is outside by
    # Vieta: r1*r2 = -1).  Select it.
    for r in (r1, r2):
        if abs(r) < 1.0:
            return r
    raise ValueError(f"no physical root in (-1,1) for dt={dt!r}")
