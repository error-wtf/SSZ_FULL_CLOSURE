"""Sagnac reference transport - PDE first-passage route (Route 2).

Advection equation on the ring (units c = L = 1), direction s in {+1,-1}:

    d/dt psi + s * (1 - s*beta) * d/dx psi = 0,    x in [0, 1), periodic.

Emitter/receiver sit at x = 0.  The initial condition is a narrow Gaussian
packet centered at x = 1/2 (so the first arrival at x=0 is unambiguous —
the packet must travel the full half-ring in EITHER direction... no: for a
fair t_s measurement the packet starts at x=0+eps and wraps the whole ring).
Here: initial packet at x = 1 - w (just "behind" the receiver), so the
co-rotating packet must wrap almost the entire ring to reach the receiver.

Arrival = first time the receiver amplitude exceeds half of the packet's
initial peak (declared BEFORE the run; sub-sample refinement by linear
interpolation between the two bracketing time steps).

Numerics: method-of-lines with first-order upwind flux, CFL = 0.4.
Grid convergence to the closed form is gate G107_Sag (observed order >= 1).
"""
from __future__ import annotations

import numpy as np


def first_passage_time(beta: float, s: int = +1, n_grid: int = 4096,
                       cfl: float = 0.4, t_max: float = 8.0,
                       threshold_frac: float = 0.05) -> float:
    if s not in (+1, -1):
        raise ValueError("s must be +1 or -1")
    if not -1.0 < beta < 1.0:
        raise ValueError(f"beta must satisfy |beta| < 1, got {beta}")

    dx = 1.0 / n_grid
    speed = s * (1.0 - s * beta)          # signed phase speed
    dt_step = cfl * dx / abs(speed)
    x = (np.arange(n_grid) + 0.5) * dx    # cell centers
    w = 4.0 * dx
    # initial packet just behind the receiver (x slightly below 1), wrapping
    psi = np.exp(-0.5 * ((x - 4.0 * dx) / w) ** 2)
    peak0 = psi.max()

    t = 0.0
    prev_amp = psi[0]                      # receiver cell contains x=0
    while t < t_max:
        # upwind flux (periodic): d psi/dt = -speed * d psi/dx
        dp = np.empty_like(psi)
        if speed > 0:
            dp[1:] = psi[1:] - psi[:-1]
            dp[0] = psi[0] - psi[-1]
            psi = psi - speed * dt_step / dx * dp
        else:
            dp[:-1] = psi[1:] - psi[:-1]
            dp[-1] = psi[0] - psi[-1]
            psi = psi - speed * dt_step / dx * dp
        t += dt_step
        amp = psi[0]
        if amp >= threshold_frac * peak0 > prev_amp or (
                amp >= threshold_frac * peak0 and prev_amp < threshold_frac * peak0):
            # linear interpolation between prev_amp and amp
            frac = (threshold_frac * peak0 - prev_amp) / (amp - prev_amp)
            return t - (1.0 - frac) * dt_step
        prev_amp = amp
    raise RuntimeError(f"no arrival within t_max={t_max} (beta={beta}, s={s})")
