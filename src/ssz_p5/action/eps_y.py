from __future__ import annotations

import numpy as np

from ..types import Coefficients41, P5Background


def kappa(bg: P5Background) -> np.ndarray:
    return np.asarray(bg.h, float) * np.asarray(bg.phi_r, float) ** 2


def za(bg: P5Background, epsilon_y: float) -> np.ndarray:
    return 1.0 - 2.0 * float(epsilon_y) * kappa(bg)


def apply_epsilon_y_deformation(c: Coefficients41, epsilon_y: float) -> Coefficients41:
    """Apply the frozen production deformation epsilon_Y * Y on A0'=0.

    For the production zero-vector branch the quadratic action changes only the
    Maxwell-vector kinetic/gradient slots v1 and v10 by Z_A=1-2 epsilon_Y kappa.
    """
    if (
        not np.all(np.isfinite(c.background.A0prime))
        or np.max(np.abs(c.background.A0prime)) > 1e-12
    ):
        raise ValueError("production epsilon_Y deformation requires A0prime=0")
    z = za(c.background, epsilon_y)
    if not np.all(np.isfinite(z)) or np.min(z) <= 0:
        raise ValueError("epsilon_Y deformation violates Z_A>0")
    slots = {k: np.array(v, float, copy=True) for k, v in c.slots.items()}
    for key in ("v1", "v10"):
        if key not in slots:
            raise KeyError(key)
        slots[key] *= z
    return Coefficients41(c.background, slots, c.provenance_id + f":epsY={epsilon_y:g}")
