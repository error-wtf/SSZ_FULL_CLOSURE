"""Accepted identities, shared by validation and negative regressions."""

import numpy as np

from ..jets.jet9d8 import derivative
from ..policy import numerical_policy


def a5_holonomic(r, a1, a2, a0prime, v4, v5):
    return (
        derivative(r, a2)
        - derivative(r, a1, 2)
        - derivative(r, a0prime * v4 / 2)
        + a0prime * v5 / 2
    )


def require_v12(slots, h):
    error = np.asarray(slots["v12"]) + np.asarray(slots["v6"]) / (2 * np.asarray(h))
    if (
        not np.all(np.isfinite(error))
        or np.max(np.abs(error)) > numerical_policy()["branch_identity_abs"]
    ):
        raise ValueError("Wrong v12 branch")
