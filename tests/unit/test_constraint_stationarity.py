"""Test the maps against Euler equations of the unreduced action."""

import numpy as np
import pandas as pd
import pytest

from ssz_p5.config import SLOT_NAMES
from ssz_p5.numerics import module


@pytest.mark.parametrize("L", [6, 42, 1000])
def test_h1_and_vector_constraints_stationary(L):
    n = 17
    r = np.linspace(1, 2, n)
    d = pd.DataFrame({s: np.full(n, 0.03 * (i + 1)) for i, s in enumerate(SLOT_NAMES)})
    d["x"] = r
    d["a1"] = 0.1 * r**2
    d["a3"] = -r
    d["a4"] = 0.5
    d["b1"] = 0.2
    d["v10"] = -0.4
    d["v11"] = 0.03
    m = module("ssz_hybrid_full_constraint_maps_JET9D8(1).py").maps(d, L)
    for order in (0, 1):
        h = m[f"H1{order}"]
        H = m[f"H2{order}"]
        H1 = m[f"H1dot{order}"]
        A1 = m[f"dA1dot{order}"]
        J = d.b4.to_numpy()[:, None] * H + L * d.b5.to_numpy()[:, None] * h
        J[:, 1] += d["b3" if order == 0 else "b2"]
        residual = 2 * L * d.b1.to_numpy()[:, None] * H1 + L * d.v11.to_numpy()[:, None] * A1 + J
        np.testing.assert_allclose(residual, 0, atol=1e-11)
        source = -0.5 * L * d.v6.to_numpy()[:, None] * h
        if order == 0:
            source[:, 2] += 2 * d.v1
        residual = (
            2 * L * d.v10.to_numpy()[:, None] * A1 + L * d.v11.to_numpy()[:, None] * H1 + source
        )
        np.testing.assert_allclose(residual, 0, atol=1e-11)
