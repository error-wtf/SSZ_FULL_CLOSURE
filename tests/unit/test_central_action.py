import numpy as np
import pandas as pd
import pytest

from ssz_p5.config import repo_root
from ssz_p5.production.central_action import (
    background_residuals,
    central_action_inputs,
    replay_central_action,
)
from ssz_p5.production.regional_coefficients import central_selected
from ssz_p5.reducer.kinetic_schur import kinetic_schur
from ssz_p5.reducer.zk_kinetic import zk_kinetic


def rn_action():
    r = np.linspace(3, 6, 401)
    charge = 0.4
    f = 1 - 2 / r + charge**2 / (2 * r**2)
    d = pd.DataFrame(dict(x=r, f=f, h=f, phiprime=np.ones(len(r)), A0prime=charge / r**2))
    for c in ("f2X", "f2Y", "f3", "f3X", "f4", "f4X", "f4XX", "tf4"):
        d[c] = 0.0
    d["f2F"] = 1.0
    d["f2"] = charge**2 / (2 * r**4)
    return d


def test_background_equations_reissner_nordstrom_control():
    out = background_residuals(rn_action())
    for name in ("E00", "E11", "JA_prime", "eq131_residual"):
        np.testing.assert_allclose(out[name][8:-8], 0, atol=1e-9)
    np.testing.assert_allclose(out["JA"], 0.4, atol=1e-14)


def test_background_checker_detects_action_change():
    d = rn_action()
    d["f2"] += 0.1
    out = background_residuals(d)
    np.testing.assert_allclose(out["E00"], -0.1 * d.x**2 * d.f, atol=1e-9)


def test_zero_electric_identity_remains_finite():
    d = rn_action()
    d["A0prime"] = 0.0
    assert np.isfinite(background_residuals(d)["eq131_residual"]).all()


def test_central_published_formula_matches_independent_schur():
    d = central_selected(repo_root()).sort_values("x").reset_index(drop=True)
    d = d[(d.u > 0.685) & (d.u < 0.695)].reset_index(drop=True)
    for L in (6, 42, 1000):
        np.testing.assert_allclose(
            zk_kinetic(d, L)[8:-8], kinetic_schur(d, L)["K"][8:-8], rtol=2e-7, atol=2e-7
        )
    d["a1"] = 0.1
    with pytest.raises(ValueError, match="Einstein"):
        zk_kinetic(d, 6)


def test_central_action_consumes_explicit_inputs():
    d = central_action_inputs(repo_root())
    out = replay_central_action(d)
    changed = d.copy()
    changed["f2F"] += 0.1
    other = replay_central_action(changed)
    np.testing.assert_allclose(other.v10 - out.v10, -0.05 * np.sqrt(d.f * d.h), atol=1e-12)
    ref = central_selected(repo_root())
    np.testing.assert_allclose(out.v5, ref.v5, rtol=2e-8, atol=2e-8)
    np.testing.assert_allclose(out.v6, ref.v6, rtol=1e-10, atol=1e-10)
    # d3 contains partial_phi(v6); the direct action selector must use the mixed
    # action jets rather than the total on-curve radial derivative.
    mask = (d.u >= 0.61) & (d.u < 0.71)
    err = abs(out.loc[mask, "d3"] - ref.loc[mask, "d3"]) / np.maximum(1, abs(ref.loc[mask, "d3"]))
    assert float(err.max()) < 3e-7
