"""Analytic controls for the necessary background identity, not closure certificates."""

import numpy as np
import pytest

from ssz_p5.action.onshell import identity_terms, luminal_e00_minus_e22


@pytest.mark.parametrize("mass,cosmological", [(0.0, 0.0), (0.4, 0.0), (0.0, 0.01)])
def test_exact_vacuum_controls(mass, cosmological):
    r = np.linspace(2.0, 4.0, 21)
    f = 1 - 2 * mass / r - cosmological * r**2
    fp = 2 * mass / r**2 - 2 * cosmological * r
    fpp = -4 * mass / r**3 - 2 * cosmological
    lhs, rhs = identity_terms(r, f, f, 1.0, f / 2, fp / 2, fp, fpp)
    np.testing.assert_allclose(lhs, rhs, atol=1e-14, rtol=0)
    residual = luminal_e00_minus_e22(r, f, f, 1.0, fp, 0.0, fp, fpp)
    np.testing.assert_allclose(residual, 0.0, atol=1e-14, rtol=0)


def test_reissner_nordstrom_requires_electric_term():
    r = np.linspace(2.0, 4.0, 21)
    mass, charge = 0.4, 0.7
    f = 1 - 2 * mass / r + charge**2 / (2 * r**2)
    fp = 2 * mass / r**2 - charge**2 / r**3
    fpp = -4 * mass / r**3 + 3 * charge**2 / r**4
    lhs, rhs = identity_terms(r, f, f, 1.0, f / 2, fp / 2, fp, fpp,
                              electric_prime=charge / r**2, v8=1 / (2 * f))
    np.testing.assert_allclose(lhs, rhs, atol=1e-14, rtol=0)
    lhs, rhs = identity_terms(r, f, f, 1.0, f / 2, fp / 2, fp, fpp)
    np.testing.assert_allclose(lhs - rhs, -f * charge**2 / r**3, atol=1e-14, rtol=0)


def test_independent_background_subtraction_with_varying_G4():
    # Off-shell analytic functions deliberately exercise every derivative term.
    r = np.linspace(0.5, 2.0, 15)
    f, h, H = 1 + r**2, 1 + r, 2 + r**3
    fp, fpp, hp, Hp = 2 * r, 2.0, 1.0, 3 * r**2
    a4 = np.sqrt(f * h) * H / 2
    ap = np.sqrt(f * h) / 2 * (Hp + H / 2 * (fp / f + hp / h))
    lhs, rhs = identity_terms(r, f, h, H, a4, ap, fp, fpp)
    delta = luminal_e00_minus_e22(r, f, h, H, hp, Hp, fp, fpp)
    np.testing.assert_allclose(lhs - rhs, -r * f * np.sqrt(f / h) * delta,
                               atol=1e-12, rtol=1e-13)


def test_no_division_at_light_ring():
    r, f, h, H, fp, fpp = 2.0, 1.0, 0.7, 1.0, 1.0, 0.1
    lhs, rhs = identity_terms(r, f, h, H, np.sqrt(f * h) / 2, 0.3, fp, fpp)
    assert lhs == 0.0
    assert np.isfinite(rhs)
