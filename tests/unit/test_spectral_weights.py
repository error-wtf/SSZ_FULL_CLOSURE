import numpy as np
import pytest

from ssz_p5.qnm.spectral_weights import (
    has_pairwise_reordering,
    kinetic_density_weights,
    qnm_pole_weight,
    ranking_by_radius,
)


def test_bound_mode_kinetic_density_normalizes():
    r = np.linspace(0.0, 1.0, 101)
    k = np.repeat(np.eye(2)[None, :, :], len(r), axis=0)
    modes = np.zeros((2, len(r), 2))
    modes[0, :, 0] = 1.0
    modes[1, :, 1] = 1.0 + r
    _, weights = kinetic_density_weights(r, k, modes)
    np.testing.assert_allclose(np.trapezoid(weights, r, axis=1), 1.0, atol=1e-12)


def test_bound_mode_rejects_ghost_metric():
    r = np.linspace(0.0, 1.0, 11)
    k = np.repeat(np.diag([1.0, -1.0])[None, :, :], len(r), axis=0)
    modes = np.ones((1, len(r), 2))
    with pytest.raises(ValueError, match="positive definite"):
        kinetic_density_weights(r, k, modes)


def test_qnm_residue_is_left_right_rescaling_invariant():
    left = np.array([1.0 + 0.2j, 0.4 - 0.1j])
    right = np.array([0.7 - 0.3j, 1.2 + 0.1j])
    deriv = np.array([[2.0, 0.1j], [-0.1j, 1.5]], complex)
    probe = np.array([1.0, 0.5j])
    w0 = qnm_pole_weight(left, right, deriv, probe)
    alpha = 2.3 - 0.7j
    w1 = qnm_pole_weight(left / np.conj(alpha), alpha * right, deriv, probe)
    np.testing.assert_allclose(w0, w1, rtol=1e-12, atol=1e-12)


def test_ranking_and_pairwise_reordering():
    weights = np.array([[0.9, 0.6, 0.2], [0.1, 0.4, 0.8]])
    ranking = ranking_by_radius(weights)
    assert ranking[0].tolist() == [0, 1]
    assert ranking[-1].tolist() == [1, 0]
    assert has_pairwise_reordering(weights)
