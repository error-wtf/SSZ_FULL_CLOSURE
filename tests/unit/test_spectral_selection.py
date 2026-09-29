import numpy as np
import pytest

from ssz_p5.qnm.spectral_selection import (
    evaluate_spectral_selection,
    local_kinetic_weights,
)


def _base():
    r = np.linspace(0.0, 1.0, 101)
    K = np.repeat(np.eye(2)[None, :, :], len(r), axis=0)
    return r, K


def test_k_normalization():
    r, K = _base()
    psi = np.zeros((2, len(r), 2))
    psi[0, :, 0] = 1.0
    psi[1, :, 1] = np.sqrt(3.0) * r
    normalized, weights = local_kinetic_weights(r, K, psi)
    density = np.einsum("mni,nij,mnj->mn", normalized, K, normalized)
    np.testing.assert_allclose(np.trapezoid(density, r, axis=1), 1.0, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(np.trapezoid(weights, r, axis=1), 1.0, rtol=1e-12, atol=1e-12)


def test_detects_radial_weight_reordering():
    r, K = _base()
    psi = np.zeros((2, len(r), 2))
    psi[0, :, 0] = 1.0 - 0.8 * r
    psi[1, :, 1] = 0.2 + 0.8 * r
    result = evaluate_spectral_selection(r, K, psi)
    assert result.status == "SPECTRAL_SELECTION_PASS"
    assert result.reordered


def test_reports_null_when_ranking_is_fixed():
    r, K = _base()
    psi = np.zeros((2, len(r), 2))
    psi[0, :, 0] = 2.0
    psi[1, :, 1] = 1.0
    result = evaluate_spectral_selection(r, K, psi)
    assert result.status == "SPECTRAL_SELECTION_NULL"
    assert not result.reordered


def test_rejects_nonpositive_kinetic_metric():
    r, K = _base()
    K[:, 1, 1] = -1.0
    psi = np.ones((1, len(r), 2))
    with pytest.raises(ValueError, match="positive definite"):
        evaluate_spectral_selection(r, K, psi)
