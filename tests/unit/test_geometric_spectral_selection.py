import numpy as np

from ssz_p5.qnm.principal_tracking import analyze_principal_tracking
from ssz_p5.qnm.spectral_weights import (
    has_pairwise_reordering,
    normalized_kinetic_weights,
)


def test_principal_tracking_keeps_continuous_crossing_identity():
    t = np.linspace(0.0, 1.0, 101)
    K = np.repeat(np.eye(3)[None, :, :], len(t), axis=0)
    G = np.empty_like(K)
    # Smooth rotation of two eigenvectors while eigenvalues cross.
    for i, x in enumerate(t):
        a = np.pi * x / 3
        R = np.array(
            [[np.cos(a), -np.sin(a), 0.0],
             [np.sin(a), np.cos(a), 0.0],
             [0.0, 0.0, 1.0]]
        )
        vals = np.diag([0.5 + x, 1.5 - x, 3.0])
        G[i] = R @ vals @ R.T
    out = analyze_principal_tracking(K, G)
    assert out["min_adjacent_overlap"] > 0.99


def test_bound_mode_weights_detect_reordering():
    r = np.linspace(0.0, 1.0, 101)
    K = np.repeat(np.eye(2)[None, :, :], len(r), axis=0)
    psi = np.zeros((2, len(r), 2))
    psi[0, :, 0] = 1.0 - 0.8 * r
    psi[1, :, 1] = 0.2 + 0.8 * r
    w = normalized_kinetic_weights(r, K, psi)
    assert has_pairwise_reordering(w)


def test_bound_mode_weights_reject_ghost_metric():
    r = np.linspace(0.0, 1.0, 11)
    K = np.repeat(np.eye(2)[None, :, :], len(r), axis=0)
    K[:, 1, 1] = -1.0
    psi = np.ones((1, len(r), 2))
    try:
        normalized_kinetic_weights(r, K, psi)
    except ValueError as exc:
        assert "positive definite" in str(exc)
    else:
        raise AssertionError("ghost metric was accepted")
