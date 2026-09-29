import numpy as np

from ssz_p5.benchmarks.li_cst_aah import critical_site
from ssz_p5.benchmarks.li_cst_aah import solve as solve_li
from ssz_p5.benchmarks.weisz_fk import solve as solve_fk
from ssz_p5.qnm.principal_tracking import (
    continuity_metrics,
    track_by_overlap,
    whitened_principal_modes,
)
from ssz_p5.qnm.residue import pole_residue, spectral_density_from_green


def test_li_published_critical_site():
    assert critical_site(2584, 1.5, 1.0, 1.0) == 1937


def test_li_small_problem_is_normalized():
    r = solve_li(n=144, lam=1.5, sigma=1.0)
    np.testing.assert_allclose(np.sum(np.abs(r.eigenvectors) ** 2, axis=0), 1.0, atol=1e-12)


def test_weisz_optical_residue_completeness():
    r = solve_fk(0.02, 48, 49)
    assert abs(r.optical_weights.sum() - 1.0) < 1e-12
    assert np.max(r.optical_weights) > 0.1


def test_overlap_tracking_recovers_crossing_identity():
    # Naive eigenvalue order swaps the physical basis vectors at the crossing.
    vals = np.array([[0.0, 2.0], [0.9, 1.1], [0.8, 1.2]])
    vecs = np.array([
        np.eye(2),
        np.eye(2),
        np.array([[0.0, 1.0], [1.0, 0.0]]),
    ])
    _, tracked, perms = track_by_overlap(vals, vecs)
    assert np.array_equal(perms[-1], [1, 0])
    assert continuity_metrics(tracked)["min_adjacent_overlap"] == 1.0


def test_whitened_principal_modes():
    K = np.repeat(np.eye(2)[None], 3, axis=0)
    G = np.array([np.diag([1.0, 2.0]), np.diag([1.2, 1.8]), np.diag([1.4, 1.6])])
    vals, _ = whitened_principal_modes(K, G)
    np.testing.assert_allclose(vals[:, 0], [1.0, 1.2, 1.4])


def test_qnm_residue_helper_and_green_density():
    right = np.array([1.0 + 0j])
    left = np.array([1.0 + 0j])
    d = np.array([[2.0 + 0j]])
    assert pole_residue(right, left, d) == 0.5
    g = np.array([[[1.0 - 2.0j]], [[2.0 - 4.0j]]])
    rho = spectral_density_from_green(g)
    np.testing.assert_allclose(rho, [2.0 / np.pi, 4.0 / np.pi])
