import numpy as np

from ssz_p5.production.strong_field_transition import (
    U_HISTORICAL_SEAM,
    U_INNER_LIGHT_RING,
    U_LEFT,
    U_RIGHT,
    compact_transition_bump,
    first_zero_crossing,
    track_modes,
    transition_basis,
    transition_coordinate,
)


def test_physical_interval_contains_light_ring_but_not_a_physical_seam():
    assert U_LEFT < U_INNER_LIGHT_RING < U_RIGHT
    assert U_LEFT < U_HISTORICAL_SEAM < U_RIGHT
    assert not np.isclose(U_HISTORICAL_SEAM, U_INNER_LIGHT_RING)


def test_transition_coordinate_and_basis_are_flat_at_physical_endpoints():
    u = np.array([U_LEFT, (U_LEFT + U_RIGHT) / 2.0, U_RIGHT])
    s = transition_coordinate(u)
    assert np.allclose(s, [0.0, 0.5, 1.0])
    bump = compact_transition_bump(u)
    assert bump[0] == 0.0
    assert bump[-1] == 0.0
    assert bump[1] > 0.0
    basis = transition_basis(u, order=4)
    assert basis.shape == (3, 4)
    assert np.all(basis[0] == 0.0)
    assert np.all(basis[-1] == 0.0)


def test_mode_tracking_preserves_eigenvector_identity_through_order_exchange():
    # Diagonal matrices with exchanged ordered eigenvalues are a simple exact
    # test that the overlap tracker follows vectors instead of resorting labels.
    K = np.array([
        np.diag([1.0, 2.0, 3.0]),
        np.diag([2.1, 1.1, 3.0]),
        np.diag([2.2, 1.2, 3.0]),
    ])
    vals, vecs = track_modes(K)
    # tracked branch 0 stays e_x and therefore follows 1.0 -> 2.1 -> 2.2
    assert np.allclose(vals[:, 0], [1.0, 2.1, 2.2])
    assert np.allclose(np.abs(vecs[:, :, 0]), [[1, 0, 0]] * 3)


def test_first_zero_crossing_is_interpolated():
    u = np.array([0.70, 0.71, 0.72])
    y = np.array([1.0, -1.0, -2.0])
    assert np.isclose(first_zero_crossing(u, y), 0.705)
