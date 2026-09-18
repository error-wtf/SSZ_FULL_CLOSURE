import numpy as np

from ssz_p5.production.inner_targets import left_endpoint_jets, partitioned_target


def test_central_endpoint_jet_uses_only_central_side():
    r = np.linspace(1.0, 2.0, 101)
    y = 3 + 2 * (r - 1.5)
    y[r < 1.5] = 1e20
    jets = left_endpoint_jets(r, y, 1.5, 1)
    np.testing.assert_allclose(jets, [3, 2], atol=1e-11)


def test_target_preserves_flat_partition_endpoint_jet():
    r = np.linspace(1.5, 1.4, 101)
    S = np.ones(101)
    S[-10:] = 0
    q = partitioned_target(r, S, [3, 2], 1.5)
    np.testing.assert_allclose(q[:10], 3 + 2 * (r[:10] - 1.5))
    np.testing.assert_array_equal(q[-10:], 0)
