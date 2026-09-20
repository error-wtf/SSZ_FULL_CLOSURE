import numpy as np
from ssz_p5.config import repo_root
from ssz_p5.production.hsvt_eps_y import build_region, EPSILON_Y


def test_direct_hsvt_carrier_is_finite_and_background_null_branch():
    d = build_region(repo_root(), "carrier")
    assert np.all(d.A0prime == 0)
    assert np.all(d.c3 == 0) and np.all(d.e3 == 0) and np.all(d.v5 == 0)
    assert np.min(d.Z_A) > 0
    np.testing.assert_allclose(d.Z_A, 1 - 2 * EPSILON_Y * d.kappa, rtol=0, atol=1e-14)
