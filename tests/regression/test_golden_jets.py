import importlib.util

import numpy as np

from ssz_p5.config import repo_root
from ssz_p5.jets.jet9d8 import derivative, profile_derivative


def test_shared_jet_matches_frozen_scalar_implementation():
    p = repo_root() / "archive/full_working_snapshot/ssz_p5_higher_jet_closure_2026-09-16.py"
    spec = importlib.util.spec_from_file_location("frozen_hj", p)
    ref = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ref)
    x = np.geomspace(0.5, 2.0, 51)
    y = np.sin(x)
    for order in (1, 2, 3):
        expected = ref.local_poly_deriv(x, y, order, 9, 8)
        np.testing.assert_array_equal(derivative(x, y, order), expected)
        np.testing.assert_allclose(profile_derivative(x, y, order), expected, atol=2e-7, rtol=2e-7)
