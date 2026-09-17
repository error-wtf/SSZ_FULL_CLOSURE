import numpy as np

from ssz_p5.jets.jet9d8 import derivative


def test_polynomial_exact_interior():
    x = np.linspace(0.2, 2.0, 41)
    y = x**5 - 2 * x**3 + 7 * x
    ref = 5 * x**4 - 6 * x**2 + 7
    got = derivative(x, y, window=9, degree=8)
    assert np.max(np.abs(got - ref)) < 1e-8


def test_profile_service_supports_declared_convergence_stencils():
    from ssz_p5.jets.jet9d8 import profile_derivative

    x = np.linspace(0.2, 2.0, 45) ** 1.2
    y = x**5 - 2 * x**3 + 7 * x
    for window, degree in [(7, 6), (9, 8), (11, 8), (13, 8)]:
        for order in (1, 2):
            actual = profile_derivative(x, y, order, window, degree)
            expected = derivative(x, y, order, window=window, degree=degree)
            np.testing.assert_allclose(actual, expected, rtol=2e-8, atol=2e-8)
