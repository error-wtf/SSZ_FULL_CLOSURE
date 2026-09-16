import numpy as np

from ssz_p5.jets.jet9d8 import derivative


def test_polynomial_exact_interior():
    x = np.linspace(0.2, 2.0, 41)
    y = x**5 - 2 * x**3 + 7 * x
    ref = 5 * x**4 - 6 * x**2 + 7
    got = derivative(x, y, window=9, degree=8)
    assert np.max(np.abs(got - ref)) < 1e-8
