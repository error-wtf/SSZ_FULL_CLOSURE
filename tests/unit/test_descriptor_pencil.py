import numpy as np
import pandas as pd
from scipy import sparse

from ssz_p5.qnm.descriptor_pencil import (
    local_poly_differentiation_matrix,
    semidiscrete_polynomial,
    quadratic_descriptor_linearization,
    pencil_structure,
)


def test_local_poly_matrix_differentiates_degree_eight_polynomial():
    x = np.linspace(0.3, 1.7, 19)
    D = local_poly_differentiation_matrix(x, 2)
    y = x**8 - 0.4*x**5 + 2*x**2
    want = 56*x**6 - 8*x**3 + 4
    got = D @ y
    assert np.max(np.abs(got - want)) < 2e-7


def test_quadratic_linearization_preserves_singular_mass():
    A0 = sparse.csr_matrix([[2.0, 0.0], [0.0, 3.0]])
    A1 = sparse.csr_matrix((2, 2))
    A2 = sparse.csr_matrix([[1.0, 0.0], [0.0, 0.0]])
    L, R = quadratic_descriptor_linearization(A0, A1, A2)
    assert L.shape == R.shape == (4, 4)
    assert np.linalg.matrix_rank(R.toarray()) < 4


def test_semidiscrete_descriptor_keeps_dae_singular_kinetic_rows():
    n = 17
    x = np.linspace(1.0, 2.0, n)
    P = {
        (2, 0): np.tile(np.diag([1.0, 0.0])[None, :, :], (n, 1, 1)),
        (0, 2): np.tile(np.eye(2)[None, :, :], (n, 1, 1)),
        (0, 0): np.tile(np.eye(2)[None, :, :], (n, 1, 1)),
    }
    A0, A1, A2, xs = semidiscrete_polynomial(P, x)
    s = pencil_structure(A0, A1, A2)
    assert len(xs) == n
    assert s["dimension"] == 2*n
    assert s["rank_A2"] == n
    assert s["singular_kinetic_expected"]
